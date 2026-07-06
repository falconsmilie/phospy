"""Reference-manifest loading and validation."""

from __future__ import annotations

import json
import re
from datetime import date
from hashlib import sha256
from pathlib import Path
from typing import cast

from phospy.science.references.errors import ReferenceManifestError
from phospy.science.references.manifest import (
    REFERENCE_MANIFEST_SCHEMA_VERSION,
    RedistributionStatus,
    ReferenceFileManifest,
    ReferenceManifest,
)

_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_REQUIRED_MANIFEST_FIELDS = frozenset(
    {
        "reference_id",
        "display_name",
        "organism",
        "taxonomy_id",
        "protein_namespace",
        "reference_version",
        "source_name",
        "source_url",
        "source_version",
        "retrieved_at",
        "table_sha256",
        "license_name",
        "license_url",
        "redistribution_status",
        "redistribution_notes",
        "derived_from",
        "generated_by",
        "generated_at_utc",
        "manifest_schema_version",
        "files",
    }
)
_REQUIRED_FILE_FIELDS = frozenset(
    {
        "relative_path",
        "role",
        "format",
        "sha256",
        "row_count",
        "column_names",
    }
)


def load_reference_manifest(
    manifest_path: str | Path,
    *,
    bundle_root: str | Path | None = None,
    bundled: bool = False,
    require_redistribution_allowed: bool = False,
    require_all_files_listed: bool = False,
) -> ReferenceManifest:
    """Load, parse, and validate one reference manifest JSON file."""

    resolved_manifest_path = Path(manifest_path)
    try:
        payload = json.loads(resolved_manifest_path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ReferenceManifestError(
            f"reference manifest does not exist: {resolved_manifest_path}"
        ) from exc
    except json.JSONDecodeError as exc:
        raise ReferenceManifestError(
            f"reference manifest is not valid JSON: {resolved_manifest_path}"
        ) from exc
    if not isinstance(payload, dict):
        raise ReferenceManifestError(
            f"reference manifest must decode to an object: {resolved_manifest_path}"
        )
    root = (
        Path(bundle_root) if bundle_root is not None else resolved_manifest_path.parent
    )
    manifest = parse_reference_manifest_payload(
        cast(dict[str, object], payload),
        context=str(resolved_manifest_path),
    )
    return validate_reference_manifest(
        manifest,
        bundle_root=root,
        bundled=bundled,
        require_redistribution_allowed=require_redistribution_allowed,
        require_all_files_listed=require_all_files_listed,
    )


def parse_reference_manifest_payload(
    payload: dict[str, object],
    *,
    context: str,
) -> ReferenceManifest:
    """Parse a JSON object into the typed manifest model."""

    _require_fields(payload, required_fields=_REQUIRED_MANIFEST_FIELDS, context=context)
    files = _parse_file_manifests(payload.get("files"), context=f"{context}.files")
    return ReferenceManifest(
        reference_id=_require_string(payload, key="reference_id", context=context),
        display_name=_require_string(payload, key="display_name", context=context),
        organism=_require_string(payload, key="organism", context=context),
        taxonomy_id=_optional_int(payload, key="taxonomy_id", context=context),
        protein_namespace=_require_string(
            payload,
            key="protein_namespace",
            context=context,
        ),
        reference_version=_require_string(
            payload,
            key="reference_version",
            context=context,
        ),
        source_name=_require_string(payload, key="source_name", context=context),
        source_version=_optional_string(payload, key="source_version", context=context),
        source_url=_optional_string(payload, key="source_url", context=context),
        retrieved_at=_require_date(payload, key="retrieved_at", context=context),
        table_sha256=_require_string(payload, key="table_sha256", context=context),
        source_publication=_optional_string(
            payload,
            key="source_publication",
            context=context,
        ),
        license_name=_optional_string(payload, key="license_name", context=context),
        license_url=_optional_string(payload, key="license_url", context=context),
        redistribution_status=_require_redistribution_status(
            payload,
            key="redistribution_status",
            context=context,
        ),
        redistribution_notes=_require_string(
            payload,
            key="redistribution_notes",
            context=context,
        ),
        derived_from=_require_string_tuple(
            payload,
            key="derived_from",
            context=context,
        ),
        generated_by=_require_string(payload, key="generated_by", context=context),
        generated_at_utc=_require_string(
            payload,
            key="generated_at_utc",
            context=context,
        ),
        manifest_schema_version=_require_string(
            payload,
            key="manifest_schema_version",
            context=context,
        ),
        files=files,
        sequence_context_policy=_optional_string(
            payload,
            key="sequence_context_policy",
            context=context,
        ),
        sequence_window_length=_optional_int(
            payload,
            key="sequence_window_length",
            context=context,
        ),
        sequence_center_index=_optional_int(
            payload,
            key="sequence_center_index",
            context=context,
        ),
        allowed_sequence_alphabet=_optional_string(
            payload,
            key="allowed_sequence_alphabet",
            context=context,
        ),
        organism_common_name=_optional_string(
            payload,
            key="organism_common_name",
            context=context,
        ),
        supports=_optional_string_tuple(payload, key="supports", context=context),
        limitations=_optional_string_tuple(payload, key="limitations", context=context),
    )


def validate_reference_manifest(
    manifest: ReferenceManifest,
    *,
    bundle_root: str | Path,
    bundled: bool = False,
    require_redistribution_allowed: bool = False,
    require_all_files_listed: bool = False,
) -> ReferenceManifest:
    """Validate manifest semantics and verify all declared file hashes."""

    if not isinstance(manifest, ReferenceManifest):
        raise ReferenceManifestError("reference manifest must be ReferenceManifest")
    root = Path(bundle_root)
    _validate_required_manifest_values(manifest)
    _validate_sequence_context(manifest)
    if (
        bundled or require_redistribution_allowed
    ) and manifest.redistribution_status is not RedistributionStatus.APPROVED:
        raise ReferenceManifestError(
            "bundled reference manifest redistribution_status must be 'approved' "
            f"for release-gate validation: {manifest.reference_id} declares "
            f"{manifest.redistribution_status.value!r}"
        )
    listed_files: set[Path] = set()
    for file_manifest in manifest.files:
        resolved_path = _resolve_manifest_file_path(
            root,
            file_manifest.relative_path,
            context=manifest.reference_id,
        )
        listed_files.add(resolved_path)
        _validate_file_manifest(file_manifest, path=resolved_path)
    _validate_table_sha256_matches_declared_file(manifest)
    if require_all_files_listed:
        _validate_all_bundle_files_listed(
            root=root,
            listed_files=listed_files,
            context=manifest.reference_id,
        )
    return manifest


def validate_bundled_reference_manifests(
    reference_bundles_root: str | Path,
    *,
    require_redistribution_allowed: bool = True,
    require_all_files_listed: bool = True,
) -> tuple[ReferenceManifest, ...]:
    """Validate every bundled reference manifest under a package data root."""

    root = Path(reference_bundles_root)
    if not root.is_dir():
        raise ReferenceManifestError(f"bundled reference root does not exist: {root}")
    manifests: list[ReferenceManifest] = []
    for bundle_root in _iter_reference_bundle_roots(root):
        manifest_path = bundle_root / "manifest.json"
        if not manifest_path.is_file():
            raise ReferenceManifestError(
                f"bundled reference is missing manifest: {bundle_root}"
            )
        manifests.append(
            load_reference_manifest(
                manifest_path,
                bundle_root=bundle_root,
                bundled=require_redistribution_allowed,
                require_redistribution_allowed=require_redistribution_allowed,
                require_all_files_listed=require_all_files_listed,
            )
        )
    return tuple(manifests)


def _validate_required_manifest_values(manifest: ReferenceManifest) -> None:
    for field_name in (
        "reference_id",
        "organism",
        "protein_namespace",
        "reference_version",
        "source_name",
        "table_sha256",
    ):
        value = getattr(manifest, field_name)
        if not isinstance(value, str) or not value.strip():
            raise ReferenceManifestError(
                f"reference manifest {field_name} must be a non-empty string"
            )
    if not _SHA256_PATTERN.match(manifest.table_sha256):
        raise ReferenceManifestError(
            "reference manifest table_sha256 is missing or invalid"
        )
    if not isinstance(manifest.redistribution_status, RedistributionStatus):
        raise ReferenceManifestError(
            "reference manifest redistribution_status must be a RedistributionStatus"
        )
    if manifest.redistribution_status is RedistributionStatus.APPROVED and (
        manifest.license_name is None or manifest.license_url is None
    ):
        raise ReferenceManifestError(
            "reference manifest license_name and license_url are required when "
            "redistribution_status is 'approved'"
        )
    if not manifest.files:
        raise ReferenceManifestError("reference manifest files must not be empty")


def _validate_sequence_context(manifest: ReferenceManifest) -> None:
    length = manifest.sequence_window_length
    center = manifest.sequence_center_index
    if length is not None and center is None:
        raise ReferenceManifestError(
            "sequence-aware reference manifest declares sequence_window_length "
            "without sequence_center_index"
        )
    if center is not None and length is None:
        raise ReferenceManifestError(
            "sequence-aware reference manifest declares sequence_center_index "
            "without sequence_window_length"
        )
    if length is None or center is None:
        return
    if length <= 0:
        raise ReferenceManifestError(
            "reference manifest sequence_window_length must be > 0"
        )
    if center < 0 or center >= length:
        raise ReferenceManifestError(
            "reference manifest sequence_center_index must be within "
            "sequence_window_length"
        )


def _validate_file_manifest(
    file_manifest: ReferenceFileManifest,
    *,
    path: Path,
) -> None:
    if not path.is_file():
        raise ReferenceManifestError(
            "reference manifest listed file does not exist: "
            f"{file_manifest.relative_path}"
        )
    if not file_manifest.role:
        raise ReferenceManifestError(
            f"reference manifest file role must be non-empty: {file_manifest.relative_path}"
        )
    if not file_manifest.format:
        raise ReferenceManifestError(
            "reference manifest file format must be non-empty: "
            f"{file_manifest.relative_path}"
        )
    if not _SHA256_PATTERN.match(file_manifest.sha256):
        raise ReferenceManifestError(
            f"reference manifest file sha256 is missing or invalid: {file_manifest.relative_path}"
        )
    actual_hash = _sha256_file(path)
    if actual_hash != file_manifest.sha256:
        raise ReferenceManifestError(
            "reference manifest file hash mismatch for "
            f"{file_manifest.relative_path}: expected {file_manifest.sha256}, "
            f"actual {actual_hash}"
        )


def _validate_table_sha256_matches_declared_file(manifest: ReferenceManifest) -> None:
    declared_file_hashes = {item.sha256 for item in manifest.files}
    if manifest.table_sha256 not in declared_file_hashes:
        raise ReferenceManifestError(
            "reference manifest table_sha256 must match one declared file sha256"
        )


def _validate_all_bundle_files_listed(
    *,
    root: Path,
    listed_files: set[Path],
    context: str,
) -> None:
    manifest_path = (root / "manifest.json").resolve()
    unlisted = [
        path
        for path in sorted(root.rglob("*"))
        if path.is_file()
        and path.resolve() not in listed_files
        and path.resolve() != manifest_path
    ]
    if unlisted:
        preview = ", ".join(str(path.relative_to(root)) for path in unlisted[:5])
        raise ReferenceManifestError(
            f"bundled reference manifest {context} does not list bundled file(s): "
            f"{preview}"
        )


def _resolve_manifest_file_path(
    root: Path, relative_path: str, *, context: str
) -> Path:
    raw_path = Path(relative_path)
    if raw_path.is_absolute():
        raise ReferenceManifestError(
            f"reference manifest {context} file path must be relative: {relative_path}"
        )
    resolved_root = root.resolve()
    resolved_path = (root / raw_path).resolve()
    try:
        resolved_path.relative_to(resolved_root)
    except ValueError as exc:
        raise ReferenceManifestError(
            f"reference manifest {context} file path escapes bundle root: {relative_path}"
        ) from exc
    return resolved_path


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _iter_reference_bundle_roots(root: Path) -> tuple[Path, ...]:
    bundle_roots: list[Path] = []
    for organism_root in sorted(path for path in root.iterdir() if path.is_dir()):
        bundle_roots.extend(
            sorted(path for path in organism_root.iterdir() if path.is_dir())
        )
    return tuple(bundle_roots)


def _parse_file_manifests(
    value: object, *, context: str
) -> tuple[ReferenceFileManifest, ...]:
    if not isinstance(value, list):
        raise ReferenceManifestError(f"{context} must be an array")
    if not value:
        raise ReferenceManifestError(f"{context} must not be empty")
    files: list[ReferenceFileManifest] = []
    for index, item in enumerate(value):
        if not isinstance(item, dict):
            raise ReferenceManifestError(f"{context}[{index}] must be an object")
        file_payload = cast(dict[str, object], item)
        _require_fields(
            file_payload,
            required_fields=_REQUIRED_FILE_FIELDS,
            context=f"{context}[{index}]",
        )
        files.append(
            ReferenceFileManifest(
                relative_path=_require_string(
                    file_payload,
                    key="relative_path",
                    context=f"{context}[{index}]",
                ),
                role=_require_string(
                    file_payload,
                    key="role",
                    context=f"{context}[{index}]",
                ),
                format=_require_string(
                    file_payload,
                    key="format",
                    context=f"{context}[{index}]",
                ),
                sha256=_require_string(
                    file_payload,
                    key="sha256",
                    context=f"{context}[{index}]",
                ),
                row_count=_optional_int(
                    file_payload,
                    key="row_count",
                    context=f"{context}[{index}]",
                ),
                column_names=_optional_column_names(
                    file_payload,
                    key="column_names",
                    context=f"{context}[{index}]",
                ),
            )
        )
    return tuple(files)


def _require_fields(
    payload: dict[str, object],
    *,
    required_fields: frozenset[str],
    context: str,
) -> None:
    missing = sorted(field for field in required_fields if field not in payload)
    if missing:
        raise ReferenceManifestError(
            f"reference manifest is missing required field(s) for {context}: "
            f"{', '.join(missing)}"
        )


def _require_string(payload: dict[str, object], *, key: str, context: str) -> str:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ReferenceManifestError(
            f"reference manifest {key} must be a non-empty string for {context}"
        )
    return value.strip()


def _optional_string(
    payload: dict[str, object], *, key: str, context: str
) -> str | None:
    if key not in payload or payload.get(key) is None:
        return None
    value = payload.get(key)
    if not isinstance(value, str):
        raise ReferenceManifestError(
            f"reference manifest {key} must be a string or null for {context}"
        )
    text = value.strip()
    return text if text else None


def _require_date(payload: dict[str, object], *, key: str, context: str) -> date:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ReferenceManifestError(
            f"reference manifest {key} must be YYYY-MM-DD for {context}"
        )
    try:
        return date.fromisoformat(value.strip())
    except ValueError as exc:
        raise ReferenceManifestError(
            f"reference manifest {key} must be YYYY-MM-DD for {context}"
        ) from exc


def _require_redistribution_status(
    payload: dict[str, object], *, key: str, context: str
) -> RedistributionStatus:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        allowed = ", ".join(item.value for item in RedistributionStatus)
        raise ReferenceManifestError(
            f"reference manifest {key} must be one of {allowed} for {context}"
        )
    try:
        return RedistributionStatus(value.strip())
    except ValueError as exc:
        allowed = ", ".join(item.value for item in RedistributionStatus)
        raise ReferenceManifestError(
            f"reference manifest {key} must be one of {allowed} for {context}"
        ) from exc


def _optional_int(payload: dict[str, object], *, key: str, context: str) -> int | None:
    if key not in payload or payload.get(key) is None:
        return None
    value = payload.get(key)
    if isinstance(value, bool) or not isinstance(value, int):
        raise ReferenceManifestError(
            f"reference manifest {key} must be an integer or null for {context}"
        )
    return int(value)


def _require_string_tuple(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
) -> tuple[str, ...]:
    if key not in payload:
        raise ReferenceManifestError(
            f"reference manifest {key} is required for {context}"
        )
    return _string_tuple(payload.get(key), key=key, context=context, allow_empty=False)


def _optional_string_tuple(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
) -> tuple[str, ...]:
    if key not in payload or payload.get(key) is None:
        return ()
    return _string_tuple(payload.get(key), key=key, context=context, allow_empty=True)


def _string_tuple(
    value: object,
    *,
    key: str,
    context: str,
    allow_empty: bool,
) -> tuple[str, ...]:
    if not isinstance(value, list):
        raise ReferenceManifestError(
            f"reference manifest {key} must be an array of strings for {context}"
        )
    resolved: list[str] = []
    for index, item in enumerate(value):
        if not isinstance(item, str) or not item.strip():
            raise ReferenceManifestError(
                f"reference manifest {key}[{index}] must be a non-empty string "
                f"for {context}"
            )
        resolved.append(item.strip())
    if not allow_empty and not resolved:
        raise ReferenceManifestError(
            f"reference manifest {key} must not be empty for {context}"
        )
    return tuple(resolved)


def _optional_column_names(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
) -> tuple[str, ...] | None:
    if key not in payload or payload.get(key) is None:
        return None
    value = payload.get(key)
    if not isinstance(value, list):
        raise ReferenceManifestError(
            f"reference manifest {key} must be an array of strings or null for {context}"
        )
    resolved: list[str] = []
    for index, item in enumerate(value):
        if not isinstance(item, str) or not item.strip():
            raise ReferenceManifestError(
                f"reference manifest {key}[{index}] must be a non-empty string "
                f"for {context}"
            )
        resolved.append(item.strip())
    return tuple(resolved)


__all__ = [
    "REFERENCE_MANIFEST_SCHEMA_VERSION",
    "load_reference_manifest",
    "parse_reference_manifest_payload",
    "validate_bundled_reference_manifests",
    "validate_reference_manifest",
]
