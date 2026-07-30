"""Reference bundle file-path and digest integrity validation."""

from __future__ import annotations

import re
from hashlib import sha256
from pathlib import Path, PurePosixPath

from phospy.science.references.errors import ReferenceManifestError
from phospy.science.references.manifest import ReferenceFileManifest, ReferenceManifest
from phospy.science.references.validation._diagnostics import (
    _format_release_gate_failure,
    _format_release_validation_failure,
    _ReferenceDiagnosticContext,
)

_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


def _validate_posix_relative_path(
    value: str,
    *,
    manifest: ReferenceManifest,
    field: str,
) -> None:
    reason = _posix_relative_path_rejection_reason(value)
    if reason is not None:
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field=field,
                actual_value=value,
                reason=reason,
            )
        )


def _posix_relative_path_rejection_reason(value: str) -> str | None:
    if not isinstance(value, str) or not value.strip():
        return "path must be a non-empty relative POSIX path"
    if "\\" in value:
        return "path must use POSIX '/' separators, not backslashes"
    if re.match(r"^[A-Za-z]:", value):
        return "path must be a POSIX relative path, not a drive-qualified path"
    path = PurePosixPath(value)
    if path.is_absolute():
        return "path must be relative"
    parts = path.parts
    if not parts or any(part in {"", ".", ".."} for part in parts):
        return "path must not contain empty, current-directory, or parent segments"
    return None


def _find_repository_file(bundle_root: Path, relative_path: str) -> Path | None:
    path = Path(relative_path)
    for root in (bundle_root, *bundle_root.parents):
        candidate = root / path
        if candidate.is_file():
            return candidate
    return None


def _validate_file_manifest(
    file_manifest: ReferenceFileManifest,
    *,
    path: Path,
    manifest: ReferenceManifest,
) -> None:
    if not path.is_file():
        raise ReferenceManifestError(
            _format_file_validation_failure(
                manifest,
                field=f"files[{file_manifest.relative_path!r}]",
                file_path=file_manifest.relative_path,
                actual_value="missing",
                expected_digest=file_manifest.sha256,
                actual_digest="missing",
                reason="reference manifest listed file does not exist",
            )
        )
    if not file_manifest.role:
        raise ReferenceManifestError(
            _format_file_validation_failure(
                manifest,
                field=f"files[{file_manifest.relative_path!r}].role",
                file_path=file_manifest.relative_path,
                actual_value=file_manifest.role,
                expected_digest="non-empty file role",
                actual_digest=file_manifest.role,
                reason="reference manifest file role must be non-empty",
            )
        )
    if not file_manifest.format:
        raise ReferenceManifestError(
            _format_file_validation_failure(
                manifest,
                field=f"files[{file_manifest.relative_path!r}].format",
                file_path=file_manifest.relative_path,
                actual_value=file_manifest.format,
                expected_digest="non-empty file format",
                actual_digest=file_manifest.format,
                reason="reference manifest file format must be non-empty",
            )
        )
    if not _SHA256_PATTERN.match(file_manifest.sha256):
        raise ReferenceManifestError(
            _format_file_validation_failure(
                manifest,
                field=f"files[{file_manifest.relative_path!r}].sha256",
                file_path=file_manifest.relative_path,
                actual_value=file_manifest.sha256,
                expected_digest="lowercase sha256 hex digest",
                actual_digest=file_manifest.sha256,
                reason=(
                    "reference manifest file sha256 must be a 64-character "
                    "lowercase hexadecimal SHA-256 digest"
                ),
            )
        )
    actual_hash = _sha256_file(path)
    if actual_hash != file_manifest.sha256:
        raise ReferenceManifestError(
            _format_file_validation_failure(
                manifest,
                field=f"files[{file_manifest.relative_path!r}].sha256",
                file_path=file_manifest.relative_path,
                actual_value=actual_hash,
                expected_digest=file_manifest.sha256,
                actual_digest=actual_hash,
                reason=(
                    "reference manifest file hash mismatch; "
                    f"expected {file_manifest.sha256}; actual {actual_hash}"
                ),
            )
        )


def _format_file_validation_failure(
    manifest: ReferenceManifest,
    *,
    field: str,
    file_path: str,
    actual_value: object,
    expected_digest: object,
    actual_digest: object,
    reason: str,
) -> str:
    return _format_release_validation_failure(
        _ReferenceDiagnosticContext.from_manifest(manifest),
        field=field,
        actual_value=actual_value,
        file_path=file_path,
        expected_digest=expected_digest,
        actual_digest=actual_digest,
        reason=reason,
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
