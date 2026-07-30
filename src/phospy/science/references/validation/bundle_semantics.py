"""Reference manifest semantic and bundled-resource validation."""

from __future__ import annotations

from pathlib import Path

from phospy.science.references.errors import ReferenceManifestError
from phospy.science.references.manifest import RedistributionStatus, ReferenceManifest
from phospy.science.references.validation._diagnostics import (
    _format_release_gate_failure,
)
from phospy.science.references.validation.manifest_schema import load_reference_manifest
from phospy.science.references.validation.redistribution_policy import (
    _validate_release_gate_redistribution_approval,
)
from phospy.science.references.validation.resource_integrity import (
    _SHA256_PATTERN,
    _iter_reference_bundle_roots,
    _resolve_manifest_file_path,
    _validate_all_bundle_files_listed,
    _validate_file_manifest,
    _validate_table_sha256_matches_declared_file,
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
    release_gate = bundled or require_redistribution_allowed
    _validate_required_manifest_values(manifest, release_gate=release_gate)
    _validate_redistribution_allowed_consistency(manifest, release_gate=release_gate)
    _validate_sequence_context(manifest)
    if release_gate:
        _validate_release_gate_redistribution_approval(manifest, bundle_root=root)
    listed_files: set[Path] = set()
    for file_manifest in manifest.files:
        resolved_path = _resolve_manifest_file_path(
            root,
            file_manifest.relative_path,
            context=manifest.reference_id,
        )
        listed_files.add(resolved_path)
        _validate_file_manifest(file_manifest, path=resolved_path, manifest=manifest)
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


def _validate_required_manifest_values(
    manifest: ReferenceManifest,
    *,
    release_gate: bool,
) -> None:
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
        if release_gate:
            field = "license_name" if manifest.license_name is None else "license_url"
            actual_value = (
                manifest.license_name
                if manifest.license_name is None
                else manifest.license_url
            )
            raise ReferenceManifestError(
                _format_release_gate_failure(
                    manifest,
                    field=field,
                    actual_value=actual_value,
                    reason=(
                        "approved bundled reference requires license_name and "
                        "license_url"
                    ),
                )
            )
        raise ReferenceManifestError(
            "reference manifest license_name and license_url are required when "
            "redistribution_status is 'approved'"
        )
    if not manifest.files:
        raise ReferenceManifestError("reference manifest files must not be empty")


def _validate_redistribution_allowed_consistency(
    manifest: ReferenceManifest,
    *,
    release_gate: bool,
) -> None:
    raw_allowed = manifest.raw_redistribution_allowed
    if raw_allowed is None:
        return
    expected = manifest.redistribution_status is RedistributionStatus.APPROVED
    if raw_allowed is expected:
        return
    reason = (
        "redistribution_allowed must be true when redistribution_status is 'approved'"
        if expected
        else "redistribution_allowed must be false for non-releasable "
        "redistribution_status values"
    )
    if release_gate:
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_allowed",
                actual_value=raw_allowed,
                reason=f"{reason}; got {raw_allowed!r}",
            )
        )
    raise ReferenceManifestError(
        "reference manifest redistribution_allowed contradicts "
        f"redistribution_status for {manifest.reference_id!r}: {reason}; "
        f"got {raw_allowed!r}"
    )


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
