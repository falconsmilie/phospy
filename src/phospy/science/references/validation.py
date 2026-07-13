"""Reference-manifest loading and validation."""

from __future__ import annotations

import json
import re
from datetime import date
from hashlib import sha256
from pathlib import Path, PurePosixPath
from typing import cast

from phospy.science.references.errors import ReferenceManifestError
from phospy.science.references.manifest import (
    REFERENCE_MANIFEST_SCHEMA_VERSION,
    RedistributionAttribution,
    RedistributionEvidence,
    RedistributionEvidenceType,
    RedistributionScope,
    RedistributionStatus,
    ReferenceFileManifest,
    ReferenceManifest,
    UpstreamPackageLicenseEvidence,
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
_REQUIRED_REDISTRIBUTION_EVIDENCE_FIELDS = frozenset(
    {
        "evidence_type",
        "upstream_package",
        "scope",
        "attribution",
        "independent_database_permission_claimed",
    }
)
_ALLOWED_REDISTRIBUTION_EVIDENCE_FIELDS = frozenset(
    {
        "evidence_type",
        "upstream_package",
        "scope",
        "attribution",
        "independent_database_permission_claimed",
        "evidence_url",
        "verified_at",
        "notes",
    }
)
_REQUIRED_UPSTREAM_PACKAGE_FIELDS = frozenset(
    {
        "package_name",
        "package_version",
        "license_name",
    }
)
_ALLOWED_UPSTREAM_PACKAGE_FIELDS = frozenset(
    {
        "package_name",
        "package_version",
        "license_name",
        "license_url",
    }
)
_REQUIRED_REDISTRIBUTION_SCOPE_FIELDS = frozenset(
    {
        "reference_id",
        "reference_version",
        "applies_to_exact_packaged_files",
        "packaged_files",
        "applies_to_future_bundles",
    }
)
_ALLOWED_REDISTRIBUTION_SCOPE_FIELDS = _REQUIRED_REDISTRIBUTION_SCOPE_FIELDS
_REQUIRED_REDISTRIBUTION_ATTRIBUTION_FIELDS = frozenset(
    {
        "repository_notice_path",
        "bundle_attribution_path",
    }
)
_ALLOWED_REDISTRIBUTION_ATTRIBUTION_FIELDS = _REQUIRED_REDISTRIBUTION_ATTRIBUTION_FIELDS
_PLACEHOLDER_SOURCE_VERSIONS = frozenset(
    {
        "unknown",
        "unspecified",
        "n/a",
        "na",
        "none",
        "null",
        "tbd",
        "not specified",
    }
)
_APPROVAL_CONTRADICTION_PATTERNS = (
    (
        "not legal approval",
        re.compile(r"\bnot\s+legal\s+approval\b", flags=re.IGNORECASE),
    ),
    (
        "redistribution is not approved",
        re.compile(
            r"\bredistribution\s+(?:is\s+)?not\s+approved\b",
            flags=re.IGNORECASE,
        ),
    ),
    (
        "approval has not been independently verified",
        re.compile(
            r"\bapproval\s+(?:has\s+)?not\s+been\s+independently\s+verified\b",
            flags=re.IGNORECASE,
        ),
    ),
    (
        "redistribution remains unresolved",
        re.compile(
            r"\bredistribution\s+remains\s+unresolved\b",
            flags=re.IGNORECASE,
        ),
    ),
    (
        "the exact packaged bundle has not been approved",
        re.compile(
            r"\b(?:the\s+)?exact\s+packaged\s+bundle\s+"
            r"(?:has\s+)?not\s+been\s+approved\b",
            flags=re.IGNORECASE,
        ),
    ),
)
_INDEPENDENT_DATABASE_PERMISSION_PATTERN = re.compile(
    r"\bindependent\s+direct\s+(?:redistribution\s+)?permission\s+from\s+"
    r"(?:phosphositeplus|psp|pride|kinase\s+library|"
    r"(?:an?\s+)?(?:upstream\s+)?(?:scientific\s+)?database|"
    r"other\s+upstream\s+databases?)\b[^.?!;]{0,120}\b"
    r"(?:is\s+)?(?:claimed|granted|approved|obtained|secured)\b",
    flags=re.IGNORECASE,
)
_LICENSE_NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9 .()+/\-]*$")
_LICENSE_PROSE_MARKERS = (
    "approved",
    "approval",
    "permission",
    "redistribution",
    "see ",
    "refer ",
    "unknown",
    "unspecified",
)
_RAT_L6_REFERENCE_ID = "l6_native"
_RAT_L6_REFERENCE_VERSION = "bundled-snapshot-2026-04-16"
_RAT_L6_SOURCE_VERSION = "PhosR 1.20.0"
_RAT_L6_PACKAGE_NAME = "PhosR"
_RAT_L6_PACKAGE_VERSION = "1.20.0"
_RAT_L6_REPOSITORY_NOTICE_PATH = "NOTICE.md"
_RAT_L6_BUNDLE_ATTRIBUTION_PATH = "ATTRIBUTION.md"

_MISSING = object()


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
    redistribution_status = _require_redistribution_status(
        payload,
        key="redistribution_status",
        context=context,
    )
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
        redistribution_status=redistribution_status,
        redistribution_notes=_require_string(
            payload,
            key="redistribution_notes",
            context=context,
        ),
        redistribution_evidence=_optional_redistribution_evidence(
            payload,
            key="redistribution_evidence",
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
        raw_redistribution_allowed=_optional_bool(
            payload,
            key="redistribution_allowed",
            context=context,
        ),
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


def _validate_release_gate_redistribution_approval(
    manifest: ReferenceManifest,
    *,
    bundle_root: Path,
) -> None:
    if manifest.redistribution_status is not RedistributionStatus.APPROVED:
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_status",
                actual_value=manifest.redistribution_status.value,
                reason=(
                    "bundled reference requires redistribution_status 'approved' "
                    "for release-gate validation"
                ),
            )
        )
    _require_release_source_version(manifest)
    _require_release_license_metadata(manifest)
    _validate_approved_text_has_no_contradictions(manifest)
    evidence = manifest.redistribution_evidence
    if evidence is None:
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_evidence",
                actual_value=None,
                reason=(
                    "approved bundled reference requires structured exact-file "
                    "redistribution evidence"
                ),
            )
        )
    if (
        evidence.evidence_type
        is not RedistributionEvidenceType.UPSTREAM_PACKAGE_LICENSE
    ):
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_evidence.evidence_type",
                actual_value=_enum_value(evidence.evidence_type),
                reason=(
                    "approved bundled reference requires upstream_package_license "
                    "redistribution evidence"
                ),
            )
        )
    _validate_release_evidence_upstream_package(manifest, evidence)
    _validate_release_evidence_scope(manifest, evidence)
    _validate_release_evidence_attribution(
        manifest,
        evidence,
        bundle_root=bundle_root,
    )
    _validate_release_independent_database_permission(manifest, evidence)
    _validate_rat_l6_native_policy(manifest, evidence)


def _require_release_source_version(manifest: ReferenceManifest) -> None:
    value = manifest.source_version
    if not isinstance(value, str) or not value.strip():
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="source_version",
                actual_value=value,
                reason=(
                    "bundled reference requires source_version to be a "
                    "non-empty string for release-gate validation"
                ),
            )
        )
    if value.strip().casefold() in _PLACEHOLDER_SOURCE_VERSIONS:
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="source_version",
                actual_value=value,
                reason="bundled reference source_version must not be a placeholder",
            )
        )


def _require_release_license_metadata(manifest: ReferenceManifest) -> None:
    if manifest.license_name is None:
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="license_name",
                actual_value=manifest.license_name,
                reason="approved bundled reference requires license_name",
            )
        )
    if not _is_machine_readable_license_name(manifest.license_name):
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="license_name",
                actual_value=manifest.license_name,
                reason=(
                    "approved bundled reference requires a machine-readable "
                    "license_name"
                ),
            )
        )
    if manifest.license_url is None:
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="license_url",
                actual_value=manifest.license_url,
                reason="approved bundled reference requires license_url",
            )
        )


def _validate_release_evidence_upstream_package(
    manifest: ReferenceManifest,
    evidence: RedistributionEvidence,
) -> None:
    upstream = evidence.upstream_package
    for field_name in ("package_name", "package_version", "license_name"):
        value = getattr(upstream, field_name)
        if not isinstance(value, str) or not value.strip():
            raise ReferenceManifestError(
                _format_release_gate_failure(
                    manifest,
                    field=f"redistribution_evidence.upstream_package.{field_name}",
                    actual_value=value,
                    reason=(
                        "approved bundled reference requires non-empty "
                        f"upstream package {field_name}"
                    ),
                )
            )
    if not _is_machine_readable_license_name(upstream.license_name):
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_evidence.upstream_package.license_name",
                actual_value=upstream.license_name,
                reason=(
                    "approved bundled reference requires machine-readable "
                    "evidence-level license metadata"
                ),
            )
        )
    if upstream.license_name != manifest.license_name:
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_evidence.upstream_package.license_name",
                actual_value=upstream.license_name,
                reason=(
                    f"expected manifest license_name {manifest.license_name!r}; "
                    f"got {upstream.license_name!r}"
                ),
            )
        )
    if (
        upstream.license_url is not None
        and manifest.license_url is not None
        and upstream.license_url != manifest.license_url
    ):
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_evidence.upstream_package.license_url",
                actual_value=upstream.license_url,
                reason=(
                    f"expected manifest license_url {manifest.license_url!r}; "
                    f"got {upstream.license_url!r}"
                ),
            )
        )


def _validate_release_evidence_scope(
    manifest: ReferenceManifest,
    evidence: RedistributionEvidence,
) -> None:
    scope = evidence.scope
    if scope.reference_id != manifest.reference_id:
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_evidence.scope.reference_id",
                actual_value=scope.reference_id,
                reason=(
                    f"expected {manifest.reference_id!r}; got {scope.reference_id!r}"
                ),
            )
        )
    if scope.reference_version != manifest.reference_version:
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_evidence.scope.reference_version",
                actual_value=scope.reference_version,
                reason=(
                    f"expected {manifest.reference_version!r}; "
                    f"got {scope.reference_version!r}"
                ),
            )
        )
    if not scope.applies_to_exact_packaged_files:
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_evidence.scope.applies_to_exact_packaged_files",
                actual_value=scope.applies_to_exact_packaged_files,
                reason=(
                    "approved bundled reference evidence must apply to exact "
                    "packaged files"
                ),
            )
        )
    if scope.applies_to_future_bundles:
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_evidence.scope.applies_to_future_bundles",
                actual_value=scope.applies_to_future_bundles,
                reason="approved bundled reference evidence must not cover future bundles",
            )
        )
    packaged_files = scope.packaged_files
    duplicates = sorted(
        path for path in set(packaged_files) if packaged_files.count(path) > 1
    )
    if duplicates:
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_evidence.scope.packaged_files",
                actual_value=duplicates,
                reason="scope.packaged_files must not contain duplicates",
            )
        )
    for index, packaged_file in enumerate(packaged_files):
        _validate_posix_relative_path(
            packaged_file,
            manifest=manifest,
            field=f"redistribution_evidence.scope.packaged_files[{index}]",
        )
    manifest_files = {item.relative_path for item in manifest.files}
    scoped_files = set(packaged_files)
    if scoped_files != manifest_files:
        missing = sorted(manifest_files - scoped_files)
        extra = sorted(scoped_files - manifest_files)
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_evidence.scope.packaged_files",
                actual_value={
                    "missing": missing,
                    "extra": extra,
                },
                reason=(
                    "scope.packaged_files must exactly match manifest file paths; "
                    f"missing={missing!r}; extra={extra!r}"
                ),
            )
        )


def _validate_release_evidence_attribution(
    manifest: ReferenceManifest,
    evidence: RedistributionEvidence,
    *,
    bundle_root: Path,
) -> None:
    attribution = evidence.attribution
    _validate_posix_relative_path(
        attribution.repository_notice_path,
        manifest=manifest,
        field="redistribution_evidence.attribution.repository_notice_path",
    )
    _validate_posix_relative_path(
        attribution.bundle_attribution_path,
        manifest=manifest,
        field="redistribution_evidence.attribution.bundle_attribution_path",
    )
    manifest_files = {item.relative_path for item in manifest.files}
    if attribution.bundle_attribution_path not in manifest_files:
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_evidence.attribution.bundle_attribution_path",
                actual_value=attribution.bundle_attribution_path,
                reason="bundle attribution path must be included in manifest.files",
            )
        )
    bundle_attribution_file = bundle_root / Path(attribution.bundle_attribution_path)
    if not bundle_attribution_file.is_file():
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_evidence.attribution.bundle_attribution_path",
                actual_value=attribution.bundle_attribution_path,
                reason="bundle attribution file does not exist",
            )
        )
    if _find_repository_file(bundle_root, attribution.repository_notice_path) is None:
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_evidence.attribution.repository_notice_path",
                actual_value=attribution.repository_notice_path,
                reason="repository notice file does not exist",
            )
        )


def _validate_release_independent_database_permission(
    manifest: ReferenceManifest,
    evidence: RedistributionEvidence,
) -> None:
    if evidence.independent_database_permission_claimed:
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field="redistribution_evidence.independent_database_permission_claimed",
                actual_value=evidence.independent_database_permission_claimed,
                reason=(
                    "approved bundled reference must not claim independent direct "
                    "permission from upstream scientific databases"
                ),
            )
        )
    claim = _find_independent_database_permission_claim(manifest)
    if claim is not None:
        field, value = claim
        raise ReferenceManifestError(
            _format_release_gate_failure(
                manifest,
                field=field,
                actual_value=value,
                reason=(
                    "manifest text claims independent direct database permission "
                    "while structured evidence says none is claimed"
                ),
            )
        )


def _validate_rat_l6_native_policy(
    manifest: ReferenceManifest,
    evidence: RedistributionEvidence,
) -> None:
    if not _is_rat_l6_native_snapshot(manifest):
        return
    expectations: tuple[tuple[str, object, object], ...] = (
        ("source_version", manifest.source_version, _RAT_L6_SOURCE_VERSION),
        (
            "redistribution_status",
            manifest.redistribution_status,
            RedistributionStatus.APPROVED,
        ),
        ("redistribution_allowed", manifest.redistribution_allowed, True),
        (
            "redistribution_evidence.evidence_type",
            evidence.evidence_type,
            RedistributionEvidenceType.UPSTREAM_PACKAGE_LICENSE,
        ),
        (
            "redistribution_evidence.upstream_package.package_name",
            evidence.upstream_package.package_name,
            _RAT_L6_PACKAGE_NAME,
        ),
        (
            "redistribution_evidence.upstream_package.package_version",
            evidence.upstream_package.package_version,
            _RAT_L6_PACKAGE_VERSION,
        ),
        (
            "redistribution_evidence.scope.applies_to_exact_packaged_files",
            evidence.scope.applies_to_exact_packaged_files,
            True,
        ),
        (
            "redistribution_evidence.scope.applies_to_future_bundles",
            evidence.scope.applies_to_future_bundles,
            False,
        ),
        (
            "redistribution_evidence.attribution.repository_notice_path",
            evidence.attribution.repository_notice_path,
            _RAT_L6_REPOSITORY_NOTICE_PATH,
        ),
        (
            "redistribution_evidence.attribution.bundle_attribution_path",
            evidence.attribution.bundle_attribution_path,
            _RAT_L6_BUNDLE_ATTRIBUTION_PATH,
        ),
        (
            "redistribution_evidence.independent_database_permission_claimed",
            evidence.independent_database_permission_claimed,
            False,
        ),
    )
    for field, actual, expected in expectations:
        if actual != expected:
            raise ReferenceManifestError(
                _format_release_gate_failure(
                    manifest,
                    field=field,
                    actual_value=_enum_value(actual),
                    reason=f"expected {_enum_value(expected)!r}; got {_enum_value(actual)!r}",
                )
            )


def _validate_approved_text_has_no_contradictions(
    manifest: ReferenceManifest,
) -> None:
    for field, value in _iter_manifest_string_values(manifest.to_payload()):
        phrase = _find_approval_contradiction_phrase(value)
        if phrase is not None:
            raise ReferenceManifestError(
                _format_release_gate_failure(
                    manifest,
                    field=field,
                    actual_value=value,
                    reason=(
                        "approved bundled reference contains contradictory "
                        f"approval text: {phrase!r}"
                    ),
                )
            )


def _find_approval_contradiction_phrase(value: str) -> str | None:
    normalized = " ".join(value.split())
    for phrase, pattern in _APPROVAL_CONTRADICTION_PATTERNS:
        if pattern.search(normalized):
            return phrase
    return None


def _find_independent_database_permission_claim(
    manifest: ReferenceManifest,
) -> tuple[str, str] | None:
    for field, value in _iter_manifest_string_values(manifest.to_payload()):
        if _contains_affirmative_independent_database_permission_claim(value):
            return field, value
    return None


def _contains_affirmative_independent_database_permission_claim(value: str) -> bool:
    normalized = " ".join(value.lower().split())
    for match in _INDEPENDENT_DATABASE_PERMISSION_PATTERN.finditer(normalized):
        prefix = normalized[max(0, match.start() - 60) : match.start()]
        if (
            prefix.strip().endswith("no")
            or "no independent direct" in prefix
            or "not claim" in prefix
            or "does not claim" in prefix
            or "without independent direct" in prefix
        ):
            continue
        return True
    return False


def _iter_manifest_string_values(
    value: object,
    *,
    path: str = "",
) -> tuple[tuple[str, str], ...]:
    fields: list[tuple[str, str]] = []
    if isinstance(value, str):
        fields.append((path, value))
    elif isinstance(value, dict):
        for key, item in value.items():
            child_path = str(key) if not path else f"{path}.{key}"
            fields.extend(_iter_manifest_string_values(item, path=child_path))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            child_path = f"{path}[{index}]"
            fields.extend(_iter_manifest_string_values(item, path=child_path))
    return tuple(fields)


def _is_machine_readable_license_name(value: str) -> bool:
    normalized = " ".join(value.split())
    if not normalized or len(normalized) > 120:
        return False
    if not _LICENSE_NAME_PATTERN.match(normalized):
        return False
    lowered = normalized.casefold()
    return not any(marker in lowered for marker in _LICENSE_PROSE_MARKERS)


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


def _is_rat_l6_native_snapshot(manifest: ReferenceManifest) -> bool:
    organism = (manifest.organism_common_name or manifest.organism).casefold()
    return (
        organism in {"rat", "rattus norvegicus"}
        and manifest.reference_id == _RAT_L6_REFERENCE_ID
        and manifest.reference_version == _RAT_L6_REFERENCE_VERSION
    )


def _enum_value(value: object) -> object:
    if isinstance(value, (RedistributionStatus, RedistributionEvidenceType)):
        return value.value
    return value


def _format_release_gate_failure(
    manifest: ReferenceManifest,
    *,
    field: str,
    actual_value: object = _MISSING,
    reason: str,
) -> str:
    actual = "" if actual_value is _MISSING else f", actual_value={actual_value!r}"
    return (
        "Reference release validation failed: "
        f"reference_id={manifest.reference_id!r}, "
        f"display_name={manifest.display_name!r}, "
        f"organism={manifest.organism!r}, "
        f"namespace={manifest.protein_namespace!r}, "
        f"field={field!r}, "
        f"redistribution_status={manifest.redistribution_status.value!r}"
        f"{actual}: "
        f"{reason}"
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
                reason="reference manifest file sha256 is missing or invalid",
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
    reason: str,
) -> str:
    return (
        "Reference manifest file validation failed: "
        f"reference_id={manifest.reference_id!r}, "
        f"display_name={manifest.display_name!r}, "
        f"organism={manifest.organism!r}, "
        f"namespace={manifest.protein_namespace!r}, "
        f"field={field!r}, "
        f"redistribution_status={manifest.redistribution_status.value!r}, "
        f"file={file_path!r}, "
        f"actual_value={actual_value!r}: "
        f"{reason}"
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


def _optional_date(
    payload: dict[str, object], *, key: str, context: str
) -> date | None:
    if key not in payload or payload.get(key) is None:
        return None
    return _require_date(payload, key=key, context=context)


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


def _optional_redistribution_evidence(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
) -> RedistributionEvidence | None:
    if key not in payload or payload.get(key) is None:
        return None
    value = payload.get(key)
    if not isinstance(value, dict):
        raise ReferenceManifestError(
            f"reference manifest {key} must be an object or null for {context}"
        )
    evidence_payload = cast(dict[str, object], value)
    evidence_context = f"{context}.{key}"
    _require_fields(
        evidence_payload,
        required_fields=_REQUIRED_REDISTRIBUTION_EVIDENCE_FIELDS,
        context=evidence_context,
    )
    _reject_unrecognized_fields(
        evidence_payload,
        allowed_fields=_ALLOWED_REDISTRIBUTION_EVIDENCE_FIELDS,
        context=evidence_context,
    )
    return RedistributionEvidence(
        evidence_type=_require_redistribution_evidence_type(
            evidence_payload,
            key="evidence_type",
            context=evidence_context,
        ),
        upstream_package=_require_upstream_package_license_evidence(
            evidence_payload,
            context=evidence_context,
        ),
        scope=_require_redistribution_scope(
            evidence_payload,
            context=evidence_context,
        ),
        attribution=_require_redistribution_attribution(
            evidence_payload,
            context=evidence_context,
        ),
        independent_database_permission_claimed=_require_bool(
            evidence_payload,
            key="independent_database_permission_claimed",
            context=evidence_context,
        ),
        evidence_url=_optional_string(
            evidence_payload,
            key="evidence_url",
            context=evidence_context,
        ),
        verified_at=_optional_date(
            evidence_payload,
            key="verified_at",
            context=evidence_context,
        ),
        notes=_optional_string(
            evidence_payload,
            key="notes",
            context=evidence_context,
        ),
    )


def _require_upstream_package_license_evidence(
    payload: dict[str, object],
    *,
    context: str,
) -> UpstreamPackageLicenseEvidence:
    upstream_context = f"{context}.upstream_package"
    upstream_payload = _require_object(
        payload,
        key="upstream_package",
        context=context,
    )
    _require_fields(
        upstream_payload,
        required_fields=_REQUIRED_UPSTREAM_PACKAGE_FIELDS,
        context=upstream_context,
    )
    _reject_unrecognized_fields(
        upstream_payload,
        allowed_fields=_ALLOWED_UPSTREAM_PACKAGE_FIELDS,
        context=upstream_context,
    )
    return UpstreamPackageLicenseEvidence(
        package_name=_require_string(
            upstream_payload,
            key="package_name",
            context=upstream_context,
        ),
        package_version=_require_string(
            upstream_payload,
            key="package_version",
            context=upstream_context,
        ),
        license_name=_require_string(
            upstream_payload,
            key="license_name",
            context=upstream_context,
        ),
        license_url=_optional_string(
            upstream_payload,
            key="license_url",
            context=upstream_context,
        ),
    )


def _require_redistribution_scope(
    payload: dict[str, object],
    *,
    context: str,
) -> RedistributionScope:
    scope_context = f"{context}.scope"
    scope_payload = _require_object(payload, key="scope", context=context)
    _require_fields(
        scope_payload,
        required_fields=_REQUIRED_REDISTRIBUTION_SCOPE_FIELDS,
        context=scope_context,
    )
    _reject_unrecognized_fields(
        scope_payload,
        allowed_fields=_ALLOWED_REDISTRIBUTION_SCOPE_FIELDS,
        context=scope_context,
    )
    return RedistributionScope(
        reference_id=_require_string(
            scope_payload,
            key="reference_id",
            context=scope_context,
        ),
        reference_version=_require_string(
            scope_payload,
            key="reference_version",
            context=scope_context,
        ),
        applies_to_exact_packaged_files=_require_bool(
            scope_payload,
            key="applies_to_exact_packaged_files",
            context=scope_context,
        ),
        packaged_files=_require_string_tuple(
            scope_payload,
            key="packaged_files",
            context=scope_context,
        ),
        applies_to_future_bundles=_require_bool(
            scope_payload,
            key="applies_to_future_bundles",
            context=scope_context,
        ),
    )


def _require_redistribution_attribution(
    payload: dict[str, object],
    *,
    context: str,
) -> RedistributionAttribution:
    attribution_context = f"{context}.attribution"
    attribution_payload = _require_object(
        payload,
        key="attribution",
        context=context,
    )
    _require_fields(
        attribution_payload,
        required_fields=_REQUIRED_REDISTRIBUTION_ATTRIBUTION_FIELDS,
        context=attribution_context,
    )
    _reject_unrecognized_fields(
        attribution_payload,
        allowed_fields=_ALLOWED_REDISTRIBUTION_ATTRIBUTION_FIELDS,
        context=attribution_context,
    )
    return RedistributionAttribution(
        repository_notice_path=_require_string(
            attribution_payload,
            key="repository_notice_path",
            context=attribution_context,
        ),
        bundle_attribution_path=_require_string(
            attribution_payload,
            key="bundle_attribution_path",
            context=attribution_context,
        ),
    )


def _require_object(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
) -> dict[str, object]:
    value = payload.get(key)
    if not isinstance(value, dict):
        raise ReferenceManifestError(
            f"reference manifest {key} must be an object for {context}"
        )
    return cast(dict[str, object], value)


def _reject_unrecognized_fields(
    payload: dict[str, object],
    *,
    allowed_fields: frozenset[str],
    context: str,
) -> None:
    extra = sorted(field for field in payload if field not in allowed_fields)
    if extra:
        raise ReferenceManifestError(
            f"reference manifest has unrecognized field(s) for {context}: "
            f"{', '.join(extra)}"
        )


def _require_redistribution_evidence_type(
    payload: dict[str, object],
    *,
    key: str,
    context: str,
) -> RedistributionEvidenceType:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        allowed = ", ".join(item.value for item in RedistributionEvidenceType)
        raise ReferenceManifestError(
            f"reference manifest {key} must be one of {allowed} for {context}"
        )
    try:
        return RedistributionEvidenceType(value.strip())
    except ValueError as exc:
        allowed = ", ".join(item.value for item in RedistributionEvidenceType)
        raise ReferenceManifestError(
            f"reference manifest {key} must be one of {allowed} for {context}"
        ) from exc


def _require_bool(payload: dict[str, object], *, key: str, context: str) -> bool:
    value = payload.get(key)
    if not isinstance(value, bool):
        raise ReferenceManifestError(
            f"reference manifest {key} must be true or false for {context}"
        )
    return value


def _optional_bool(
    payload: dict[str, object], *, key: str, context: str
) -> bool | None:
    if key not in payload or payload.get(key) is None:
        return None
    return _require_bool(payload, key=key, context=context)


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
