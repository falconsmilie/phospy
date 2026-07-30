"""Reference redistribution and bundled-release policy validation."""

from __future__ import annotations

import re
from datetime import date
from pathlib import Path

from phospy.science.references.errors import ReferenceManifestError
from phospy.science.references.manifest import (
    RedistributionEvidence,
    RedistributionEvidenceType,
    RedistributionStatus,
    ReferenceManifest,
)
from phospy.science.references.validation._diagnostics import (
    _enum_value,
    _format_release_gate_failure,
)
from phospy.science.references.validation.resource_integrity import (
    _find_repository_file,
    _validate_posix_relative_path,
)

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
    _validate_release_evidence_verified_at(manifest, evidence)
    _validate_release_evidence_scope(manifest, evidence)
    _validate_release_evidence_attribution(
        manifest,
        evidence,
        bundle_root=bundle_root,
    )
    _validate_release_independent_database_permission(manifest, evidence)
    _validate_rat_l6_native_policy(manifest, evidence)


def _validate_release_evidence_verified_at(
    manifest: ReferenceManifest,
    evidence: RedistributionEvidence,
) -> None:
    if isinstance(evidence.verified_at, date):
        return
    raise ReferenceManifestError(
        _format_release_gate_failure(
            manifest,
            field="redistribution_evidence.verified_at",
            actual_value=evidence.verified_at,
            reason=(
                "missing or null verified_at; approved bundled evidence requires "
                "an explicit verification date"
            ),
        )
    )


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


def _is_rat_l6_native_snapshot(manifest: ReferenceManifest) -> bool:
    organism = (manifest.organism_common_name or manifest.organism).casefold()
    return (
        organism in {"rat", "rattus norvegicus"}
        and manifest.reference_id == _RAT_L6_REFERENCE_ID
        and manifest.reference_version == _RAT_L6_REFERENCE_VERSION
    )
