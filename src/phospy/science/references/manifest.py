"""Validated reference-provenance manifest models."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date
from enum import Enum

from phospy.provenance.models import JsonValue

REFERENCE_MANIFEST_SCHEMA_VERSION = "1.1"


class RedistributionStatus(str, Enum):
    """Machine-readable redistribution status for reference manifests."""

    APPROVED = "approved"
    EXTERNAL_ONLY = "external_only"
    UNRESOLVED = "unresolved"


class RedistributionEvidenceType(str, Enum):
    """Machine-readable evidence category for exact-file redistribution approval."""

    UPSTREAM_PACKAGE_LICENSE = "upstream_package_license"


@dataclass(frozen=True, slots=True)
class SequenceWindowDefinition:
    """Reference sequence-window definition for centralized site sequences."""

    upstream_residues: int
    downstream_residues: int
    central_residue_required: bool

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "upstream_residues": int(self.upstream_residues),
            "downstream_residues": int(self.downstream_residues),
            "central_residue_required": bool(self.central_residue_required),
        }


@dataclass(frozen=True, slots=True)
class UpstreamPackageLicenseEvidence:
    """Upstream package metadata that supplies the redistribution basis."""

    package_name: str
    package_version: str
    license_name: str
    license_url: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "package_name",
            _required_string(
                self.package_name,
                "redistribution_evidence.upstream_package.package_name",
            ),
        )
        object.__setattr__(
            self,
            "package_version",
            _required_string(
                self.package_version,
                "redistribution_evidence.upstream_package.package_version",
            ),
        )
        object.__setattr__(
            self,
            "license_name",
            _required_string(
                self.license_name,
                "redistribution_evidence.upstream_package.license_name",
            ),
        )
        object.__setattr__(
            self,
            "license_url",
            _optional_string(self.license_url),
        )

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "package_name": self.package_name,
            "package_version": self.package_version,
            "license_name": self.license_name,
            "license_url": self.license_url,
        }


@dataclass(frozen=True, slots=True)
class RedistributionScope:
    """Exact PhosPy bundle snapshot and file scope covered by the evidence."""

    reference_id: str
    reference_version: str
    applies_to_exact_packaged_files: bool
    packaged_files: tuple[str, ...]
    applies_to_future_bundles: bool

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "reference_id",
            _required_string(
                self.reference_id,
                "redistribution_evidence.scope.reference_id",
            ),
        )
        object.__setattr__(
            self,
            "reference_version",
            _required_string(
                self.reference_version,
                "redistribution_evidence.scope.reference_version",
            ),
        )
        if not isinstance(self.applies_to_exact_packaged_files, bool):
            raise ValueError(
                "reference manifest "
                "redistribution_evidence.scope.applies_to_exact_packaged_files "
                "must be bool"
            )
        object.__setattr__(
            self,
            "applies_to_exact_packaged_files",
            self.applies_to_exact_packaged_files,
        )
        packaged_files = tuple(
            _required_string(
                item,
                "redistribution_evidence.scope.packaged_files",
            )
            for item in self.packaged_files
        )
        object.__setattr__(
            self,
            "packaged_files",
            packaged_files,
        )
        if not isinstance(self.applies_to_future_bundles, bool):
            raise ValueError(
                "reference manifest "
                "redistribution_evidence.scope.applies_to_future_bundles "
                "must be bool"
            )
        object.__setattr__(
            self,
            "applies_to_future_bundles",
            self.applies_to_future_bundles,
        )

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "reference_id": self.reference_id,
            "reference_version": self.reference_version,
            "applies_to_exact_packaged_files": bool(
                self.applies_to_exact_packaged_files
            ),
            "packaged_files": list(self.packaged_files),
            "applies_to_future_bundles": bool(self.applies_to_future_bundles),
        }


@dataclass(frozen=True, slots=True)
class RedistributionAttribution:
    """Repository and bundle-local attribution locations for a bundled reference."""

    repository_notice_path: str
    bundle_attribution_path: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "repository_notice_path",
            _required_string(
                self.repository_notice_path,
                "redistribution_evidence.attribution.repository_notice_path",
            ),
        )
        object.__setattr__(
            self,
            "bundle_attribution_path",
            _required_string(
                self.bundle_attribution_path,
                "redistribution_evidence.attribution.bundle_attribution_path",
            ),
        )

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "repository_notice_path": self.repository_notice_path,
            "bundle_attribution_path": self.bundle_attribution_path,
        }


@dataclass(frozen=True, slots=True)
class RedistributionEvidence:
    """Typed exact-snapshot redistribution evidence for release validation."""

    evidence_type: RedistributionEvidenceType | str
    upstream_package: UpstreamPackageLicenseEvidence
    scope: RedistributionScope
    attribution: RedistributionAttribution
    independent_database_permission_claimed: bool
    evidence_url: str | None = None
    verified_at: date | None = None
    notes: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "evidence_type",
            _coerce_redistribution_evidence_type(self.evidence_type),
        )
        if not isinstance(self.upstream_package, UpstreamPackageLicenseEvidence):
            raise ValueError(
                "reference manifest redistribution_evidence.upstream_package "
                "must be UpstreamPackageLicenseEvidence"
            )
        if not isinstance(self.scope, RedistributionScope):
            raise ValueError(
                "reference manifest redistribution_evidence.scope "
                "must be RedistributionScope"
            )
        if not isinstance(self.attribution, RedistributionAttribution):
            raise ValueError(
                "reference manifest redistribution_evidence.attribution "
                "must be RedistributionAttribution"
            )
        if not isinstance(self.independent_database_permission_claimed, bool):
            raise ValueError(
                "reference manifest "
                "redistribution_evidence.independent_database_permission_claimed "
                "must be bool"
            )
        object.__setattr__(
            self,
            "independent_database_permission_claimed",
            self.independent_database_permission_claimed,
        )
        object.__setattr__(self, "evidence_url", _optional_string(self.evidence_url))
        object.__setattr__(
            self,
            "verified_at",
            None if self.verified_at is None else _coerce_date(self.verified_at),
        )
        object.__setattr__(self, "notes", _optional_string(self.notes))

    def to_payload(self) -> dict[str, JsonValue]:
        evidence_type = _coerce_redistribution_evidence_type(self.evidence_type)
        return {
            "evidence_type": evidence_type.value,
            "upstream_package": self.upstream_package.to_payload(),
            "scope": self.scope.to_payload(),
            "attribution": self.attribution.to_payload(),
            "independent_database_permission_claimed": bool(
                self.independent_database_permission_claimed
            ),
            "evidence_url": self.evidence_url,
            "verified_at": (
                None if self.verified_at is None else self.verified_at.isoformat()
            ),
            "notes": self.notes,
        }


@dataclass(frozen=True, slots=True)
class ReferenceFileManifest:
    """Hash-verifiable file metadata for one file in a reference bundle."""

    relative_path: str
    role: str
    format: str
    sha256: str
    row_count: int | None = None
    column_names: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "relative_path", str(self.relative_path).strip())
        object.__setattr__(self, "role", str(self.role).strip())
        object.__setattr__(self, "format", str(self.format).strip())
        object.__setattr__(self, "sha256", str(self.sha256).strip())
        if self.row_count is not None:
            object.__setattr__(self, "row_count", int(self.row_count))
        if self.column_names is not None:
            object.__setattr__(
                self,
                "column_names",
                tuple(str(item) for item in self.column_names),
            )

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "relative_path": self.relative_path,
            "role": self.role,
            "format": self.format,
            "sha256": self.sha256,
            "row_count": self.row_count,
            "column_names": (
                None if self.column_names is None else list(self.column_names)
            ),
        }

    def to_source_file_payload(self) -> dict[str, JsonValue]:
        return {
            "path": self.relative_path,
            "relative_path": self.relative_path,
            "role": self.role,
            "format": self.format,
            "sha256": self.sha256,
            "row_count": self.row_count,
            "column_names": (
                None if self.column_names is None else list(self.column_names)
            ),
        }


@dataclass(frozen=True, slots=True)
class ReferenceManifest:
    """Machine-readable metadata describing one logical reference bundle."""

    reference_id: str
    display_name: str
    organism: str
    taxonomy_id: int | None
    protein_namespace: str
    reference_version: str
    source_name: str
    source_url: str | None
    source_version: str | None
    retrieved_at: date
    table_sha256: str
    license_name: str | None
    license_url: str | None
    redistribution_status: RedistributionStatus
    redistribution_notes: str
    derived_from: tuple[str, ...]
    generated_by: str
    generated_at_utc: str
    manifest_schema_version: str
    files: tuple[ReferenceFileManifest, ...]
    source_publication: str | None = None
    sequence_context_policy: str | None = None
    sequence_window_length: int | None = None
    sequence_center_index: int | None = None
    allowed_sequence_alphabet: str | None = None
    organism_common_name: str | None = None
    supports: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()
    redistribution_evidence: RedistributionEvidence | None = None
    raw_redistribution_allowed: bool | None = None

    def __post_init__(self) -> None:
        for field_name in (
            "reference_id",
            "display_name",
            "organism",
            "protein_namespace",
            "reference_version",
            "source_name",
            "table_sha256",
            "redistribution_notes",
            "generated_by",
            "generated_at_utc",
            "manifest_schema_version",
        ):
            object.__setattr__(self, field_name, str(getattr(self, field_name)).strip())
        object.__setattr__(
            self,
            "derived_from",
            tuple(str(item).strip() for item in self.derived_from),
        )
        object.__setattr__(self, "files", tuple(self.files))
        object.__setattr__(
            self, "source_version", _optional_string(self.source_version)
        )
        object.__setattr__(self, "source_url", _optional_string(self.source_url))
        object.__setattr__(self, "retrieved_at", _coerce_date(self.retrieved_at))
        object.__setattr__(self, "license_name", _optional_string(self.license_name))
        object.__setattr__(self, "license_url", _optional_string(self.license_url))
        object.__setattr__(
            self,
            "redistribution_status",
            _coerce_redistribution_status(self.redistribution_status),
        )
        if self.redistribution_evidence is not None and not isinstance(
            self.redistribution_evidence,
            RedistributionEvidence,
        ):
            raise ValueError(
                "reference manifest redistribution_evidence must be "
                "RedistributionEvidence or None"
            )
        if self.raw_redistribution_allowed is not None and not isinstance(
            self.raw_redistribution_allowed,
            bool,
        ):
            raise ValueError(
                "reference manifest raw_redistribution_allowed must be bool or None"
            )
        object.__setattr__(
            self,
            "source_publication",
            _optional_string(self.source_publication),
        )
        object.__setattr__(
            self,
            "sequence_context_policy",
            _optional_string(self.sequence_context_policy),
        )
        object.__setattr__(
            self,
            "allowed_sequence_alphabet",
            _optional_string(self.allowed_sequence_alphabet),
        )
        object.__setattr__(
            self,
            "organism_common_name",
            _optional_string(self.organism_common_name),
        )
        object.__setattr__(
            self,
            "supports",
            tuple(str(item).strip() for item in self.supports),
        )
        object.__setattr__(
            self,
            "limitations",
            tuple(str(item).strip() for item in self.limitations),
        )

    @property
    def bundle_id(self) -> str:
        """Compatibility alias for older reference-bundle metadata."""

        return self.reference_id

    @property
    def identifier_namespace(self) -> str:
        """Compatibility alias for older reference-bundle metadata."""

        return self.protein_namespace

    @property
    def license(self) -> str:
        """Compatibility alias for older reference-bundle metadata."""

        return "" if self.license_name is None else self.license_name

    @property
    def source_license(self) -> str:
        """Compatibility alias for older reference-bundle metadata."""

        return self.license

    @property
    def source_license_url(self) -> str | None:
        """Compatibility alias for older reference-bundle metadata."""

        return self.license_url

    @property
    def redistribution_allowed(self) -> bool:
        """Compatibility boolean derived from structured redistribution status."""

        return self.redistribution_status is RedistributionStatus.APPROVED

    @property
    def retrieval_method(self) -> str:
        """Compatibility alias for older reference-bundle metadata."""

        return self.generated_by

    @property
    def redistribution_basis(self) -> str:
        """Compatibility alias for older reference-bundle metadata."""

        return self.redistribution_notes

    @property
    def provenance_notes(self) -> tuple[str, ...]:
        """Compatibility alias for older reference-bundle metadata."""

        return (*self.derived_from, self.redistribution_notes)

    @property
    def sequence_window(self) -> SequenceWindowDefinition:
        """Compatibility sequence-window object derived from index fields."""

        if self.sequence_window_length is None or self.sequence_center_index is None:
            return SequenceWindowDefinition(
                upstream_residues=0,
                downstream_residues=0,
                central_residue_required=False,
            )
        return SequenceWindowDefinition(
            upstream_residues=int(self.sequence_center_index),
            downstream_residues=(
                int(self.sequence_window_length) - int(self.sequence_center_index) - 1
            ),
            central_residue_required=True,
        )

    @property
    def source_files(self) -> dict[str, JsonValue]:
        """Compatibility source-file mapping derived from file manifests."""

        return {
            _source_file_key(item): item.to_source_file_payload() for item in self.files
        }

    def to_payload(self) -> dict[str, JsonValue]:
        payload: dict[str, JsonValue] = {
            "reference_id": self.reference_id,
            "display_name": self.display_name,
            "organism": self.organism,
            "taxonomy_id": self.taxonomy_id,
            "protein_namespace": self.protein_namespace,
            "reference_version": self.reference_version,
            "source_name": self.source_name,
            "source_version": self.source_version,
            "source_url": self.source_url,
            "retrieved_at": self.retrieved_at.isoformat(),
            "table_sha256": self.table_sha256,
            "source_publication": self.source_publication,
            "license_name": self.license_name,
            "license_url": self.license_url,
            "redistribution_status": self.redistribution_status.value,
            "redistribution_allowed": self.redistribution_allowed,
            "redistribution_notes": self.redistribution_notes,
            "derived_from": list(self.derived_from),
            "generated_by": self.generated_by,
            "generated_at_utc": self.generated_at_utc,
            "manifest_schema_version": self.manifest_schema_version,
            "files": [item.to_payload() for item in self.files],
            "sequence_context_policy": self.sequence_context_policy,
            "sequence_window_length": self.sequence_window_length,
            "sequence_center_index": self.sequence_center_index,
            "allowed_sequence_alphabet": self.allowed_sequence_alphabet,
            "organism_common_name": self.organism_common_name,
            "supports": list(self.supports),
            "limitations": list(self.limitations),
        }
        if self.redistribution_evidence is not None:
            payload["redistribution_evidence"] = (
                self.redistribution_evidence.to_payload()
            )
        payload.update(
            {
                "bundle_id": self.bundle_id,
                "identifier_namespace": self.identifier_namespace,
                "license": self.license,
                "license_url": self.license_url,
                "source_license": self.source_license,
                "source_license_url": self.source_license_url,
                "retrieval_method": self.retrieval_method,
                "redistribution_basis": self.redistribution_basis,
                "provenance_notes": list(self.provenance_notes),
                "sequence_window": self.sequence_window.to_payload(),
                "source_files": self.source_files,
            }
        )
        return payload


def _optional_string(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text if text else None


def _required_string(value: object, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"reference manifest {field_name} must be non-empty")
    return value.strip()


def _coerce_date(value: date | str) -> date:
    if isinstance(value, date):
        return value
    if not isinstance(value, str) or not value.strip():
        raise ValueError("reference manifest retrieved_at must be YYYY-MM-DD")
    return date.fromisoformat(value.strip())


def _coerce_redistribution_status(
    value: RedistributionStatus | str,
) -> RedistributionStatus:
    if isinstance(value, RedistributionStatus):
        return value
    try:
        return RedistributionStatus(str(value).strip())
    except ValueError as exc:
        allowed = ", ".join(item.value for item in RedistributionStatus)
        raise ValueError(
            f"reference manifest redistribution_status must be one of: {allowed}"
        ) from exc


def _coerce_redistribution_evidence_type(
    value: RedistributionEvidenceType | str,
) -> RedistributionEvidenceType:
    if isinstance(value, RedistributionEvidenceType):
        return value
    try:
        return RedistributionEvidenceType(str(value).strip())
    except ValueError as exc:
        allowed = ", ".join(item.value for item in RedistributionEvidenceType)
        raise ValueError(
            "reference manifest redistribution_evidence.evidence_type must be "
            f"one of: {allowed}"
        ) from exc


def _source_file_key(file_manifest: ReferenceFileManifest) -> str:
    normalized = file_manifest.role.strip().lower().replace("-", "_").replace(" ", "_")
    return normalized or file_manifest.relative_path


__all__ = [
    "REFERENCE_MANIFEST_SCHEMA_VERSION",
    "ReferenceFileManifest",
    "ReferenceManifest",
    "RedistributionAttribution",
    "RedistributionEvidence",
    "RedistributionEvidenceType",
    "RedistributionScope",
    "RedistributionStatus",
    "SequenceWindowDefinition",
    "UpstreamPackageLicenseEvidence",
]
