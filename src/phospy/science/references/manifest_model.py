"""Reference manifest aggregate value model and compatibility aliases."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date

from phospy.provenance.models import JsonValue
from phospy.science.references.manifest_common import _coerce_date, _optional_string
from phospy.science.references.manifest_files import (
    ReferenceFileManifest,
    SequenceWindowDefinition,
    _source_file_key,
)
from phospy.science.references.manifest_policy import RedistributionStatus
from phospy.science.references.redistribution import RedistributionEvidence


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


__all__ = [
    "ReferenceManifest",
]
