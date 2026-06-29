"""Validated reference-provenance manifest models."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date

from phospy.provenance.models import JsonValue

REFERENCE_MANIFEST_SCHEMA_VERSION = "1.0"


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
        object.__setattr__(self, "sha256", str(self.sha256).strip().lower())
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
    source_license: str
    redistribution_allowed: bool
    redistribution_notes: str
    derived_from: tuple[str, ...]
    generated_by: str
    generated_at_utc: str
    manifest_schema_version: str
    files: tuple[ReferenceFileManifest, ...]
    source_version: str | None = None
    source_url: str | None = None
    source_publication: str | None = None
    source_license_url: str | None = None
    sequence_context_policy: str | None = None
    sequence_window_length: int | None = None
    sequence_center_index: int | None = None
    allowed_sequence_alphabet: str | None = None
    organism_common_name: str | None = None
    supports: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for field_name in (
            "reference_id",
            "display_name",
            "organism",
            "protein_namespace",
            "reference_version",
            "source_name",
            "source_license",
            "redistribution_notes",
            "generated_by",
            "generated_at_utc",
            "manifest_schema_version",
        ):
            object.__setattr__(self, field_name, str(getattr(self, field_name)).strip())
        object.__setattr__(
            self, "redistribution_allowed", bool(self.redistribution_allowed)
        )
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
        object.__setattr__(
            self,
            "source_publication",
            _optional_string(self.source_publication),
        )
        object.__setattr__(
            self,
            "source_license_url",
            _optional_string(self.source_license_url),
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

        return self.source_license

    @property
    def license_url(self) -> str | None:
        """Compatibility alias for older reference-bundle metadata."""

        return self.source_license_url

    @property
    def redistribution_status(self) -> str:
        """Compatibility alias for older reference-bundle metadata."""

        return self.redistribution_notes

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
    def retrieved_at(self) -> date:
        """Compatibility date derived from the manifest generation timestamp."""

        return date.fromisoformat(self.generated_at_utc[:10])

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
            "source_publication": self.source_publication,
            "source_license": self.source_license,
            "source_license_url": self.source_license_url,
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
        payload.update(
            {
                "bundle_id": self.bundle_id,
                "identifier_namespace": self.identifier_namespace,
                "license": self.license,
                "license_url": self.license_url,
                "redistribution_status": self.redistribution_status,
                "retrieved_at": self.retrieved_at.isoformat(),
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


def _source_file_key(file_manifest: ReferenceFileManifest) -> str:
    normalized = file_manifest.role.strip().lower().replace("-", "_").replace(" ", "_")
    return normalized or file_manifest.relative_path


__all__ = [
    "REFERENCE_MANIFEST_SCHEMA_VERSION",
    "ReferenceFileManifest",
    "ReferenceManifest",
    "SequenceWindowDefinition",
]
