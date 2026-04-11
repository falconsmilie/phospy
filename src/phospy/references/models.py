from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from ..motifs import build_validated_motif_library
from ..types import KinaseMotifSequenceMap, KinaseSubstrateMap
from ..validation.errors import InputCompatibilityError
from ..validation.values.collections import normalize_sequence_mapping


@dataclass(frozen=True, slots=True)
class ReferenceBundleSourceMetadata:
    """Source metadata for one kinase reference bundle."""

    source: str
    reference: str
    version: str | None = None

    def __post_init__(self) -> None:
        if not self.source.strip():
            msg = "ReferenceBundle source_metadata.source must not be empty"
            raise InputCompatibilityError(msg)
        if not self.reference.strip():
            msg = "ReferenceBundle source_metadata.reference must not be empty"
            raise InputCompatibilityError(msg)
        if self.version is not None and not self.version.strip():
            msg = "ReferenceBundle source_metadata.version must not be empty when provided"
            raise InputCompatibilityError(msg)


@dataclass(frozen=True, slots=True)
class ReferenceBundleProvenance:
    """Provenance describing how a kinase reference bundle was resolved."""

    provider: str
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.provider.strip():
            msg = "ReferenceBundle provenance.provider must not be empty"
            raise InputCompatibilityError(msg)
        normalized_notes = tuple(str(note) for note in self.notes)
        if any(not note.strip() for note in normalized_notes):
            msg = "ReferenceBundle provenance.notes must not contain empty entries"
            raise InputCompatibilityError(msg)
        object.__setattr__(self, "notes", normalized_notes)


@dataclass(frozen=True, slots=True, init=False)
class ReferenceBundle:
    """Typed kinase-prior contract between reference resolution and workflow setup."""

    substrate_map: dict[str, tuple[str, ...]]
    motif_sequences: dict[str, tuple[str, ...]]
    species: str
    source_metadata: ReferenceBundleSourceMetadata
    provenance: ReferenceBundleProvenance

    def __init__(
        self,
        *,
        substrate_map: KinaseSubstrateMap,
        motif_sequences: KinaseMotifSequenceMap,
        species: str,
        source_metadata: ReferenceBundleSourceMetadata,
        provenance: ReferenceBundleProvenance,
    ) -> None:
        normalized_substrate_map = normalize_sequence_mapping(
            substrate_map,
            field_name="substrate_map",
            empty_message="ReferenceBundle substrate_map must not be empty",
        )
        normalized_motif_sequences = normalize_sequence_mapping(
            motif_sequences,
            field_name="motif_sequences",
            empty_message="ReferenceBundle motif_sequences must not be empty",
        )
        resolved_species = str(species).strip()
        if not resolved_species:
            msg = "ReferenceBundle species must not be empty"
            raise InputCompatibilityError(msg)

        _validate_reference_mapping_values(
            normalized_substrate_map,
            field_name="substrate_map",
        )
        _validate_reference_mapping_values(
            normalized_motif_sequences,
            field_name="motif_sequences",
        )

        substrate_kinases = set(normalized_substrate_map)
        motif_kinases = set(normalized_motif_sequences)
        if substrate_kinases != motif_kinases:
            missing_in_motifs = sorted(substrate_kinases - motif_kinases)
            missing_in_substrates = sorted(motif_kinases - substrate_kinases)
            parts: list[str] = []
            if missing_in_motifs:
                parts.append(
                    "missing from motif_sequences: " + ", ".join(missing_in_motifs)
                )
            if missing_in_substrates:
                parts.append(
                    "missing from substrate_map: " + ", ".join(missing_in_substrates)
                )
            msg = "ReferenceBundle kinase sets must match exactly"
            if parts:
                msg = f"{msg} ({'; '.join(parts)})"
            raise InputCompatibilityError(msg)

        build_validated_motif_library(
            motif_sequences=normalized_motif_sequences,
            context="ReferenceBundle motif_sequences",
        )

        object.__setattr__(self, "substrate_map", dict(normalized_substrate_map))
        object.__setattr__(self, "motif_sequences", dict(normalized_motif_sequences))
        object.__setattr__(self, "species", resolved_species)
        object.__setattr__(self, "source_metadata", source_metadata)
        object.__setattr__(self, "provenance", provenance)


@runtime_checkable
class ReferenceProvider(Protocol):
    """Protocol for resolving kinase prior inputs into a ReferenceBundle."""

    def resolve(
        self,
        *,
        species: str,
        reference: str = "auto",
    ) -> ReferenceBundle: ...


def _validate_reference_mapping_values(
    mapping: dict[str, tuple[str, ...]],
    *,
    field_name: str,
) -> None:
    empty_kinases = sorted(kinase for kinase, values in mapping.items() if not values)
    if empty_kinases:
        msg = (
            f"ReferenceBundle {field_name} entries must not be empty: "
            f"{', '.join(empty_kinases)}"
        )
        raise InputCompatibilityError(msg)
