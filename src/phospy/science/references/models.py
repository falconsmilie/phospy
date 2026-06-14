"""Reference domain models."""

from __future__ import annotations

from dataclasses import InitVar, dataclass
from datetime import date
from enum import Enum
from os import PathLike
from pathlib import Path
from typing import cast

import pandas as pd

from phospy.errors.validation import ReferenceValidationError
from phospy.frames.ownership import export_dataframe, own_dataframe
from phospy.provenance.hashing import fingerprint_table
from phospy.provenance.models import JsonValue, ReferenceProvenance
from phospy.science.references.identifiers import (
    merge_reference_identifier_normalisation_reports,
)
from phospy.tables.references import KinaseSubstrateReference, SiteSequenceReference


class Organism(str, Enum):
    """Public organism identifiers used in dataset/reference contracts.

    Enum membership defines contract syntax. Bundled runtime scientific support
    may be narrower in a given release.
    """

    HUMAN = "human"
    MOUSE = "mouse"
    RAT = "rat"


class ReferencePreset(str, Enum):
    """Built-in organism presets for bundled-reference resolution.

    Enum values define public organism lanes accepted by request contracts.
    Bundled runtime references may cover only a subset in a given release.
    """

    AUTO = "auto"
    HUMAN = "human"
    MOUSE = "mouse"
    RAT = "rat"


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


ReferenceBuildPath = str | Path | PathLike[str]


@dataclass(frozen=True, slots=True)
class ReferenceBundleBuildRequest:
    """Request for building a local-source reference bundle.

    Construction stores caller intent only. File existence, source metadata,
    column mapping, organism compatibility, and reference validity are enforced
    by ``ReferenceBundleBuilder.run(...)``.
    """

    organism: Organism
    kinase_substrate_path: ReferenceBuildPath
    site_sequence_path: ReferenceBuildPath
    source_name: str
    source_version: str
    retrieved_at: date | str
    license: str
    redistribution_status: str
    identifier_namespace: str
    sequence_window: SequenceWindowDefinition | None = None
    bundle_id: str | None = None
    organism_common_name: str | None = None
    supports: tuple[str, ...] = (
        "kinase_workflow",
        "site_sequence_derivation",
    )
    limitations: tuple[str, ...] = (
        "caller-supplied local source files; redistribution governed by request metadata",
    )


@dataclass(frozen=True, slots=True)
class ReferenceManifest:
    """Machine-readable metadata describing one runtime reference bundle."""

    bundle_id: str
    organism: str
    organism_common_name: str | None
    identifier_namespace: str
    source_name: str
    source_version: str
    retrieved_at: date
    license: str
    redistribution_status: str
    sequence_window: SequenceWindowDefinition
    supports: tuple[str, ...]
    limitations: tuple[str, ...]
    source_url: str | None = None
    license_url: str | None = None
    retrieval_method: str | None = None
    redistribution_basis: str | None = None
    source_files: dict[str, JsonValue] | None = None
    provenance_notes: tuple[str, ...] | None = None

    def to_payload(self) -> dict[str, JsonValue]:
        payload: dict[str, JsonValue] = {
            "bundle_id": self.bundle_id,
            "organism": self.organism,
            "organism_common_name": self.organism_common_name,
            "identifier_namespace": self.identifier_namespace,
            "source_name": self.source_name,
            "source_version": self.source_version,
            "retrieved_at": self.retrieved_at.isoformat(),
            "license": self.license,
            "redistribution_status": self.redistribution_status,
            "sequence_window": self.sequence_window.to_payload(),
            "supports": self.supports,
            "limitations": self.limitations,
        }
        if self.source_url is not None:
            payload["source_url"] = self.source_url
        if self.license_url is not None:
            payload["license_url"] = self.license_url
        if self.retrieval_method is not None:
            payload["retrieval_method"] = self.retrieval_method
        if self.redistribution_basis is not None:
            payload["redistribution_basis"] = self.redistribution_basis
        if self.source_files is not None:
            payload["source_files"] = self.source_files
        if self.provenance_notes is not None:
            payload["provenance_notes"] = self.provenance_notes
        return payload


@dataclass(frozen=True, slots=True)
class BundledReferenceLane:
    """Inventory metadata for one packaged bundled reference lane."""

    organism: Organism
    bundle_id: str
    source_name: str
    source_version: str
    retrieved_at: date
    redistribution_status: str
    supports: tuple[str, ...]
    limitations: tuple[str, ...]

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "organism": self.organism.value,
            "bundle_id": self.bundle_id,
            "source_name": self.source_name,
            "source_version": self.source_version,
            "retrieved_at": self.retrieved_at.isoformat(),
            "redistribution_status": self.redistribution_status,
            "supports": self.supports,
            "limitations": self.limitations,
        }


@dataclass(frozen=True, slots=True)
class ReferenceBundle:
    """Resolved workflow reference resources.

    Large kinase-substrate maps are supported. Runtime in downstream workflows
    is primarily controlled by dataset/reference overlap after interpreter and
    scoring-lane filtering, not only by raw map row count.
    """

    organism: Organism
    kinase_substrate_map: pd.DataFrame
    site_sequences: pd.DataFrame
    provenance: ReferenceProvenance | None = None
    manifest: ReferenceManifest | None = None
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        if not isinstance(cast(object, self.organism), Organism):
            raise ReferenceValidationError(
                "references.organism must be an Organism enum value"
            )
        kinase_substrate_map = own_dataframe(
            self.kinase_substrate_map,
            field_name="references.kinase_substrate_map",
            error_type=ReferenceValidationError,
            assume_owned=_assume_owned,
        )
        site_sequences = own_dataframe(
            self.site_sequences,
            field_name="references.site_sequences",
            error_type=ReferenceValidationError,
            assume_owned=_assume_owned,
        )
        kinase_substrate_reference = KinaseSubstrateReference(
            frame=kinase_substrate_map,
            _assume_owned=True,
        )
        site_sequence_reference = SiteSequenceReference(
            frame=site_sequences,
            _assume_owned=True,
        )
        substrate_sites = {
            str(value)
            for value in kinase_substrate_reference.frame["substrate_site"].tolist()
        }
        known_sites = set(site_sequence_reference.frame.index.tolist())
        missing_sequences = sorted(substrate_sites.difference(known_sites))
        if missing_sequences:
            missing_sample = ", ".join(missing_sequences[:10])
            raise ReferenceValidationError(
                "references.site_sequences is missing sequence entries for "
                f"substrate sites in references.kinase_substrate_map: {missing_sample}"
            )
        identifier_normalisation = merge_reference_identifier_normalisation_reports(
            report
            for report in (
                kinase_substrate_reference.identifier_normalisation,
                site_sequence_reference.identifier_normalisation,
            )
            if report is not None
        )
        provenance = self.provenance
        if provenance is None:
            provenance = ReferenceProvenance(
                source_type="explicit",
                organism=self.organism.value,
                bundle_id=None,
                table_fingerprints=(
                    fingerprint_table(
                        kinase_substrate_reference.frame,
                        name="references.kinase_substrate_map",
                    ),
                    fingerprint_table(
                        site_sequence_reference.frame,
                        name="references.site_sequences",
                    ),
                ),
                identifier_normalisation=identifier_normalisation,
            )
        elif not isinstance(cast(object, provenance), ReferenceProvenance):
            raise ReferenceValidationError(
                "references.provenance must be ReferenceProvenance or None"
            )
        elif (
            provenance.source_type == "explicit"
            and provenance.identifier_normalisation is None
        ):
            provenance = ReferenceProvenance(
                source_type=provenance.source_type,
                organism=provenance.organism,
                bundle_id=provenance.bundle_id,
                source_name=provenance.source_name,
                source_version=provenance.source_version,
                retrieved_at=provenance.retrieved_at,
                identifier_namespace=provenance.identifier_namespace,
                sequence_window=provenance.sequence_window,
                manifest=provenance.manifest,
                table_fingerprints=provenance.table_fingerprints,
                identifier_normalisation=identifier_normalisation,
            )
        manifest = self.manifest
        if manifest is not None and not isinstance(
            cast(object, manifest), ReferenceManifest
        ):
            raise ReferenceValidationError(
                "references.manifest must be ReferenceManifest or None"
            )
        object.__setattr__(
            self,
            "kinase_substrate_map",
            kinase_substrate_reference.frame,
        )
        object.__setattr__(
            self,
            "site_sequences",
            site_sequence_reference.frame,
        )
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(self, "manifest", manifest)

    @classmethod
    def _from_owned(
        cls,
        *,
        organism: Organism,
        kinase_substrate_map: pd.DataFrame,
        site_sequences: pd.DataFrame,
        provenance: ReferenceProvenance | None = None,
        manifest: ReferenceManifest | None = None,
    ) -> ReferenceBundle:
        return cls(
            organism=organism,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            provenance=provenance,
            manifest=manifest,
            _assume_owned=True,
        )

    def kinase_substrate_map_dataframe(self) -> pd.DataFrame:
        """Return a kinase-substrate map snapshot isolated from this bundle."""

        return export_dataframe(self.kinase_substrate_map)

    def site_sequences_dataframe(self) -> pd.DataFrame:
        """Return a site-sequence snapshot isolated from this bundle."""

        return export_dataframe(self.site_sequences)
