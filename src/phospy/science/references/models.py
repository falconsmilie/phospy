"""Reference domain models."""

from __future__ import annotations

from dataclasses import InitVar, dataclass, field
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
from phospy.science.references.manifest import (
    ReferenceFileManifest as ReferenceFileManifest,
)
from phospy.science.references.manifest import (
    ReferenceManifest,
    SequenceWindowDefinition,
)


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


ReferenceBuildPath = str | Path | PathLike[str]


def _empty_json_dict() -> dict[str, JsonValue]:
    return {}


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
class ReferenceBundleMissingValueCount:
    """Missing-value count for one important reference field."""

    table_name: str
    field_name: str
    missing_count: int

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "table_name": self.table_name,
            "field_name": self.field_name,
            "missing_count": int(self.missing_count),
        }


@dataclass(frozen=True, slots=True)
class ReferenceBundleTableValidationReport:
    """Validation status for one required reference table."""

    table_name: str
    present: bool
    row_count: int | None
    required_columns: tuple[str, ...]
    present_columns: tuple[str, ...]
    missing_required_columns: tuple[str, ...] = ()
    missing_values: tuple[ReferenceBundleMissingValueCount, ...] = ()

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "table_name": self.table_name,
            "present": bool(self.present),
            "row_count": self.row_count,
            "required_columns": self.required_columns,
            "present_columns": self.present_columns,
            "missing_required_columns": self.missing_required_columns,
            "missing_values": tuple(item.to_payload() for item in self.missing_values),
        }


@dataclass(frozen=True, slots=True)
class ReferenceBundleSourceFileValidationReport:
    """Manifest-declared source-file status for one expected reference role."""

    role: str
    present: bool
    path: str | None = None

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "role": self.role,
            "present": bool(self.present),
            "path": self.path,
        }


@dataclass(frozen=True, slots=True)
class ReferenceBundleValidationReport:
    """Structured reference-bundle validation and provenance summary."""

    bundle_name: str | None
    bundle_version: str | None
    organism: str | None
    organism_common_name: str | None
    identifier_namespace: str | None
    required_tables: tuple[ReferenceBundleTableValidationReport, ...]
    required_source_files: tuple[ReferenceBundleSourceFileValidationReport, ...]
    kinase_substrate_record_count: int
    duplicate_record_count: int
    duplicate_records: tuple[tuple[str, str], ...] = ()
    provenance_fields: dict[str, JsonValue] = field(default_factory=_empty_json_dict)
    compatibility_warnings: tuple[str, ...] = ()

    @property
    def warnings(self) -> tuple[str, ...]:
        """Return compatibility warnings in the report."""

        return self.compatibility_warnings

    def to_payload(self) -> dict[str, JsonValue]:
        required_tables: list[JsonValue] = [
            item.to_payload() for item in self.required_tables
        ]
        required_source_files: list[JsonValue] = [
            item.to_payload() for item in self.required_source_files
        ]
        duplicate_records: list[JsonValue] = []
        for kinase, substrate_site in self.duplicate_records:
            duplicate_records.append(
                {
                    "kinase": kinase,
                    "substrate_site": substrate_site,
                }
            )
        payload: dict[str, JsonValue] = {
            "bundle_name": self.bundle_name,
            "bundle_version": self.bundle_version,
            "organism": self.organism,
            "organism_common_name": self.organism_common_name,
            "identifier_namespace": self.identifier_namespace,
            "required_tables": required_tables,
            "required_source_files": required_source_files,
            "kinase_substrate_record_count": int(self.kinase_substrate_record_count),
            "duplicate_record_count": int(self.duplicate_record_count),
            "duplicate_records": duplicate_records,
            "provenance_fields": dict(self.provenance_fields),
            "compatibility_warnings": self.compatibility_warnings,
        }
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
    _validation_report: ReferenceBundleValidationReport = field(
        init=False,
        repr=False,
        compare=False,
    )

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
        provenance = self.provenance
        manifest = self.manifest
        from phospy.validation.references.bundle import ReferenceBundleValidator

        validation = ReferenceBundleValidator().run(
            organism=self.organism,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            provenance=provenance,
            manifest=manifest,
        )
        identifier_normalisation = validation.identifier_normalisation
        if provenance is None:
            provenance = ReferenceProvenance(
                source_type="explicit",
                organism=self.organism.value,
                bundle_id=None,
                table_fingerprints=(
                    fingerprint_table(
                        validation.kinase_substrate_map,
                        name="references.kinase_substrate_map",
                    ),
                    fingerprint_table(
                        validation.site_sequences,
                        name="references.site_sequences",
                    ),
                ),
                identifier_normalisation=identifier_normalisation,
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
        object.__setattr__(
            self,
            "kinase_substrate_map",
            validation.kinase_substrate_map,
        )
        object.__setattr__(
            self,
            "site_sequences",
            validation.site_sequences,
        )
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(self, "manifest", manifest)
        object.__setattr__(self, "_validation_report", validation.report)

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

    @property
    def validation_report(self) -> ReferenceBundleValidationReport:
        """Return the structured validation report for this bundle."""

        return self._validation_report
