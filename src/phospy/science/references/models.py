"""Reference domain models."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from datetime import date
from enum import Enum
from os import PathLike
from pathlib import Path
from typing import cast

import pandas as pd

from phospy.errors.validation import ReferenceValidationError
from phospy.frames.comparison import dataframe_equals
from phospy.frames.ownership import export_dataframe, own_dataframe
from phospy.provenance.hashing import fingerprint_table
from phospy.provenance.models import (
    JsonValue,
    ReferenceContextProtocol,
    ReferenceProvenance,
    validate_reference_source_version_agreement,
)
from phospy.provenance.organisms import Organism, normalize_organism, organism_value
from phospy.provenance.reference_context import ReferenceContext
from phospy.science.references.identifiers import (
    ReferenceIdentifierNormalisationReport,
    merge_reference_identifier_normalisation_reports,
)
from phospy.science.references.manifest import (
    RedistributionStatus as RedistributionStatus,
)
from phospy.science.references.manifest import (
    ReferenceFileManifest as ReferenceFileManifest,
)
from phospy.science.references.manifest import (
    ReferenceManifest,
    SequenceWindowDefinition,
)
from phospy.science.tables.references import (
    KinaseSubstrateReference,
    SiteSequenceReference,
)


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


def reference_context_from_provenance(
    provenance: ReferenceProvenance,
    *,
    manifest: ReferenceManifest | None = None,
) -> ReferenceContextProtocol | None:
    """Resolve a reference context from provenance or richer manifest metadata."""

    if provenance.reference_context is not None:
        return provenance.reference_context
    if manifest is not None:
        return reference_context_from_manifest_if_complete(manifest)
    if provenance.manifest is not None:
        return _reference_context_from_manifest_payload(provenance.manifest)
    if (
        provenance.source_name is None
        or provenance.source_version is None
        or provenance.identifier_namespace is None
    ):
        return None
    return ReferenceContext(
        organism=provenance.organism,
        protein_namespace=provenance.identifier_namespace,
        source_name=provenance.source_name,
        source_version=provenance.source_version,
        proteome_version=None,
        reference_table_sha256=None,
    )


def reference_context_from_manifest_if_complete(
    manifest: ReferenceManifest,
) -> ReferenceContext | None:
    """Build reference context only when manifest identity metadata is complete."""

    source_version = manifest.source_version
    if not isinstance(source_version, str) or not source_version.strip():
        return None
    return ReferenceContext.from_manifest(manifest)


def _reference_context_from_manifest_payload(
    payload: Mapping[str, JsonValue],
) -> ReferenceContext | None:
    organism = _payload_optional_text(
        payload.get("organism_common_name")
    ) or _payload_optional_text(payload.get("organism"))
    protein_namespace = _payload_optional_text(
        payload.get("protein_namespace")
    ) or _payload_optional_text(payload.get("identifier_namespace"))
    source_name = _payload_optional_text(payload.get("source_name"))
    source_version = _payload_optional_text(payload.get("source_version"))
    if (
        organism is None
        or protein_namespace is None
        or source_name is None
        or source_version is None
    ):
        return None
    return ReferenceContext(
        organism=organism,
        protein_namespace=protein_namespace,
        source_name=source_name,
        source_version=source_version,
        proteome_version=_payload_optional_text(payload.get("proteome_version")),
        reference_table_sha256=_payload_optional_text(payload.get("table_sha256")),
    )


def _payload_optional_text(value: object) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        return None
    text = value.strip()
    return text or None


def _reference_context_source_version(
    reference_context: ReferenceContextProtocol | None,
) -> str | None:
    if reference_context is None:
        return None
    return reference_context.source_version


def _require_reference_bundle_organism_coherence(
    *,
    organism: Organism,
    provenance: ReferenceProvenance | None,
    manifest: ReferenceManifest | None,
) -> None:
    values: list[tuple[str, object]] = [("references.organism", organism)]
    if provenance is not None:
        values.append(("references.provenance.organism", provenance.organism))
        if provenance.reference_context is not None:
            values.append(
                (
                    "references.provenance.reference_context.organism",
                    provenance.reference_context.organism,
                )
            )
    if manifest is not None:
        values.append(("references.manifest.organism", manifest.organism))
        if manifest.organism_common_name is not None:
            values.append(
                (
                    "references.manifest.organism_common_name",
                    manifest.organism_common_name,
                )
            )
    _require_same_organism_identity(
        values=values,
        conflict_prefix="reference bundle organism identity conflict",
    )


def _require_same_organism_identity(
    *,
    values: list[tuple[str, object]],
    conflict_prefix: str,
) -> None:
    if not values:
        return
    normalized = [
        (
            field_name,
            normalize_organism(
                value,
                field_name=field_name,
                error_type=ReferenceValidationError,
            ),
            value,
        )
        for field_name, value in values
    ]
    expected_field, expected_organism, _ = normalized[0]
    conflicts = [
        (field_name, organism, raw_value)
        for field_name, organism, raw_value in normalized[1:]
        if organism is not expected_organism
    ]
    if not conflicts:
        return
    conflict_text = "; ".join(
        f"{field_name}={_format_organism_value(raw_value)!r}"
        f" resolved_to={organism.value!r}"
        for field_name, organism, raw_value in conflicts
    )
    raise ReferenceValidationError(
        f"{conflict_prefix}; "
        f"{expected_field}={expected_organism.value!r}; {conflict_text}"
    )


def _format_organism_value(value: object) -> str:
    if isinstance(value, Organism):
        return value.value
    return str(value)


@dataclass(frozen=True, slots=True)
class ReferenceBundleBuildRequest:
    """Request for building a local-source reference bundle.

    Construction stores caller intent only. File existence, source metadata,
    column mapping, organism compatibility, and reference validity are enforced
    by ``ReferenceBundleBuilder.run(...)``. ``source_version`` identifies the
    upstream package, database, or caller source. ``reference_version`` is an
    optional local PhosPy snapshot identity consumed by the builder.
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
    reference_version: str | None = None


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


_KINASE_SUBSTRATE_TABLE = "kinase_substrate_map"
_SITE_SEQUENCES_TABLE = "site_sequences"
_KINASE_SUBSTRATE_REQUIRED_COLUMNS = ("kinase", "substrate_site")
_SITE_SEQUENCE_REQUIRED_COLUMNS = ("site_sequence",)
_SOURCE_FILE_ROLE_ALIASES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("kinase_substrate", ("kinase_substrate", "kinase_substrate_map", "substrate_map")),
    ("site_sequences", ("site_sequences", "site_sequence", "site_sequence_table")),
)


@dataclass(frozen=True, slots=True, eq=False)
class ReferenceBundleValidationResult:
    """Validated reference tables plus the public validation report.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit content comparison.
    """

    __hash__ = object.__hash__

    kinase_substrate_map: pd.DataFrame
    site_sequences: pd.DataFrame
    identifier_normalisation: ReferenceIdentifierNormalisationReport | None
    report: ReferenceBundleValidationReport

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another validation result has the same content."""

        if not isinstance(other, ReferenceBundleValidationResult):
            return False
        return (
            dataframe_equals(
                self.kinase_substrate_map,
                other.kinase_substrate_map,
            )
            and dataframe_equals(self.site_sequences, other.site_sequences)
            and self.identifier_normalisation == other.identifier_normalisation
            and self.report == other.report
        )


class ReferenceBundleValidator:
    """Validate the stable `ReferenceBundle` contract."""

    def run(
        self,
        *,
        organism: Organism,
        kinase_substrate_map: object,
        site_sequences: object,
        provenance: ReferenceProvenance | None = None,
        manifest: ReferenceManifest | None = None,
    ) -> ReferenceBundleValidationResult:
        if not isinstance(cast(object, organism), Organism):
            raise ReferenceValidationError(
                "references.organism must be an Organism enum value"
            )
        if not isinstance(kinase_substrate_map, pd.DataFrame):
            raise ReferenceValidationError(
                "references.kinase_substrate_map must be a pandas DataFrame"
            )
        if not isinstance(site_sequences, pd.DataFrame):
            raise ReferenceValidationError(
                "references.site_sequences must be a pandas DataFrame"
            )
        if provenance is not None and not isinstance(
            cast(object, provenance), ReferenceProvenance
        ):
            raise ReferenceValidationError(
                "references.provenance must be ReferenceProvenance or None"
            )
        if manifest is not None and not isinstance(
            cast(object, manifest), ReferenceManifest
        ):
            raise ReferenceValidationError(
                "references.manifest must be ReferenceManifest or None"
            )
        _require_reference_bundle_organism_coherence(
            organism=organism,
            provenance=provenance,
            manifest=manifest,
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
        report = _build_validation_report(
            organism=organism,
            kinase_substrate_map=kinase_substrate_reference.frame,
            site_sequences=site_sequence_reference.frame,
            provenance=provenance,
            manifest=manifest,
        )
        return ReferenceBundleValidationResult(
            kinase_substrate_map=kinase_substrate_reference.frame,
            site_sequences=site_sequence_reference.frame,
            identifier_normalisation=identifier_normalisation,
            report=report,
        )


def _build_validation_report(
    *,
    organism: Organism,
    kinase_substrate_map: pd.DataFrame,
    site_sequences: pd.DataFrame,
    provenance: ReferenceProvenance | None,
    manifest: ReferenceManifest | None,
) -> ReferenceBundleValidationReport:
    duplicate_records = _duplicate_kinase_substrate_records(kinase_substrate_map)
    provenance_fields = _build_provenance_fields(
        organism=organism,
        provenance=provenance,
        manifest=manifest,
    )
    required_source_files = _build_source_file_reports(
        source_files=_resolve_source_files(
            provenance=provenance,
            manifest=manifest,
        )
    )
    warnings = _build_compatibility_warnings(
        organism_common_name=_resolve_organism_common_name(manifest=manifest),
        identifier_namespace=_resolve_identifier_namespace(
            provenance=provenance,
            manifest=manifest,
        ),
        bundle_version=_resolve_bundle_version(
            provenance=provenance,
            manifest=manifest,
        ),
        required_source_files=required_source_files,
        manifest=manifest,
        provenance=provenance,
    )
    return ReferenceBundleValidationReport(
        bundle_name=_resolve_bundle_name(provenance=provenance, manifest=manifest),
        bundle_version=_resolve_bundle_version(
            provenance=provenance,
            manifest=manifest,
        ),
        organism=_resolve_organism(
            organism=organism,
            provenance=provenance,
            manifest=manifest,
        ),
        organism_common_name=_resolve_organism_common_name(manifest=manifest),
        identifier_namespace=_resolve_identifier_namespace(
            provenance=provenance,
            manifest=manifest,
        ),
        required_tables=(
            _table_report(
                table_name=_KINASE_SUBSTRATE_TABLE,
                frame=kinase_substrate_map,
                required_columns=_KINASE_SUBSTRATE_REQUIRED_COLUMNS,
                important_fields=_KINASE_SUBSTRATE_REQUIRED_COLUMNS,
            ),
            _table_report(
                table_name=_SITE_SEQUENCES_TABLE,
                frame=site_sequences,
                required_columns=_SITE_SEQUENCE_REQUIRED_COLUMNS,
                important_fields=("index", "site_sequence"),
            ),
        ),
        required_source_files=required_source_files,
        kinase_substrate_record_count=int(kinase_substrate_map.shape[0]),
        duplicate_record_count=len(duplicate_records),
        duplicate_records=duplicate_records,
        provenance_fields=provenance_fields,
        compatibility_warnings=warnings,
    )


def _table_report(
    *,
    table_name: str,
    frame: pd.DataFrame,
    required_columns: tuple[str, ...],
    important_fields: tuple[str, ...],
) -> ReferenceBundleTableValidationReport:
    present_columns = tuple(str(column) for column in frame.columns.tolist())
    missing_required_columns = tuple(
        column for column in required_columns if column not in present_columns
    )
    return ReferenceBundleTableValidationReport(
        table_name=table_name,
        present=True,
        row_count=int(frame.shape[0]),
        required_columns=required_columns,
        present_columns=present_columns,
        missing_required_columns=missing_required_columns,
        missing_values=_missing_value_counts(
            table_name=table_name,
            frame=frame,
            fields=important_fields,
        ),
    )


def _missing_value_counts(
    *,
    table_name: str,
    frame: pd.DataFrame,
    fields: tuple[str, ...],
) -> tuple[ReferenceBundleMissingValueCount, ...]:
    counts: list[ReferenceBundleMissingValueCount] = []
    for field_name in fields:
        if field_name == "index":
            values = pd.Series(frame.index.tolist(), dtype="object")
        elif field_name in frame.columns:
            values = frame[field_name]
        else:
            continue
        missing_count = int(values.isna().sum())
        if missing_count:
            counts.append(
                ReferenceBundleMissingValueCount(
                    table_name=table_name,
                    field_name=field_name,
                    missing_count=missing_count,
                )
            )
    return tuple(counts)


def _duplicate_kinase_substrate_records(
    frame: pd.DataFrame,
) -> tuple[tuple[str, str], ...]:
    kinases = frame["kinase"].astype(str).tolist()
    substrate_sites = frame["substrate_site"].astype(str).tolist()
    seen: set[tuple[str, str]] = set()
    duplicates: list[tuple[str, str]] = []
    duplicate_set: set[tuple[str, str]] = set()
    for kinase, substrate_site in zip(kinases, substrate_sites, strict=True):
        pair = (str(kinase), str(substrate_site))
        if pair in seen:
            duplicate_set.add(pair)
            if pair not in duplicates:
                duplicates.append(pair)
            continue
        seen.add(pair)
    if not duplicate_set:
        return ()
    return tuple(pair for pair in duplicates if pair in duplicate_set)


def _build_source_file_reports(
    *,
    source_files: Mapping[str, JsonValue] | None,
) -> tuple[ReferenceBundleSourceFileValidationReport, ...]:
    reports: list[ReferenceBundleSourceFileValidationReport] = []
    for role, aliases in _SOURCE_FILE_ROLE_ALIASES:
        source_file = _find_source_file(source_files=source_files, aliases=aliases)
        reports.append(
            ReferenceBundleSourceFileValidationReport(
                role=role,
                present=source_file is not None,
                path=_extract_source_file_path(source_file),
            )
        )
    return tuple(reports)


def _find_source_file(
    *,
    source_files: Mapping[str, JsonValue] | None,
    aliases: tuple[str, ...],
) -> JsonValue | None:
    if source_files is None:
        return None
    for alias in aliases:
        value = source_files.get(alias)
        if value is not None:
            return value
    return None


def _extract_source_file_path(value: JsonValue | None) -> str | None:
    if isinstance(value, str):
        return value
    if isinstance(value, Mapping):
        path = value.get("path")
        if isinstance(path, str) and path.strip():
            return path.strip()
    return None


def _resolve_source_files(
    *,
    provenance: ReferenceProvenance | None,
    manifest: ReferenceManifest | None,
) -> Mapping[str, JsonValue] | None:
    if manifest is not None:
        return manifest.source_files
    if provenance is None or provenance.manifest is None:
        return None
    raw_source_files = provenance.manifest.get("source_files")
    if isinstance(raw_source_files, Mapping):
        return {
            str(key): cast(JsonValue, value) for key, value in raw_source_files.items()
        }
    return None


def _build_provenance_fields(
    *,
    organism: Organism,
    provenance: ReferenceProvenance | None,
    manifest: ReferenceManifest | None,
) -> dict[str, JsonValue]:
    fields: dict[str, JsonValue] = {
        "source_type": provenance.source_type if provenance is not None else "explicit",
        "organism": _resolve_organism(
            organism=organism,
            provenance=provenance,
            manifest=manifest,
        ),
    }
    for key, value in (
        ("bundle_id", _resolve_bundle_name(provenance=provenance, manifest=manifest)),
        (
            "source_name",
            _resolve_source_name(provenance=provenance, manifest=manifest),
        ),
        (
            "source_version",
            _resolve_source_version(provenance=provenance, manifest=manifest),
        ),
        (
            "retrieved_at",
            _resolve_retrieved_at(provenance=provenance, manifest=manifest),
        ),
        (
            "identifier_namespace",
            _resolve_identifier_namespace(
                provenance=provenance,
                manifest=manifest,
            ),
        ),
    ):
        if value is not None:
            fields[key] = value
    reference_context = (
        None if provenance is None else provenance.reference_context
    ) or (
        reference_context_from_manifest_if_complete(manifest)
        if manifest is not None
        else None
    )
    if reference_context is not None:
        fields["reference_context_id"] = reference_context.reference_context_id
    if manifest is not None:
        fields["license"] = manifest.license
        fields["redistribution_status"] = manifest.redistribution_status.value
        fields["supports"] = manifest.supports
        fields["limitations"] = manifest.limitations
        if manifest.source_url is not None:
            fields["source_url"] = manifest.source_url
        if manifest.license_url is not None:
            fields["license_url"] = manifest.license_url
        fields["retrieval_method"] = manifest.retrieval_method
        fields["redistribution_basis"] = manifest.redistribution_basis
    fields["source_files_available"] = (
        _resolve_source_files(
            provenance=provenance,
            manifest=manifest,
        )
        is not None
    )
    return fields


def _build_compatibility_warnings(
    *,
    organism_common_name: str | None,
    identifier_namespace: str | None,
    bundle_version: str | None,
    required_source_files: tuple[ReferenceBundleSourceFileValidationReport, ...],
    manifest: ReferenceManifest | None,
    provenance: ReferenceProvenance | None,
) -> tuple[str, ...]:
    warnings: list[str] = []
    if manifest is None and (provenance is None or provenance.manifest is None):
        warnings.append(
            "reference manifest metadata is not available; provenance review is limited"
        )
    if organism_common_name is None:
        warnings.append(
            "reference organism metadata is limited to ReferenceBundle.organism"
        )
    if identifier_namespace is None:
        warnings.append("reference identifier namespace metadata is not available")
    if bundle_version is None:
        warnings.append("reference bundle version metadata is not available")
    missing_source_roles = tuple(
        item.role for item in required_source_files if not item.present
    )
    if missing_source_roles:
        roles = ", ".join(missing_source_roles)
        warnings.append(
            "reference source-file metadata is incomplete for required role(s): "
            f"{roles}"
        )
    return tuple(dict.fromkeys(warnings))


def _resolve_bundle_name(
    *,
    provenance: ReferenceProvenance | None,
    manifest: ReferenceManifest | None,
) -> str | None:
    if manifest is not None:
        return manifest.bundle_id
    if provenance is not None:
        return provenance.bundle_id
    return None


def _resolve_bundle_version(
    *,
    provenance: ReferenceProvenance | None,
    manifest: ReferenceManifest | None,
) -> str | None:
    if manifest is not None:
        return manifest.reference_version
    if provenance is not None:
        return provenance.source_version
    return None


def _resolve_source_version(
    *,
    provenance: ReferenceProvenance | None,
    manifest: ReferenceManifest | None,
) -> str | None:
    if manifest is not None:
        return manifest.source_version
    if provenance is not None:
        return provenance.source_version
    return None


def _resolve_organism(
    *,
    organism: Organism,
    provenance: ReferenceProvenance | None,
    manifest: ReferenceManifest | None,
) -> str:
    if manifest is not None:
        return manifest.organism
    if provenance is not None:
        return organism_value(provenance.organism)
    return organism_value(organism)


def _resolve_organism_common_name(
    *,
    manifest: ReferenceManifest | None,
) -> str | None:
    if manifest is None:
        return None
    return manifest.organism_common_name


def _resolve_identifier_namespace(
    *,
    provenance: ReferenceProvenance | None,
    manifest: ReferenceManifest | None,
) -> str | None:
    if manifest is not None:
        return manifest.identifier_namespace
    if provenance is not None:
        return provenance.identifier_namespace
    return None


def _resolve_source_name(
    *,
    provenance: ReferenceProvenance | None,
    manifest: ReferenceManifest | None,
) -> str | None:
    if manifest is not None:
        return manifest.source_name
    if provenance is not None:
        return provenance.source_name
    return None


def _resolve_retrieved_at(
    *,
    provenance: ReferenceProvenance | None,
    manifest: ReferenceManifest | None,
) -> str | None:
    if manifest is not None:
        return manifest.retrieved_at.isoformat()
    if provenance is not None:
        return provenance.retrieved_at
    return None


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


@dataclass(frozen=True, slots=True, init=False, eq=False)
class ReferenceBundle:
    """Resolved workflow reference resources.

    Large kinase-substrate maps are supported. Runtime in downstream workflows
    is primarily controlled by dataset/reference overlap after interpreter and
    scoring-lane filtering, not only by raw map row count.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit reference-content comparison.
    """

    __hash__ = object.__hash__

    organism: Organism
    kinase_substrate_map: pd.DataFrame
    site_sequences: pd.DataFrame
    provenance: ReferenceProvenance | None = None
    manifest: ReferenceManifest | None = None
    _validation_report: ReferenceBundleValidationReport = field(
        init=False,
        repr=False,
    )

    def __init__(
        self,
        organism: Organism,
        kinase_substrate_map: pd.DataFrame,
        site_sequences: pd.DataFrame,
        provenance: ReferenceProvenance | None = None,
        manifest: ReferenceManifest | None = None,
    ) -> None:
        self._init_reference_bundle(
            organism=organism,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            provenance=provenance,
            manifest=manifest,
            assume_owned=False,
        )

    def _init_reference_bundle(
        self,
        *,
        organism: Organism,
        kinase_substrate_map: pd.DataFrame,
        site_sequences: pd.DataFrame,
        provenance: ReferenceProvenance | None = None,
        manifest: ReferenceManifest | None = None,
        assume_owned: bool,
    ) -> None:
        if not isinstance(cast(object, organism), Organism):
            raise ReferenceValidationError(
                "references.organism must be an Organism enum value"
            )
        kinase_substrate_map = own_dataframe(
            kinase_substrate_map,
            field_name="references.kinase_substrate_map",
            error_type=ReferenceValidationError,
            assume_owned=assume_owned,
        )
        site_sequences = own_dataframe(
            site_sequences,
            field_name="references.site_sequences",
            error_type=ReferenceValidationError,
            assume_owned=assume_owned,
        )
        reference_context = (
            reference_context_from_manifest_if_complete(manifest)
            if manifest is not None
            else None
        )
        if provenance is not None:
            validate_reference_source_version_agreement(
                (
                    ("provenance.source_version", provenance.source_version),
                    (
                        "reference_context.source_version",
                        _reference_context_source_version(provenance.reference_context)
                        or _reference_context_source_version(reference_context),
                    ),
                    (
                        "manifest.source_version",
                        None if manifest is None else manifest.source_version,
                    ),
                )
            )
        validation = ReferenceBundleValidator().run(
            organism=organism,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            provenance=provenance,
            manifest=manifest,
        )
        identifier_normalisation = validation.identifier_normalisation
        if provenance is None:
            provenance = ReferenceProvenance(
                source_type="explicit",
                organism=organism,
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
                reference_context=reference_context,
            )
        else:
            resolved_reference_context = reference_context_from_provenance(
                provenance,
                manifest=manifest,
            )
            should_attach_identifier_normalisation = (
                provenance.source_type == "explicit"
                and provenance.identifier_normalisation is None
            )
            should_attach_reference_context = (
                provenance.reference_context is None
                and resolved_reference_context is not None
            )
            if (
                should_attach_identifier_normalisation
                or should_attach_reference_context
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
                    identifier_normalisation=(
                        identifier_normalisation
                        if should_attach_identifier_normalisation
                        else provenance.identifier_normalisation
                    ),
                    reference_context=(
                        resolved_reference_context
                        if should_attach_reference_context
                        else provenance.reference_context
                    ),
                )
        object.__setattr__(
            self,
            "organism",
            organism,
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
        bundle = object.__new__(cls)
        ReferenceBundle._init_reference_bundle(
            bundle,
            organism=organism,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            provenance=provenance,
            manifest=manifest,
            assume_owned=True,
        )
        return bundle

    @classmethod
    def from_trusted_owned(
        cls,
        *,
        organism: Organism,
        kinase_substrate_map: pd.DataFrame,
        site_sequences: pd.DataFrame,
        provenance: ReferenceProvenance | None = None,
        manifest: ReferenceManifest | None = None,
    ) -> ReferenceBundle:
        """Construct from already-owned tables at trusted reference boundaries."""

        return cls._from_owned(
            organism=organism,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            provenance=provenance,
            manifest=manifest,
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

    def scientifically_equals(
        self,
        other: object,
        *,
        include_provenance: bool = True,
    ) -> bool:
        """Return ``True`` when another reference bundle has the same content."""

        if not isinstance(other, ReferenceBundle):
            return False
        same_content = (
            self.organism == other.organism
            and dataframe_equals(
                self.kinase_substrate_map,
                other.kinase_substrate_map,
            )
            and dataframe_equals(self.site_sequences, other.site_sequences)
            and self._validation_report == other._validation_report
        )
        if not same_content:
            return False
        if include_provenance:
            return (
                self.provenance == other.provenance and self.manifest == other.manifest
            )
        return True
