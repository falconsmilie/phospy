"""Reference bundle validator."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import cast

import pandas as pd

from phospy.errors.validation import ReferenceValidationError
from phospy.provenance.models import JsonValue, ReferenceProvenance
from phospy.science.references.identifiers import (
    ReferenceIdentifierNormalisationReport,
    merge_reference_identifier_normalisation_reports,
)
from phospy.science.references.models import (
    Organism,
    ReferenceBundleMissingValueCount,
    ReferenceBundleSourceFileValidationReport,
    ReferenceBundleTableValidationReport,
    ReferenceBundleValidationReport,
    ReferenceManifest,
    reference_context_from_manifest_if_complete,
)
from phospy.tables.references import KinaseSubstrateReference, SiteSequenceReference

_KINASE_SUBSTRATE_TABLE = "kinase_substrate_map"
_SITE_SEQUENCES_TABLE = "site_sequences"
_KINASE_SUBSTRATE_REQUIRED_COLUMNS = ("kinase", "substrate_site")
_SITE_SEQUENCE_REQUIRED_COLUMNS = ("site_sequence",)
_SOURCE_FILE_ROLE_ALIASES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("kinase_substrate", ("kinase_substrate", "kinase_substrate_map", "substrate_map")),
    ("site_sequences", ("site_sequences", "site_sequence", "site_sequence_table")),
)


@dataclass(frozen=True, slots=True)
class ReferenceBundleValidationResult:
    """Validated reference tables plus the public validation report."""

    kinase_substrate_map: pd.DataFrame
    site_sequences: pd.DataFrame
    identifier_normalisation: ReferenceIdentifierNormalisationReport | None
    report: ReferenceBundleValidationReport


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
            _resolve_identifier_namespace(provenance=provenance, manifest=manifest),
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
        return provenance.organism
    return organism.value


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
