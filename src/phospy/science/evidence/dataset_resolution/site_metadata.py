"""Site metadata and localisation aggregation for resolved peptide evidence."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.evidence.dataset_resolution.allocation import AllocatedEvidence
from phospy.science.evidence.dataset_resolution.contracts import (
    DATASET_PEPTIDE_LOCALISATION_COMPATIBILITY_ALIAS_COLUMN,
    DATASET_PEPTIDE_LOCALISATION_SUMMARY_COLUMN,
    DATASET_PEPTIDE_LOCALISATION_SUMMARY_SEMANTICS,
    DATASET_PEPTIDE_LOCALISATION_SUMMARY_SEMANTICS_COLUMN,
    SITE_SEQUENCE_SOURCE_MISSING,
    SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT,
    SITE_SEQUENCE_SOURCE_PROVIDED,
)
from phospy.science.evidence.dataset_resolution.site_sequence import (
    SiteSequenceResolutionDiagnostics,
    non_empty_strings,
    resolve_site_sequence_for_resolved_site,
)
from phospy.science.sites.identifiers import parse_canonical_site_identifier


@dataclass(frozen=True, slots=True)
class SiteMetadataResolution:
    """Site metadata and sequence diagnostics produced from allocated evidence."""

    site_metadata: pd.DataFrame
    sequence_diagnostics: SiteSequenceResolutionDiagnostics


def aggregate_site_metadata_and_localisation(
    *,
    allocated_evidence: AllocatedEvidence,
    site_ids: pd.Index,
    multi_site_policy: str,
) -> SiteMetadataResolution:
    """Aggregate resolved-site metadata and localisation confidence."""

    mapped_rows = allocated_evidence.rows
    grouped = mapped_rows.groupby("site_id", sort=True)
    include_localisation_confidence = "localisation_confidence" in mapped_rows.columns
    site_rows: list[dict[str, object]] = []
    provided_site_sequence_used_count = 0
    peptide_context_derived_site_sequence_count = 0
    missing_site_sequence_count = 0
    rejected_provided_context_count = 0
    for site_id in site_ids.astype(str).tolist():
        group = grouped.get_group(site_id)
        gene_symbol, site = parse_canonical_site_identifier(
            site_id,
            field_name="dataset peptide evidence site_id",
            error_type=PhosPyInputError,
        )
        protein_accession = single_non_empty_string_or_error(
            group.loc[:, "protein_accession"],
            field_name="protein_accession",
            site_id=site_id,
        )
        site_sequence_resolution = resolve_site_sequence_for_resolved_site(
            group=group,
            site_id=site_id,
            resolved_site_token=site,
            multi_site_policy=multi_site_policy,
        )
        if site_sequence_resolution.source == SITE_SEQUENCE_SOURCE_PROVIDED:
            provided_site_sequence_used_count += 1
        elif site_sequence_resolution.source == SITE_SEQUENCE_SOURCE_PEPTIDE_CONTEXT:
            peptide_context_derived_site_sequence_count += 1
        elif site_sequence_resolution.source == SITE_SEQUENCE_SOURCE_MISSING:
            missing_site_sequence_count += 1
        rejected_provided_context_count += (
            site_sequence_resolution.rejected_provided_context_count
        )
        site_rows.append(
            {
                "site_id": site_id,
                "gene_symbol": gene_symbol,
                "site": site,
                "site_sequence": site_sequence_resolution.site_sequence,
                "protein_accession": protein_accession,
                "protein_namespace": (
                    "protein_accession" if protein_accession is not None else None
                ),
                "protein_identifier": protein_accession,
            }
        )
        if include_localisation_confidence:
            localisation_summary = aggregate_localisation_confidence(
                group.loc[:, "localisation_confidence"]
            )
            site_rows[-1][DATASET_PEPTIDE_LOCALISATION_SUMMARY_COLUMN] = (
                localisation_summary
            )
            site_rows[-1][DATASET_PEPTIDE_LOCALISATION_SUMMARY_SEMANTICS_COLUMN] = (
                DATASET_PEPTIDE_LOCALISATION_SUMMARY_SEMANTICS
            )
            site_rows[-1][DATASET_PEPTIDE_LOCALISATION_COMPATIBILITY_ALIAS_COLUMN] = (
                localisation_summary
            )
    site_metadata = pd.DataFrame(site_rows).set_index("site_id", drop=True)
    site_metadata.index = pd.Index(site_metadata.index.astype(str), name="site_id")
    return SiteMetadataResolution(
        site_metadata=site_metadata,
        sequence_diagnostics=SiteSequenceResolutionDiagnostics(
            rejected_provided_context_count=int(rejected_provided_context_count),
            provided_site_sequence_used_count=int(provided_site_sequence_used_count),
            peptide_context_derived_site_sequence_count=int(
                peptide_context_derived_site_sequence_count
            ),
            missing_site_sequence_count=int(missing_site_sequence_count),
        ),
    )


def aggregate_localisation_confidence(values: pd.Series) -> float | None:
    numeric = pd.to_numeric(values, errors="coerce")
    finite = numeric.loc[numeric.notna()]
    if finite.empty:
        return None
    return float(finite.mean())


def single_non_empty_string_or_error(
    values: pd.Series,
    *,
    field_name: str,
    site_id: str,
) -> str | None:
    distinct = tuple(dict.fromkeys(non_empty_strings(values)))
    if len(distinct) <= 1:
        return distinct[0] if distinct else None
    preview = ", ".join(repr(value) for value in distinct[:5])
    suffix = "" if len(distinct) <= 5 else " ..."
    raise PhosPyInputError(
        f"{field_name} must contain at most one distinct non-empty value per "
        f"resolved site_id; site_id={site_id!r}, conflicting_values=["
        f"{preview}{suffix}]. Suggested fix: disambiguate peptide-site mapping "
        "or split rows before building."
    )


__all__ = [
    "SiteMetadataResolution",
    "aggregate_localisation_confidence",
    "aggregate_site_metadata_and_localisation",
    "single_non_empty_string_or_error",
]
