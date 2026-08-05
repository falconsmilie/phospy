"""Summary assembly for peptide-evidence dataset resolution."""

from __future__ import annotations

from collections.abc import Sequence
from typing import cast

from phospy.science.evidence.dataset_resolution.allocation import SiteSignalSummary
from phospy.science.evidence.dataset_resolution.contracts import (
    DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR,
    DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
    MAPPING_FRACTION_COLUMN,
    MAPPING_WEIGHT_SUM_TOLERANCE,
    PeptideEvidenceResolutionInputMetrics,
    PeptideEvidenceResolutionSummary,
    PeptideToSiteAggregationPolicy,
)
from phospy.science.evidence.dataset_resolution.mapping import ResolvedMappingFractions
from phospy.science.evidence.dataset_resolution.site_metadata import (
    SiteMetadataResolution,
)
from phospy.science.evidence.dataset_resolution.site_sequence import (
    count_non_empty_strings,
)


def build_resolution_summary(
    *,
    multi_site_policy: str,
    input_metrics: PeptideEvidenceResolutionInputMetrics,
    resolved_mapping: ResolvedMappingFractions,
    aggregation_policy: PeptideToSiteAggregationPolicy,
    site_signals: SiteSignalSummary,
    site_metadata_resolution: SiteMetadataResolution,
) -> PeptideEvidenceResolutionSummary:
    """Assemble count and policy provenance without performing resolution logic."""

    policy_payload = aggregation_policy.to_payload()
    accepted_site_sequence_count = count_non_empty_strings(
        site_metadata_resolution.site_metadata.loc[:, "site_sequence"]
    )
    sequence_diagnostics = site_metadata_resolution.sequence_diagnostics
    mapping_row_count = int(resolved_mapping.rows.shape[0])
    mapped_peptide_observations = int(
        resolved_mapping.rows.loc[:, "peptide_row_id"].astype(str).nunique()
    )
    mapping_fractions = resolved_mapping.rows.loc[:, MAPPING_FRACTION_COLUMN].astype(
        float
    )
    fractional_mapping_rows = int(
        mapping_fractions.sub(1.0).abs().gt(MAPPING_WEIGHT_SUM_TOLERANCE).sum()
    )
    return PeptideEvidenceResolutionSummary(
        input_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
        multi_site_policy=multi_site_policy,
        peptide_to_site_aggregation_policy_id=str(
            policy_payload["peptide_to_site_aggregation_policy_id"]
        ),
        supported_input_scales=_payload_string_tuple(
            policy_payload["supported_input_scales"]
        ),
        supported_input_quantitative_meanings=_payload_string_tuple(
            policy_payload["supported_input_quantitative_meanings"]
        ),
        input_intensity_scale=str(policy_payload["input_intensity_scale"]),
        input_quantitative_meaning=str(policy_payload["input_quantitative_meaning"]),
        output_intensity_scale=str(policy_payload["output_intensity_scale"]),
        output_quantitative_meaning=str(policy_payload["output_quantitative_meaning"]),
        allocation_domain=str(policy_payload["allocation_domain"]),
        fractional_mapping_present=bool(policy_payload["fractional_mapping_present"]),
        peptide_observations_received=(input_metrics.peptide_observations_received),
        mapped_peptide_observations=mapped_peptide_observations,
        site_mapping_rows=mapping_row_count,
        allocated_evidence_rows=mapping_row_count,
        unique_site_ids_produced=int(site_signals.phospho.shape[0]),
        ambiguous_observations=input_metrics.ambiguous_observations,
        unambiguous_observations=input_metrics.unambiguous_observations,
        excluded_observations=input_metrics.excluded_observations,
        split_observations=input_metrics.split_observations,
        fractional_mapping_rows=fractional_mapping_rows,
        unit_mapping_rows=int(mapping_row_count - fractional_mapping_rows),
        mapping_weight_source_policy=str(
            policy_payload["mapping_weight_source_policy"]
        ),
        mapping_weight_normalization_policy=(
            str(policy_payload["mapping_weight_normalization_policy"])
        ),
        mapping_weight_semantics=str(policy_payload["mapping_weight_semantics"]),
        signal_allocation_policy=str(policy_payload["signal_allocation_policy"]),
        site_summarisation_policy=str(policy_payload["site_summarisation_policy"]),
        missing_value_policy=str(policy_payload["missing_value_policy"]),
        duplicate_evidence_policy=str(policy_payload["duplicate_evidence_policy"]),
        mixed_ambiguity_policy=str(policy_payload["mixed_ambiguity_policy"]),
        localisation_aggregation_policy=(
            str(policy_payload["localisation_aggregation_policy"])
        ),
        localisation_summary_policy=str(policy_payload["localisation_summary_policy"]),
        localisation_summary_semantics=str(
            policy_payload["localisation_summary_semantics"]
        ),
        localisation_output_column=str(policy_payload["localisation_output_column"]),
        localisation_compatibility_alias_column=str(
            policy_payload["localisation_compatibility_alias_column"]
        ),
        signal_conservation_policy=str(policy_payload["signal_conservation_policy"]),
        uncertainty_limitations=_payload_string_tuple(
            policy_payload["uncertainty_limitations"]
        ),
        aggregation_policy=str(policy_payload["aggregation_policy"]),
        aggregation_formula=str(policy_payload["aggregation_formula"]),
        mapping_weight_source=resolved_mapping.mapping_weight_source,
        mapping_weight_normalisation=(
            str(policy_payload["mapping_weight_normalization_policy"])
        ),
        duplicate_peptide_policy=str(policy_payload["duplicate_evidence_policy"]),
        duplicate_peptide_rows=input_metrics.duplicate_peptide_rows,
        site_sequence_column_present=input_metrics.site_sequence_column_present,
        provided_site_sequence_count=input_metrics.provided_site_sequence_count,
        accepted_site_sequence_count=accepted_site_sequence_count,
        rejected_site_sequence_count=(
            sequence_diagnostics.rejected_provided_context_count
        ),
        provided_site_sequence_used_count=(
            sequence_diagnostics.provided_site_sequence_used_count
        ),
        peptide_context_derived_site_sequence_count=(
            sequence_diagnostics.peptide_context_derived_site_sequence_count
        ),
        missing_site_sequence_count=sequence_diagnostics.missing_site_sequence_count,
        site_sequence_policy=(
            DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR
        ),
    )


__all__ = [
    "build_resolution_summary",
]


def _payload_string_tuple(value: object) -> tuple[str, ...]:
    if isinstance(value, str):
        return (value,)
    if isinstance(value, Sequence):
        return tuple(str(item) for item in cast(Sequence[object], value))
    return (str(value),)
