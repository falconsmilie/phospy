"""Summary assembly for peptide-evidence dataset resolution."""

from __future__ import annotations

from phospy.science.evidence.dataset_resolution.allocation import SiteSignalSummary
from phospy.science.evidence.dataset_resolution.contracts import (
    CURRENT_RESOLUTION_POLICIES,
    DATASET_PEPTIDE_SITE_SEQUENCE_POLICY_VALIDATE_WITHOUT_REPAIR,
    DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
    PeptideEvidenceResolutionInputMetrics,
    PeptideEvidenceResolutionSummary,
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
    site_signals: SiteSignalSummary,
    site_metadata_resolution: SiteMetadataResolution,
) -> PeptideEvidenceResolutionSummary:
    """Assemble count and policy provenance without performing resolution logic."""

    policy_payload = CURRENT_RESOLUTION_POLICIES.to_payload()
    accepted_site_sequence_count = count_non_empty_strings(
        site_metadata_resolution.site_metadata.loc[:, "site_sequence"]
    )
    sequence_diagnostics = site_metadata_resolution.sequence_diagnostics
    return PeptideEvidenceResolutionSummary(
        input_mode=DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE,
        multi_site_policy=multi_site_policy,
        peptide_observations_received=(input_metrics.peptide_observations_received),
        unique_site_ids_produced=int(site_signals.phospho.shape[0]),
        ambiguous_observations=input_metrics.ambiguous_observations,
        excluded_observations=input_metrics.excluded_observations,
        split_observations=input_metrics.split_observations,
        mapping_weight_source_policy=policy_payload["mapping_weight_source_policy"],
        mapping_weight_normalization_policy=(
            policy_payload["mapping_weight_normalization_policy"]
        ),
        signal_allocation_policy=policy_payload["signal_allocation_policy"],
        site_summarisation_policy=policy_payload["site_summarisation_policy"],
        duplicate_evidence_policy=policy_payload["duplicate_evidence_policy"],
        mixed_ambiguity_policy=policy_payload["mixed_ambiguity_policy"],
        localisation_aggregation_policy=(
            policy_payload["localisation_aggregation_policy"]
        ),
        aggregation_policy=policy_payload["aggregation_policy"],
        aggregation_formula=policy_payload["aggregation_formula"],
        mapping_weight_source=resolved_mapping.mapping_weight_source,
        mapping_weight_normalisation=(
            policy_payload["mapping_weight_normalization_policy"]
        ),
        duplicate_peptide_policy=policy_payload["duplicate_evidence_policy"],
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
