"""Coordinator for peptide-evidence dataset resolution."""

from __future__ import annotations

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.evidence.dataset_resolution.allocation import (
    allocate_peptide_signals_to_resolved_sites,
    summarise_allocated_site_signals,
)
from phospy.science.evidence.dataset_resolution.contracts import (
    DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING,
    DATASET_MULTI_SITE_POLICY_SPLIT,
    PeptideEvidenceResolutionInputMetrics,
    PeptideEvidenceResolutionResult,
    validate_dataset_multi_site_policy,
)
from phospy.science.evidence.dataset_resolution.mapping import (
    join_peptide_rows_to_site_mapping,
    resolve_and_validate_mapping_fractions,
)
from phospy.science.evidence.dataset_resolution.site_metadata import (
    aggregate_site_metadata_and_localisation,
)
from phospy.science.evidence.dataset_resolution.site_sequence import (
    count_non_empty_strings,
)
from phospy.science.evidence.dataset_resolution.summary import (
    build_resolution_summary,
)
from phospy.science.evidence.models import PeptideEvidenceTable


class PeptideEvidenceDatasetResolver:
    """Resolve peptide-level evidence into site-level dataset-builder tables."""

    def run(
        self,
        *,
        evidence: PeptideEvidenceTable,
        multi_site_policy: str,
    ) -> PeptideEvidenceResolutionResult:
        if not isinstance(evidence, PeptideEvidenceTable):
            raise PhosPyInputError(
                "dataset peptide evidence resolution requires a PeptideEvidenceTable"
            )
        validate_dataset_multi_site_policy(
            multi_site_policy,
            field_name="dataset build request multi_site_policy",
        )
        evidence_frame = evidence.to_dataframe()
        mapping = evidence.site_mapping.to_dataframe()
        input_metrics = _collect_input_metrics(
            evidence_frame=evidence_frame,
            multi_site_policy=multi_site_policy,
        )

        if mapping.empty:
            raise PhosPyInputError(
                "dataset build request peptide_evidence resolved to zero mapped "
                "site rows after applying multi_site_policy"
            )
        joined_mapping = join_peptide_rows_to_site_mapping(
            evidence_frame=evidence_frame,
            mapping=mapping,
            sample_columns=evidence.sample_intensity_columns,
        )
        if joined_mapping.rows.empty:
            raise PhosPyInputError(
                "dataset build request peptide_evidence resolved to zero mapped "
                "site rows after joining peptide evidence and site mapping"
            )
        resolved_mapping = resolve_and_validate_mapping_fractions(
            joined_mapping=joined_mapping
        )
        allocated_evidence = allocate_peptide_signals_to_resolved_sites(
            resolved_mapping=resolved_mapping,
            sample_columns=evidence.sample_intensity_columns,
        )
        site_signals = summarise_allocated_site_signals(
            allocated_evidence=allocated_evidence
        )
        site_metadata_resolution = aggregate_site_metadata_and_localisation(
            allocated_evidence=allocated_evidence,
            site_ids=site_signals.phospho.index,
            multi_site_policy=multi_site_policy,
        )
        summary = build_resolution_summary(
            multi_site_policy=multi_site_policy,
            input_metrics=input_metrics,
            resolved_mapping=resolved_mapping,
            site_signals=site_signals,
            site_metadata_resolution=site_metadata_resolution,
        )
        return PeptideEvidenceResolutionResult(
            phospho=site_signals.phospho,
            site_metadata=site_metadata_resolution.site_metadata,
            summary=summary,
        )


def _collect_input_metrics(
    *,
    evidence_frame: pd.DataFrame,
    multi_site_policy: str,
) -> PeptideEvidenceResolutionInputMetrics:
    site_sequence_column_present = "site_sequence" in evidence_frame.columns
    provided_site_sequence_count = (
        count_non_empty_strings(evidence_frame.loc[:, "site_sequence"])
        if site_sequence_column_present
        else 0
    )
    peptide_observations_received = int(evidence_frame.shape[0])
    ambiguous_observations = int(evidence_frame.loc[:, "multi_site"].astype(bool).sum())
    excluded_observations = (
        ambiguous_observations
        if multi_site_policy == DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING
        else 0
    )
    split_observations = (
        ambiguous_observations
        if multi_site_policy == DATASET_MULTI_SITE_POLICY_SPLIT
        else 0
    )
    duplicate_peptide_rows = int(
        evidence_frame.loc[:, "peptide_sequence"]
        .astype(str)
        .duplicated(keep=False)
        .sum()
    )
    return PeptideEvidenceResolutionInputMetrics(
        peptide_observations_received=peptide_observations_received,
        ambiguous_observations=ambiguous_observations,
        excluded_observations=excluded_observations,
        split_observations=split_observations,
        duplicate_peptide_rows=duplicate_peptide_rows,
        site_sequence_column_present=site_sequence_column_present,
        provided_site_sequence_count=provided_site_sequence_count,
    )


__all__ = [
    "PeptideEvidenceDatasetResolver",
]
