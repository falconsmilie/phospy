"""Compose user-facing kinase workflow eligibility diagnostics."""

from __future__ import annotations

import pandas as pd

from phospy.api.results import KinaseEligibilityReport
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_LOCALISATION,
)
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)

_KINASE_COLUMN = "kinase"
_SUBSTRATE_COLUMN = "substrate_site"


class KinaseEligibilityReportComposer:
    """Build compact kinase workflow eligibility counts for result diagnostics."""

    def run(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        config: ResolvedKinaseExecutionConfig,
    ) -> KinaseEligibilityReport:
        dataset_site_index = request.dataset.phospho.index.astype(str)
        dataset_sites = set(dataset_site_index.tolist())
        total_dataset_sites = int(len(dataset_sites))
        sequence_complete_sites = int(len(request.scoring_site_index))

        localisation_eligible_sites, excluded_low_localisation = (
            self._resolve_localisation_counts(
                request,
                final_dataset_sites=total_dataset_sites,
            )
        )

        substrate_series = request.kinase_substrate_map.loc[
            :, _SUBSTRATE_COLUMN
        ].astype(str)
        reference_sites = set(substrate_series.tolist())
        sequence_complete_site_ids = set(
            request.scoring_site_index.astype(str).tolist()
        )
        reference_overlap_site_ids = sequence_complete_site_ids.intersection(
            reference_sites
        )
        reference_overlap_sites = int(len(reference_overlap_site_ids))
        excluded_no_reference_match = max(
            sequence_complete_sites - reference_overlap_sites,
            0,
        )

        overlap_map = request.kinase_substrate_map.loc[
            substrate_series.isin(reference_overlap_site_ids),
            [_KINASE_COLUMN, _SUBSTRATE_COLUMN],
        ]
        if overlap_map.empty:
            per_kinase_quantified = pd.Series(dtype="int64")
        else:
            per_kinase_quantified = (
                overlap_map.groupby(_KINASE_COLUMN, sort=False)[_SUBSTRATE_COLUMN]
                .nunique()
                .astype("int64")
            )
        eligible_mask = per_kinase_quantified >= int(config.scoring_min_substrates)
        eligible_kinases = int(eligible_mask.sum())
        excluded_kinases_below_min_substrates = int(
            per_kinase_quantified.size - eligible_kinases
        )

        return KinaseEligibilityReport(
            total_dataset_sites=total_dataset_sites,
            sequence_complete_sites=sequence_complete_sites,
            localisation_eligible_sites=localisation_eligible_sites,
            reference_overlap_sites=reference_overlap_sites,
            excluded_no_reference_match=excluded_no_reference_match,
            excluded_low_localisation=excluded_low_localisation,
            eligible_kinases=eligible_kinases,
            excluded_kinases_below_min_substrates=(
                excluded_kinases_below_min_substrates
            ),
        )

    @staticmethod
    def _resolve_localisation_counts(
        request: ResolvedKinaseWorkflowRequest,
        *,
        final_dataset_sites: int,
    ) -> tuple[int | None, int | None]:
        report = request.dataset.preprocessing_report
        if report is None:
            return (None, None)
        operations = report.operations
        if operations.empty or "stage" not in operations.columns:
            return (None, None)
        stage_mask = (
            operations.loc[:, "stage"].astype(str)
            == DATASET_PREPROCESSING_STAGE_LOCALISATION
        )
        if not bool(stage_mask.any()):
            return (None, None)
        stage_row = operations.loc[stage_mask, :].iloc[-1]
        input_rows = int(stage_row.loc["input_rows"])
        output_rows = int(stage_row.loc["output_rows"])
        return (min(output_rows, final_dataset_sites), max(input_rows - output_rows, 0))


__all__ = ["KinaseEligibilityReportComposer"]
