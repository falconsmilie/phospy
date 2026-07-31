"""Backward-compatible import route for peptide-to-site estimate combination.

The implementation is no longer the old experimental raw-table combiner. It
delegates to the supported typed estimate model exposed from
``phospy.science.differential.aggregation``.
"""

from __future__ import annotations

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.differential.aggregation.executor import (
    PeptideToSiteAggregationExecutor,
)
from phospy.science.differential.aggregation.models import (
    PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS,
    PeptideDifferentialEstimateTable,
    PeptideToSiteAggregationConfig,
    PeptideToSiteAggregationResult,
)

EXPERIMENTAL_INTERNAL_API = False
EXPERIMENTAL_INTERNAL_REASON = (
    "The old raw-table post-hoc combiner was replaced by a supported typed "
    "estimate-combination contract with estimate identity and moderated-t "
    "row-consistency validation. Use phospy.science.differential.aggregation "
    "for supported advanced code."
)


class PeptideToSiteAggregator:
    """Supported shell for typed post-hoc peptide-to-site estimate combination."""

    experimental_internal_api: bool = False
    scientific_support_status: str = PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS

    def __init__(
        self,
        *,
        executor: PeptideToSiteAggregationExecutor | None = None,
    ) -> None:
        self._executor = executor or PeptideToSiteAggregationExecutor()

    def run(
        self,
        estimates: PeptideDifferentialEstimateTable | pd.DataFrame,
        *,
        config: PeptideToSiteAggregationConfig | None = None,
        contrast_name: str = "aggregated",
    ) -> PeptideToSiteAggregationResult:
        resolved_estimates = _coerce_estimates(estimates)
        resolved_config = config or PeptideToSiteAggregationConfig()
        return self._executor.run_estimates(
            estimates=resolved_estimates,
            config=resolved_config,
            contrast_name=contrast_name,
        )

    def run_estimates(
        self,
        *,
        estimates: PeptideDifferentialEstimateTable | pd.DataFrame,
        config: PeptideToSiteAggregationConfig | None = None,
        contrast_name: str = "aggregated",
    ) -> PeptideToSiteAggregationResult:
        return self.run(
            estimates,
            config=config,
            contrast_name=contrast_name,
        )

    def run_table(
        self,
        *,
        estimate_table: PeptideDifferentialEstimateTable | pd.DataFrame | None = None,
        estimates: PeptideDifferentialEstimateTable | pd.DataFrame | None = None,
        peptide_differential_table: pd.DataFrame | None = None,
        evidence: object | None = None,
        config: PeptideToSiteAggregationConfig | None = None,
        contrast_name: str = "aggregated",
    ) -> PeptideToSiteAggregationResult:
        if peptide_differential_table is not None or evidence is not None:
            raise PhosPyInputError(
                "raw peptide_differential_table + evidence aggregation is no "
                "longer supported because it would infer uncertainty from "
                "logFC/t. Build a PeptideDifferentialEstimateTable with typed "
                "standard_error, degrees of freedom, source_experiment_id, "
                "dependence_policy, and peptide_to_site_mapping_policy."
            )
        resolved = estimate_table if estimate_table is not None else estimates
        if resolved is None:
            raise PhosPyInputError(
                "peptide-to-site aggregation requires estimate_table or estimates"
            )
        return self.run(
            resolved,
            config=config,
            contrast_name=contrast_name,
        )


def _coerce_estimates(
    value: PeptideDifferentialEstimateTable | pd.DataFrame,
) -> PeptideDifferentialEstimateTable:
    if isinstance(value, PeptideDifferentialEstimateTable):
        return value
    if isinstance(value, pd.DataFrame):
        return PeptideDifferentialEstimateTable(value)
    raise PhosPyInputError(
        "peptide-to-site aggregation estimates must be a "
        "PeptideDifferentialEstimateTable or pandas DataFrame"
    )


__all__ = [
    "EXPERIMENTAL_INTERNAL_API",
    "EXPERIMENTAL_INTERNAL_REASON",
    "PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS",
    "PeptideDifferentialEstimateTable",
    "PeptideToSiteAggregationConfig",
    "PeptideToSiteAggregationResult",
    "PeptideToSiteAggregator",
]
