"""Experimental compatibility route for peptide-to-site differential aggregation.

This module is internal/experimental. It is intentionally not re-exported from
``phospy.api``, ``phospy.science.differential``, or the
``phospy.science.differential.aggregation`` package root.

The retained implementation performs post-hoc aggregation of peptide-level
differential statistics from the same experiment. It is not a supported
production inferential lane while the statistical model is being corrected.
"""

from __future__ import annotations

import pandas as pd

from phospy.science.differential.aggregation.executor import (
    PeptideToSiteAggregationExecutor,
)
from phospy.science.differential.aggregation.models import (
    PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS,
    PEPTIDE_TO_SITE_STRATEGY_COMPAT_BEST_P_VALUE,
    PEPTIDE_TO_SITE_STRATEGY_FIXED_EFFECT_META,
    PEPTIDE_TO_SITE_STRATEGY_INVERSE_VARIANCE_WEIGHTED,
    PEPTIDE_TO_SITE_STRATEGY_RANDOM_EFFECT_META,
    PEPTIDE_TO_SITE_STRATEGY_STOUFFER_Z,
    PeptideToSiteAggregationConfig,
    PeptideToSiteAggregationResult,
)
from phospy.science.differential.models import DifferentialAnalysisResult
from phospy.science.evidence.models import PeptideEvidenceTable

EXPERIMENTAL_INTERNAL_API = True
EXPERIMENTAL_INTERNAL_REASON = (
    "Post-hoc same-experiment peptide meta-analysis is not a supported "
    "site-level inferential lane while the statistical model is being corrected."
)


class PeptideToSiteAggregator:
    """Experimental/internal compatibility shell for post-hoc site summaries."""

    experimental_internal_api: bool = True
    scientific_support_status: str = PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS

    def __init__(
        self,
        *,
        executor: PeptideToSiteAggregationExecutor | None = None,
    ) -> None:
        self._executor = executor or PeptideToSiteAggregationExecutor()

    def run_table(
        self,
        *,
        peptide_differential_table: pd.DataFrame,
        evidence: PeptideEvidenceTable,
        config: PeptideToSiteAggregationConfig | None = None,
        contrast_name: str = "aggregated",
    ) -> PeptideToSiteAggregationResult:
        resolved_config = config or PeptideToSiteAggregationConfig()
        return self._executor.run_table(
            peptide_differential_table=peptide_differential_table,
            evidence=evidence,
            config=resolved_config,
            contrast_name=contrast_name,
        )

    def run_differential_result(
        self,
        *,
        differential_result: DifferentialAnalysisResult,
        evidence: PeptideEvidenceTable,
        config: PeptideToSiteAggregationConfig | None = None,
    ) -> dict[str, PeptideToSiteAggregationResult]:
        resolved_config = config or PeptideToSiteAggregationConfig()
        return self._executor.run_differential_result(
            differential_result=differential_result,
            evidence=evidence,
            config=resolved_config,
        )


__all__ = [
    "EXPERIMENTAL_INTERNAL_API",
    "EXPERIMENTAL_INTERNAL_REASON",
    "PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS",
    "PEPTIDE_TO_SITE_STRATEGY_COMPAT_BEST_P_VALUE",
    "PEPTIDE_TO_SITE_STRATEGY_FIXED_EFFECT_META",
    "PEPTIDE_TO_SITE_STRATEGY_INVERSE_VARIANCE_WEIGHTED",
    "PEPTIDE_TO_SITE_STRATEGY_RANDOM_EFFECT_META",
    "PEPTIDE_TO_SITE_STRATEGY_STOUFFER_Z",
    "PeptideToSiteAggregationConfig",
    "PeptideToSiteAggregationResult",
    "PeptideToSiteAggregator",
]
