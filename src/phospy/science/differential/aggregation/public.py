"""Public peptide-to-site differential aggregation shell."""

from __future__ import annotations

import pandas as pd

from phospy.science.differential.aggregation.executor import (
    PeptideToSiteAggregationExecutor,
)
from phospy.science.differential.aggregation.models import (
    PeptideToSiteAggregationConfig,
    PeptideToSiteAggregationResult,
)
from phospy.science.differential.models import DifferentialAnalysisResult
from phospy.science.evidence.models import PeptideEvidenceTable


class PeptideToSiteAggregator:
    """Aggregate peptide-level differential outputs to site-level summaries."""

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
