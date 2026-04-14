from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from ..internal.defaults import (
    DEFAULT_KINASE_ACTIVITY_MIN_SUBSTRATES,
    DEFAULT_KINASE_ACTIVITY_THRESHOLD,
    DEFAULT_KINASE_ACTIVITY_TOP_N_SUBSTRATES,
)
from ..io import load_pred_mat
from ..io.writers import KinaseActivityResultWriter, KinaseActivityWriter
from ..validation.requests.analysis import AnalysisInputs, validate_analysis_request
from .results import KinaseActivityResult
from .scoring import (
    build_kinase_target_table,
    compute_activity_from_inputs,
    count_predicted_targets,
)

__all__ = ["KinaseActivityAnalyzer"]


if TYPE_CHECKING:
    from ..prediction.results import PredMatResult


PredMatLoader = Callable[[str | Path], pd.DataFrame]


class _ActivityRunner:
    def execute(
        self,
        inputs: AnalysisInputs,
    ) -> KinaseActivityResult:
        weighted_activity, ksea_scores, ksea_counts = compute_activity_from_inputs(
            inputs
        )
        target_counts = count_predicted_targets(
            inputs.pred_mat,
            threshold=inputs.threshold,
        )
        target_table = build_kinase_target_table(
            inputs.pred_mat,
            threshold=inputs.threshold,
        )

        return KinaseActivityResult(
            weighted_activity=weighted_activity,
            ksea_scores=ksea_scores,
            ksea_counts=ksea_counts,
            target_counts=target_counts,
            target_table=target_table,
            overlap_summary=inputs.overlap_summary,
        )


@dataclass(slots=True)
class KinaseActivityAnalyzer:
    """Run downstream kinase activity analysis from validated tabular inputs."""

    pred_mat_loader: PredMatLoader = load_pred_mat
    result_writer: KinaseActivityResultWriter = field(
        default_factory=KinaseActivityWriter
    )
    runner: _ActivityRunner = field(default_factory=_ActivityRunner)

    def load_pred_mat(self, pred_mat_path: str | Path) -> pd.DataFrame:
        """Load and validate a kinase prediction matrix from disk."""

        return self.pred_mat_loader(pred_mat_path)

    def _validate_request(
        self,
        *,
        pred_mat: pd.DataFrame | PredMatResult,
        phospho_matrix: pd.DataFrame,
        threshold: float = DEFAULT_KINASE_ACTIVITY_THRESHOLD,
        min_substrates: int = DEFAULT_KINASE_ACTIVITY_MIN_SUBSTRATES,
        top_n_substrates: int = DEFAULT_KINASE_ACTIVITY_TOP_N_SUBSTRATES,
    ) -> AnalysisInputs:
        return validate_analysis_request(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=threshold,
            min_substrates=min_substrates,
            top_n_substrates=top_n_substrates,
        )

    def run(
        self,
        pred_mat: pd.DataFrame | PredMatResult,
        phospho_matrix: pd.DataFrame,
        threshold: float = DEFAULT_KINASE_ACTIVITY_THRESHOLD,
        min_substrates: int = DEFAULT_KINASE_ACTIVITY_MIN_SUBSTRATES,
        top_n_substrates: int = DEFAULT_KINASE_ACTIVITY_TOP_N_SUBSTRATES,
    ) -> KinaseActivityResult:
        """Compute downstream kinase summaries from raw public inputs."""

        inputs = self._validate_request(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=threshold,
            min_substrates=min_substrates,
            top_n_substrates=top_n_substrates,
        )
        return self.run_validated(inputs)

    def run_validated(
        self,
        inputs: AnalysisInputs,
    ) -> KinaseActivityResult:
        return self.runner.execute(inputs)

    def write_outputs(self, result: KinaseActivityResult, outdir: str | Path) -> None:
        self.result_writer.write(result, outdir)
