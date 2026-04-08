from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from .activities import (
    build_kinase_target_table,
    compute_ksea_scores,
    compute_weighted_kinase_activity,
    count_predicted_targets,
)
from .io import load_pred_mat
from .validation.analysis import (
    ValidatedAnalysisRequest,
    validate_analysis_request,
)
from .writers import KinaseActivityResultWriter, KinaseActivityWriter

__all__ = ["KinaseActivityAnalyzer", "KinaseActivityResult"]


if TYPE_CHECKING:
    from .prediction.models import PredMatResult


PredMatLoader = Callable[[str | Path], pd.DataFrame]


@dataclass(slots=True)
class KinaseActivityResult:
    """Kinase activity tables produced by one analyzer run."""

    weighted_activity: pd.DataFrame
    ksea_scores: pd.DataFrame
    ksea_counts: pd.Series
    target_counts: pd.Series
    target_table: pd.DataFrame


class _ActivityRunner:
    def execute(
        self,
        request: ValidatedAnalysisRequest,
    ) -> KinaseActivityResult:
        weighted_activity = compute_weighted_kinase_activity(
            pred_mat=request.pred_mat,
            phospho_matrix=request.phospho_matrix,
            top_n_substrates=request.request.top_n_substrates,
            min_substrates=request.request.min_substrates,
        )
        ksea_scores, ksea_counts = compute_ksea_scores(
            pred_mat=request.pred_mat,
            phospho_matrix=request.phospho_matrix,
            threshold=request.request.threshold,
            min_substrates=request.request.min_substrates,
        )
        target_counts = count_predicted_targets(
            request.pred_mat,
            threshold=request.request.threshold,
        )
        target_table = build_kinase_target_table(
            request.pred_mat,
            threshold=request.request.threshold,
        )

        return KinaseActivityResult(
            weighted_activity=weighted_activity,
            ksea_scores=ksea_scores,
            ksea_counts=ksea_counts,
            target_counts=target_counts,
            target_table=target_table,
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
        threshold: float = 0.6,
        min_substrates: int = 3,
        top_n_substrates: int = 20,
    ) -> ValidatedAnalysisRequest:
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
        threshold: float = 0.6,
        min_substrates: int = 3,
        top_n_substrates: int = 20,
    ) -> KinaseActivityResult:
        """Compute downstream kinase summaries from raw public inputs."""

        request = self._validate_request(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=threshold,
            min_substrates=min_substrates,
            top_n_substrates=top_n_substrates,
        )
        return self._run_request(request)

    def _run_request(
        self,
        request: ValidatedAnalysisRequest,
    ) -> KinaseActivityResult:
        return self.runner.execute(request)

    def write_outputs(self, result: KinaseActivityResult, outdir: str | Path) -> None:
        self.result_writer.write(result, outdir)
