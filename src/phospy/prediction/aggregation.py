from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

from ..validation.prediction import PredictionRequest
from .models import KinasePredictionResult

if TYPE_CHECKING:
    from .execution import KinasePredictionBatch, PredictionTraceState


class PredictionAggregator:
    @staticmethod
    def empty_result(
        *,
        request: PredictionRequest,
    ) -> KinasePredictionResult:
        empty = pd.DataFrame(index=request.combined_scores.index.copy(), dtype=float)
        return KinasePredictionResult(
            pred_matrix=empty,
            substrate_list={},
            debug_traces={} if request.trace_level != "none" else None,
            trace_level=request.trace_level,
            trace_sink=request.trace_sink,
        )

    @staticmethod
    def initialize_prediction_matrix(
        *,
        feature_mat: pd.DataFrame,
        substrate_list: dict[str, list[str]],
    ) -> pd.DataFrame:
        return pd.DataFrame(
            0.0,
            index=feature_mat.index.copy(),
            columns=list(substrate_list),
        )

    @staticmethod
    def add_kinase_scores(
        *,
        pred_matrix: pd.DataFrame,
        batch: KinasePredictionBatch,
    ) -> None:
        pred_matrix.loc[:, batch.kinase] += batch.scores

    @staticmethod
    def finalize(
        *,
        pred_matrix: pd.DataFrame,
        substrate_list: dict[str, list[str]],
        request: PredictionRequest,
        trace_state: PredictionTraceState,
    ) -> KinasePredictionResult:
        pred_matrix /= float(request.ensemble_size)
        return KinasePredictionResult(
            pred_matrix=pred_matrix,
            substrate_list=substrate_list,
            debug_traces=trace_state.debug_traces,
            trace_level=request.trace_level,
            trace_sink=request.trace_sink,
        )


__all__ = [
    "PredictionAggregator",
]
