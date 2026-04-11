from __future__ import annotations

from typing import TYPE_CHECKING, NoReturn

import pandas as pd

from ..errors import NoCandidateKinasesError
from ..validation.requests import PredictionRequest
from .results import KinasePredictionResult

if TYPE_CHECKING:
    from .execution import KinasePredictionBatch, PredictionTraceState


class PredictionAggregator:
    @staticmethod
    def raise_no_candidate_kinases(
        *,
        request: PredictionRequest,
    ) -> NoReturn:
        msg = (
            "No candidate kinases qualified for prediction from combined_scores "
            f"using top={request.top}, score_threshold={request.score_threshold}, "
            f"and inclusion={request.inclusion}. Lower score_threshold or inclusion, "
            "or increase top."
        )
        raise NoCandidateKinasesError(msg)

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
