from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, NoReturn

import numpy as np
import pandas as pd

from ..errors import NoCandidateKinasesError, format_no_candidate_kinases_message
from ..validation.requests.prediction import PredictionRequest
from .candidates import summarize_candidate_shortfall
from .results import KinasePredictionResult

if TYPE_CHECKING:
    from .execution import KinasePredictionBatch, PredictionTraceState


@dataclass(slots=True)
class PredictionScoreBuffer:
    values: np.ndarray
    index: pd.Index
    columns: pd.Index
    column_positions: dict[str, int]


class PredictionAggregator:
    @staticmethod
    def raise_no_candidate_kinases(
        *,
        request: PredictionRequest,
    ) -> NoReturn:
        diagnostics = summarize_candidate_shortfall(
            request.combined_scores,
            top=request.top,
            score_threshold=request.score_threshold,
            inclusion=request.inclusion,
        )
        msg = format_no_candidate_kinases_message(
            source_name="combined_scores",
            top=request.top,
            score_threshold=request.score_threshold,
            inclusion=request.inclusion,
            kinase_count=diagnostics.kinase_count,
            site_count=diagnostics.site_count,
            effective_top=diagnostics.effective_top,
            qualifying_kinases=diagnostics.qualifying_kinases,
            max_qualifying_sites=diagnostics.max_qualifying_sites,
            near_miss_kinases=diagnostics.near_miss_kinases,
        )
        raise NoCandidateKinasesError(msg)

    @staticmethod
    def initialize_prediction_matrix(
        *,
        feature_mat: pd.DataFrame,
        substrate_list: dict[str, list[str]],
    ) -> PredictionScoreBuffer:
        columns = pd.Index(list(substrate_list))
        values = np.zeros((len(feature_mat.index), len(columns)), dtype=float)
        return PredictionScoreBuffer(
            values=values,
            index=feature_mat.index.copy(),
            columns=columns,
            column_positions={kinase: i for i, kinase in enumerate(columns)},
        )

    @staticmethod
    def add_kinase_scores(
        *,
        pred_matrix: PredictionScoreBuffer,
        batch: KinasePredictionBatch,
    ) -> None:
        score_values = getattr(batch, "score_values", None)
        if score_values is None:
            scores = batch.scores
            if isinstance(scores, pd.Series):
                score_values = scores.to_numpy(dtype=float, copy=False)
            else:
                score_values = np.asarray(scores, dtype=float)
        pred_matrix.values[:, pred_matrix.column_positions[batch.kinase]] += np.asarray(
            score_values,
            dtype=float,
        )

    @staticmethod
    def finalize(
        *,
        pred_matrix: PredictionScoreBuffer,
        substrate_list: dict[str, list[str]],
        request: PredictionRequest,
        trace_state: PredictionTraceState,
    ) -> KinasePredictionResult:
        pred_matrix.values /= float(request.ensemble_size)
        pred_matrix_frame = pd.DataFrame(
            pred_matrix.values,
            index=pred_matrix.index,
            columns=pred_matrix.columns,
            copy=False,
        )
        return KinasePredictionResult(
            pred_matrix=pred_matrix_frame,
            substrate_list=substrate_list,
            debug_traces=trace_state.debug_traces,
            trace_level=request.trace_level,
            trace_sink=request.trace_sink,
        )


__all__ = [
    "PredictionScoreBuffer",
    "PredictionAggregator",
]
