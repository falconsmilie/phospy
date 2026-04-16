from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, NoReturn

import numpy as np
import pandas as pd

from ..errors import (
    InputCompatibilityError,
    NoCandidateKinasesError,
    format_no_candidate_kinases_message,
)
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
        """Add one kinase batch after strict score-index validation.

        Contract:
        - batch scores are aligned by score_index labels to pred_matrix.index.
        - positional aggregation is never used when score_index is provided.
        - malformed batches raise InputCompatibilityError with seam diagnostics.
        """

        score_values = PredictionAggregator._resolve_batch_scores(
            pred_matrix=pred_matrix,
            batch=batch,
        )
        pred_matrix.values[:, pred_matrix.column_positions[batch.kinase]] += (
            score_values
        )

    @staticmethod
    def _resolve_batch_scores(
        *,
        pred_matrix: PredictionScoreBuffer,
        batch: KinasePredictionBatch,
    ) -> np.ndarray:
        score_values = getattr(batch, "score_values", None)
        score_index = getattr(batch, "score_index", None)
        if score_values is None:
            scores = batch.scores
            if isinstance(scores, pd.Series):
                score_values = scores.to_numpy(copy=False)
                if score_index is None:
                    score_index = scores.index
            else:
                score_values = scores

        if score_index is None:
            msg = (
                "Prediction batch for kinase "
                f"{batch.kinase!r} must provide score_index for index-safe aggregation"
            )
            raise InputCompatibilityError(msg)

        resolved_values = PredictionAggregator._coerce_score_values(
            score_values=score_values,
            kinase=batch.kinase,
        )
        resolved_index = pd.Index(score_index)
        PredictionAggregator._validate_score_index(
            pred_index=pred_matrix.index,
            score_index=resolved_index,
            score_count=len(resolved_values),
            kinase=batch.kinase,
        )
        if resolved_index.equals(pred_matrix.index):
            return resolved_values
        return (
            pd.Series(resolved_values, index=resolved_index, dtype=float)
            .reindex(pred_matrix.index)
            .to_numpy(dtype=float, copy=False)
        )

    @staticmethod
    def _coerce_score_values(
        *,
        score_values: object,
        kinase: str,
    ) -> np.ndarray:
        try:
            resolved = np.asarray(score_values, dtype=float)
        except (TypeError, ValueError) as error:
            msg = (
                f"Prediction batch for kinase {kinase!r} contains non-numeric "
                "score values"
            )
            raise InputCompatibilityError(msg) from error
        if resolved.ndim != 1:
            msg = (
                f"Prediction batch for kinase {kinase!r} must provide exactly one "
                "score value per phosphosite"
            )
            raise InputCompatibilityError(msg)
        return resolved

    @staticmethod
    def _validate_score_index(
        *,
        pred_index: pd.Index,
        score_index: pd.Index,
        score_count: int,
        kinase: str,
    ) -> None:
        if len(score_index) != score_count:
            msg = (
                f"Prediction batch for kinase {kinase!r} has score_index length "
                f"{len(score_index)} but score_values length {score_count}"
            )
            raise InputCompatibilityError(msg)
        if len(score_index) != len(pred_index):
            msg = (
                f"Prediction batch for kinase {kinase!r} has {len(score_index)} "
                f"score labels but prediction matrix has {len(pred_index)} rows"
            )
            raise InputCompatibilityError(msg)
        if not score_index.is_unique:
            duplicate_labels = pd.unique(score_index[score_index.duplicated()])
            msg = (
                f"Prediction batch for kinase {kinase!r} contains duplicate "
                f"score_index labels: {PredictionAggregator._preview_labels(duplicate_labels)}"
            )
            raise InputCompatibilityError(msg)
        if not pred_index.is_unique:
            msg = (
                "prediction matrix index contains duplicate phosphosite labels; "
                "cannot safely aggregate kinase batches by index"
            )
            raise InputCompatibilityError(msg)

        missing_labels = pred_index.difference(score_index)
        unexpected_labels = score_index.difference(pred_index)
        if len(missing_labels) > 0 or len(unexpected_labels) > 0:
            diagnostics: list[str] = []
            if len(missing_labels) > 0:
                diagnostics.append(
                    f"missing labels: {PredictionAggregator._preview_labels(missing_labels)}"
                )
            if len(unexpected_labels) > 0:
                diagnostics.append(
                    "unexpected labels: "
                    f"{PredictionAggregator._preview_labels(unexpected_labels)}"
                )
            msg = (
                f"Prediction batch for kinase {kinase!r} has score_index labels "
                "that do not match prediction matrix index ("
                + "; ".join(diagnostics)
                + ")"
            )
            raise InputCompatibilityError(msg)

    @staticmethod
    def _preview_labels(labels: object, *, max_items: int = 5) -> str:
        label_index = pd.Index(labels)
        if len(label_index) == 0:
            return "<none>"
        preview = ", ".join(repr(label) for label in label_index[:max_items])
        if len(label_index) > max_items:
            preview = f"{preview}, ..."
        return preview

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
