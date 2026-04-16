from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ...errors import CustomPredictorOutputError, InputCompatibilityError
from ...prediction.contracts import (
    CUSTOM_PREDICTOR_SCORE_INDEX_CONTRACT,
    EnsemblePredictorContract,
)

if TYPE_CHECKING:
    from ...prediction.execution import KinasePredictionBatch


def validate_ensemble_predictor(
    ensemble_predictor: EnsemblePredictorContract | None,
) -> EnsemblePredictorContract | None:
    """Validate an injected ensemble predictor against the prediction contract."""

    if ensemble_predictor is None:
        return None
    if not isinstance(ensemble_predictor, EnsemblePredictorContract):
        msg = (
            "ensemble_predictor must implement EnsemblePredictorContract with "
            "predict_kinase(..., sampling_session=...)"
        )
        raise InputCompatibilityError(msg)
    return ensemble_predictor


def validate_kinase_prediction_batch(
    *,
    batch: object,
    requested_kinase: str,
    feature_index: pd.Index,
    score_index_contract: str = CUSTOM_PREDICTOR_SCORE_INDEX_CONTRACT,
) -> KinasePredictionBatch:
    """Validate and normalize a custom predictor output batch at the boundary."""
    from ...prediction.execution import KinasePredictionBatch

    if score_index_contract != CUSTOM_PREDICTOR_SCORE_INDEX_CONTRACT:
        msg = (
            "Unsupported score-index contract: "
            f"{score_index_contract!r}. Expected "
            f"{CUSTOM_PREDICTOR_SCORE_INDEX_CONTRACT!r}."
        )
        raise InputCompatibilityError(msg)

    if not isinstance(batch, KinasePredictionBatch):
        raise _predictor_output_error(
            requested_kinase=requested_kinase,
            reason="predict_kinase(...) must return KinasePredictionBatch",
        )

    if batch.kinase != requested_kinase:
        raise _predictor_output_error(
            requested_kinase=requested_kinase,
            reason=(
                f"returned kinase {batch.kinase!r} but request was for "
                f"{requested_kinase!r}"
            ),
        )

    resolved_values = _coerce_score_values(
        score_values=batch.score_values,
        requested_kinase=requested_kinase,
        expected_score_count=len(feature_index),
    )
    resolved_index = pd.Index(batch.score_index)
    _validate_score_index(
        feature_index=feature_index,
        score_index=resolved_index,
        score_count=len(resolved_values),
        requested_kinase=requested_kinase,
    )

    non_finite_positions = np.flatnonzero(~np.isfinite(resolved_values))
    if len(non_finite_positions) > 0:
        raise _predictor_output_error(
            requested_kinase=requested_kinase,
            reason=(
                "score_values contain non-finite values at "
                f"{_preview_non_finite(resolved_index, resolved_values, non_finite_positions)}"
            ),
        )

    if resolved_index.equals(feature_index):
        return KinasePredictionBatch(
            kinase=requested_kinase,
            score_values=resolved_values,
            score_index=resolved_index,
        )

    aligned_values = (
        pd.Series(resolved_values, index=resolved_index, dtype=float)
        .reindex(feature_index)
        .to_numpy(dtype=float, copy=False)
    )
    return KinasePredictionBatch(
        kinase=requested_kinase,
        score_values=aligned_values,
        score_index=feature_index,
    )


def _coerce_score_values(
    *,
    score_values: object,
    requested_kinase: str,
    expected_score_count: int,
) -> np.ndarray:
    try:
        resolved = np.asarray(score_values, dtype=float)
    except (TypeError, ValueError) as error:
        raise _predictor_output_error(
            requested_kinase=requested_kinase,
            reason="score_values contain non-numeric values",
        ) from error
    if resolved.ndim != 1:
        raise _predictor_output_error(
            requested_kinase=requested_kinase,
            reason="score_values must be a 1D numeric array",
        )
    if len(resolved) != expected_score_count:
        raise _predictor_output_error(
            requested_kinase=requested_kinase,
            reason=(
                f"score_values length {len(resolved)} does not match expected "
                f"feature row count {expected_score_count}"
            ),
        )
    return resolved


def _validate_score_index(
    *,
    feature_index: pd.Index,
    score_index: pd.Index,
    score_count: int,
    requested_kinase: str,
) -> None:
    if len(score_index) != score_count:
        raise _predictor_output_error(
            requested_kinase=requested_kinase,
            reason=(
                f"score_index length {len(score_index)} does not match "
                f"score_values length {score_count}"
            ),
        )
    if len(score_index) != len(feature_index):
        raise _predictor_output_error(
            requested_kinase=requested_kinase,
            reason=(
                f"score_index length {len(score_index)} does not match expected "
                f"feature row count {len(feature_index)}"
            ),
        )
    if not score_index.is_unique:
        duplicate_labels = pd.unique(score_index[score_index.duplicated()])
        raise _predictor_output_error(
            requested_kinase=requested_kinase,
            reason=(
                "score_index contains duplicate labels: "
                f"{_preview_labels(duplicate_labels)}"
            ),
        )
    if not feature_index.is_unique:
        msg = (
            "Prediction feature index contains duplicate phosphosite labels; "
            "cannot safely align custom predictor batch scores by index"
        )
        raise InputCompatibilityError(msg)

    missing_labels = feature_index.difference(score_index)
    unexpected_labels = score_index.difference(feature_index)
    if len(missing_labels) > 0 or len(unexpected_labels) > 0:
        diagnostics: list[str] = []
        if len(missing_labels) > 0:
            diagnostics.append(f"missing labels: {_preview_labels(missing_labels)}")
        if len(unexpected_labels) > 0:
            diagnostics.append(
                f"unexpected labels: {_preview_labels(unexpected_labels)}"
            )
        raise _predictor_output_error(
            requested_kinase=requested_kinase,
            reason=(
                "score_index labels do not match the requested feature index ("
                + "; ".join(diagnostics)
                + ")"
            ),
        )


def _preview_labels(labels: object, *, max_items: int = 5) -> str:
    label_index = pd.Index(labels)
    if len(label_index) == 0:
        return "<none>"
    preview = ", ".join(repr(label) for label in label_index[:max_items])
    if len(label_index) > max_items:
        preview = f"{preview}, ..."
    return preview


def _preview_non_finite(
    score_index: pd.Index,
    score_values: np.ndarray,
    positions: np.ndarray,
    *,
    max_items: int = 5,
) -> str:
    items: list[str] = []
    for position in positions[:max_items]:
        label = score_index[int(position)]
        value = score_values[int(position)]
        items.append(f"{label!r}={value!r}")
    if len(positions) > max_items:
        items.append("...")
    return ", ".join(items)


def _predictor_output_error(
    *,
    requested_kinase: str,
    reason: str,
) -> CustomPredictorOutputError:
    return CustomPredictorOutputError(
        f"Custom predictor output invalid for kinase {requested_kinase!r}: {reason}"
    )


__all__ = ["validate_ensemble_predictor", "validate_kinase_prediction_batch"]
