"""Prediction and scoring stage result models."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True, slots=True)
class KinaseScoringResult:
    """Scoring-stage outputs.

    `profile_scores` is the required scientific score matrix for the supported
    workflow route. Optional fields remain for compatibility with extended
    scoring lanes and persisted bundles.
    """

    profile_scores: pd.DataFrame
    motif_scores: pd.DataFrame | None = None
    combined_scores: pd.DataFrame | None = None
    weights: pd.DataFrame | None = None

    def __post_init__(self) -> None:
        profile_scores = _copy_required_frame(
            self.profile_scores, field_name="scoring_result.profile_scores"
        )
        motif_scores = _copy_optional_frame(
            self.motif_scores, field_name="scoring_result.motif_scores"
        )
        combined_scores = _copy_optional_frame(
            self.combined_scores, field_name="scoring_result.combined_scores"
        )
        weights = _copy_optional_frame(
            self.weights, field_name="scoring_result.weights"
        )
        object.__setattr__(self, "profile_scores", profile_scores)
        object.__setattr__(self, "motif_scores", motif_scores)
        object.__setattr__(self, "combined_scores", combined_scores)
        object.__setattr__(self, "weights", weights)


@dataclass(frozen=True, slots=True)
class KinasePredictionResult:
    """Prediction-stage outputs."""

    pred_mat: pd.DataFrame
    substrate_list: pd.DataFrame | None = None

    def __post_init__(self) -> None:
        pred_mat = _copy_required_frame(
            self.pred_mat, field_name="prediction_result.pred_mat"
        )
        substrate_list = _copy_optional_frame(
            self.substrate_list, field_name="prediction_result.substrate_list"
        )
        object.__setattr__(self, "pred_mat", pred_mat)
        object.__setattr__(self, "substrate_list", substrate_list)


def _copy_required_frame(value: object, *, field_name: str) -> pd.DataFrame:
    if not isinstance(value, pd.DataFrame):
        raise TypeError(f"{field_name} must be a pandas DataFrame")
    return value.copy(deep=True)


def _copy_optional_frame(
    value: object | None, *, field_name: str
) -> pd.DataFrame | None:
    if value is None:
        return None
    return _copy_required_frame(value, field_name=field_name)
