"""Prediction and scoring stage result models."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True, slots=True)
class KinaseScoringResult:
    """Scoring-stage outputs."""

    profile_scores: pd.DataFrame
    motif_scores: pd.DataFrame | None = None
    combined_scores: pd.DataFrame | None = None
    weights: pd.DataFrame | None = None


@dataclass(frozen=True, slots=True)
class KinasePredictionResult:
    """Prediction-stage outputs."""

    pred_mat: pd.DataFrame
    substrate_list: pd.DataFrame | None = None
