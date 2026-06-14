"""Prediction-domain internal frame views for trusted workflow collaborators."""

from __future__ import annotations

import pandas as pd

from phospy.science.prediction.models import KinasePredictionResult, KinaseScoringResult


class KinasePredictionInternalView:
    """Narrow borrowed-frame contract for prediction-result internals."""

    __slots__ = ("_result",)

    def __init__(self, result: KinasePredictionResult) -> None:
        self._result = result

    @property
    def pred_mat(self) -> pd.DataFrame:
        return self._result._borrow_pred_mat_frame()


class KinaseScoringInternalView:
    """Narrow borrowed-frame contract for scoring-result internals."""

    __slots__ = ("_result",)

    def __init__(self, result: KinaseScoringResult) -> None:
        self._result = result

    @property
    def profile_scores(self) -> pd.DataFrame:
        return self._result._borrow_profile_scores_frame()

    @property
    def rank_weighted_fusion_scores(self) -> pd.DataFrame | None:
        return self._result._borrow_rank_weighted_fusion_scores_frame()

    @property
    def authoritative_scores(self) -> pd.DataFrame:
        return self._result._borrow_authoritative_scores_frame()

    @property
    def kinase_library_site_diagnostics(self) -> pd.DataFrame | None:
        return self._result._borrow_kinase_library_site_diagnostics_frame()


__all__ = ["KinasePredictionInternalView", "KinaseScoringInternalView"]
