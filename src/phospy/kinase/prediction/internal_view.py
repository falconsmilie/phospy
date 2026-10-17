"""Kinase prediction-domain internal view for workflow collaborators."""

from __future__ import annotations

import pandas as pd

from phospy.kinase.prediction.results import KinasePredictionResult


class KinasePredictionInternalView:
    """Narrow prediction-domain access to defensive frame snapshots."""

    __slots__ = ("_result",)

    def __init__(self, result: KinasePredictionResult) -> None:
        self._result = result

    @property
    def pred_mat(self) -> pd.DataFrame:
        return self._result._borrow_pred_mat_frame()


__all__ = ["KinasePredictionInternalView"]
