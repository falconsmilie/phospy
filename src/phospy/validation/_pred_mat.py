from __future__ import annotations

import pandas as pd

from ..prediction.models import PredMatResult


def normalize_pred_mat_input(
    pred_mat: pd.DataFrame | PredMatResult | None,
) -> pd.DataFrame | None:
    """Normalise public predMat inputs to the internal DataFrame contract.

    This helper unwraps the canonical ``PredMatResult`` boundary object into its
    owned in-memory DataFrame without taking an additional defensive copy.
    Schema validation remains responsible for ownership transfer when raw public
    inputs cross into a trusted validated request.
    """

    if isinstance(pred_mat, PredMatResult):
        return pred_mat.to_frame(copy=False)
    return pred_mat


__all__ = ["normalize_pred_mat_input"]
