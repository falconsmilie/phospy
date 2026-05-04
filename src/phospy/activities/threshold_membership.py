"""Shared threshold-membership policy for activity substrate selection."""

from __future__ import annotations

import numpy as np
import pandas as pd

THRESHOLD_MEMBERSHIP_RULE = "score >= threshold"


def threshold_membership_mask_array(
    scores: np.ndarray,
    *,
    threshold: float,
) -> np.ndarray:
    """Return membership mask using the canonical activity threshold rule."""

    threshold_value = float(threshold)
    return np.isfinite(scores) & (scores >= threshold_value)


def threshold_membership_mask_frame(
    scores: pd.DataFrame,
    *,
    threshold: float,
) -> pd.DataFrame:
    """Return a DataFrame mask using the canonical activity threshold rule."""

    mask = threshold_membership_mask_array(
        scores.to_numpy(dtype=float, copy=False),
        threshold=threshold,
    )
    return pd.DataFrame(mask, index=scores.index.copy(), columns=scores.columns.copy())


def threshold_membership_filtered_frame(
    scores: pd.DataFrame,
    *,
    threshold: float,
) -> pd.DataFrame:
    """Return score values for members and NaN for non-members."""

    return scores.where(threshold_membership_mask_frame(scores, threshold=threshold))


__all__ = [
    "THRESHOLD_MEMBERSHIP_RULE",
    "threshold_membership_filtered_frame",
    "threshold_membership_mask_array",
    "threshold_membership_mask_frame",
]
