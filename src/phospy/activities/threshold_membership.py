"""Shared threshold-membership policy for activity substrate selection."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.policy_models import ThresholdMode

THRESHOLD_MEMBERSHIP_MODE = ThresholdMode.GREATER_THAN_OR_EQUAL
THRESHOLD_MEMBERSHIP_RULE = THRESHOLD_MEMBERSHIP_MODE.value


def threshold_membership_mask_array(
    scores: np.ndarray,
    *,
    threshold: float,
    threshold_mode: ThresholdMode | str = THRESHOLD_MEMBERSHIP_MODE,
) -> np.ndarray:
    """Return membership mask using the canonical activity threshold rule."""

    threshold_value = float(threshold)
    return _threshold_comparison_mask(
        scores=scores,
        threshold_value=threshold_value,
        threshold_mode=threshold_mode,
    )


def threshold_membership_mask_frame(
    scores: pd.DataFrame,
    *,
    threshold: float,
    threshold_mode: ThresholdMode | str = THRESHOLD_MEMBERSHIP_MODE,
) -> pd.DataFrame:
    """Return a DataFrame mask using the canonical activity threshold rule."""

    mask = threshold_membership_mask_array(
        scores.to_numpy(dtype=float, copy=False),
        threshold=threshold,
        threshold_mode=threshold_mode,
    )
    return pd.DataFrame(mask, index=scores.index.copy(), columns=scores.columns.copy())


def threshold_membership_filtered_frame(
    scores: pd.DataFrame,
    *,
    threshold: float,
    threshold_mode: ThresholdMode | str = THRESHOLD_MEMBERSHIP_MODE,
) -> pd.DataFrame:
    """Return score values for members and NaN for non-members."""

    return scores.where(
        threshold_membership_mask_frame(
            scores,
            threshold=threshold,
            threshold_mode=threshold_mode,
        )
    )


def _threshold_comparison_mask(
    *,
    scores: np.ndarray,
    threshold_value: float,
    threshold_mode: ThresholdMode | str,
) -> np.ndarray:
    mode = ThresholdMode.parse(
        threshold_mode,
        field_name="activity threshold membership mode",
    )
    finite = np.isfinite(scores)
    if mode is ThresholdMode.GREATER_THAN_OR_EQUAL:
        return finite & (scores >= threshold_value)
    if mode is ThresholdMode.GREATER_THAN:
        return finite & (scores > threshold_value)
    return finite & (scores >= threshold_value)


__all__ = [
    "THRESHOLD_MEMBERSHIP_MODE",
    "THRESHOLD_MEMBERSHIP_RULE",
    "threshold_membership_filtered_frame",
    "threshold_membership_mask_array",
    "threshold_membership_mask_frame",
]
