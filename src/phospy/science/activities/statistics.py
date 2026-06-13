"""Statistics helpers for activity-score methods."""

from __future__ import annotations

import math

import pandas as pd

from phospy.science.statistics.multiple_testing import benjamini_hochberg


def two_sided_normal_p_value(z_score: float) -> float:
    """Return two-sided normal-approximation p-value from a z-score."""

    return float(math.erfc(abs(float(z_score)) / math.sqrt(2.0)))


def benjamini_hochberg_q_values(p_values: pd.Series) -> pd.Series:
    """Compute Benjamini-Hochberg adjusted q-values for one p-value series."""

    adjusted = benjamini_hochberg(p_values.to_numpy(dtype=float, copy=False))
    return pd.Series(adjusted, index=p_values.index, dtype=float, name=p_values.name)


__all__ = ["benjamini_hochberg_q_values", "two_sided_normal_p_value"]
