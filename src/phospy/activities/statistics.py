"""Statistics helpers for activity-score methods."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd


def two_sided_normal_p_value(z_score: float) -> float:
    """Return two-sided normal-approximation p-value from a z-score."""

    return float(math.erfc(abs(float(z_score)) / math.sqrt(2.0)))


def benjamini_hochberg_q_values(p_values: pd.Series) -> pd.Series:
    """Compute Benjamini-Hochberg adjusted q-values for one p-value series."""

    if p_values.empty:
        return pd.Series(dtype=float, index=p_values.index, name=p_values.name)
    finite_mask = p_values.notna() & np.isfinite(p_values.to_numpy(dtype=float))
    finite_p = p_values.loc[finite_mask].astype(float)
    if finite_p.empty:
        return pd.Series(np.nan, index=p_values.index, dtype=float, name=p_values.name)

    sorted_index = finite_p.sort_values(kind="mergesort").index
    sorted_p = finite_p.loc[sorted_index].to_numpy(dtype=float, copy=False)
    m = int(sorted_p.size)
    ranks = np.arange(1, m + 1, dtype=float)
    raw = sorted_p * m / ranks
    adjusted = np.minimum.accumulate(raw[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)

    q_values = pd.Series(np.nan, index=p_values.index, dtype=float, name=p_values.name)
    q_values.loc[sorted_index] = adjusted
    return q_values


__all__ = ["benjamini_hochberg_q_values", "two_sided_normal_p_value"]
