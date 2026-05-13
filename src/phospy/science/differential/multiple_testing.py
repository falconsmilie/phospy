"""Multiple-testing helpers."""

from __future__ import annotations

import numpy as np


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Return BH-adjusted q-values."""

    p_values = np.asarray(p_values, dtype=float)
    adjusted = np.full(p_values.shape, np.nan, dtype=float)
    finite_mask = np.isfinite(p_values)
    if not np.any(finite_mask):
        return adjusted

    finite = p_values[finite_mask]
    order = np.argsort(finite)
    ranked = finite[order]
    n_tests = ranked.size
    ranks = np.arange(1, n_tests + 1, dtype=float)
    raw = ranked * float(n_tests) / ranks
    monotone = np.minimum.accumulate(raw[::-1])[::-1]
    monotone = np.minimum(np.maximum(monotone, 0.0), 1.0)
    restored = np.empty_like(monotone)
    restored[order] = monotone
    adjusted[finite_mask] = restored
    return adjusted
