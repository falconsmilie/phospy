"""Numerical transforms used by kinase profile-correlation scoring."""

from __future__ import annotations

import numpy as np


def shift_correlation_to_unit_support(correlation: np.ndarray) -> np.ndarray:
    """Apply the shifted-unit profile support transform."""

    scores = ((correlation + 1.0) / 2.0).copy()
    valid = np.isfinite(scores)
    scores[valid] = np.clip(scores[valid], 0.0, 1.0)
    return scores


__all__ = ["shift_correlation_to_unit_support"]
