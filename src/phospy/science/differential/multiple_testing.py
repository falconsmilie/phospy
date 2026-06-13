"""Multiple-testing helpers."""

from __future__ import annotations

import numpy as np

from phospy.science.statistics.multiple_testing import (
    benjamini_hochberg as _benjamini_hochberg,
)


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Return BH-adjusted q-values."""

    return _benjamini_hochberg(p_values)
