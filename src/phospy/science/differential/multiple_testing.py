"""Multiple-testing helpers."""

from __future__ import annotations

import numpy as np

from phospy.science.statistics.multiple_testing import (
    MultipleTestingCorrection,
)
from phospy.science.statistics.multiple_testing import (
    adjust_p_values as _adjust_p_values,
)
from phospy.science.statistics.multiple_testing import (
    benjamini_hochberg as _benjamini_hochberg,
)


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Return BH-adjusted q-values."""

    return _benjamini_hochberg(p_values)


def adjust_p_values(
    p_values: np.ndarray,
    *,
    method: MultipleTestingCorrection,
) -> np.ndarray:
    """Return adjusted p-values using the shared statistics helper."""

    return _adjust_p_values(p_values, method=method)


__all__ = ["adjust_p_values", "benjamini_hochberg"]
