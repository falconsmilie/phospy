from __future__ import annotations

import numpy as np
import pytest

from phospy.workflows.kinase.scoring_transforms import (
    shift_correlation_to_unit_support,
)


def test_shift_correlation_to_unit_support_respects_nan_bounds_and_input_immutability() -> (
    None
):
    correlation = np.asarray([-1.0, 0.0, 1.0, -1.001, 1.001, np.nan], dtype=float)
    original = correlation.copy()

    scores = shift_correlation_to_unit_support(correlation)

    assert scores[0] == pytest.approx(0.0)
    assert scores[1] == pytest.approx(0.5)
    assert scores[2] == pytest.approx(1.0)
    assert scores[3] == pytest.approx(0.0)
    assert scores[4] == pytest.approx(1.0)
    assert np.isnan(scores[5])
    assert np.array_equal(correlation, original, equal_nan=True)
