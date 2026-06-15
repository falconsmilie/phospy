from __future__ import annotations

import numpy as np

from phospy.science.statistics.multiple_testing import benjamini_hochberg


def test_bh_adjustment_uses_finite_p_value_count_as_denominator() -> None:
    adjusted = benjamini_hochberg(np.array([0.01, np.nan, 0.04, 0.03]))

    np.testing.assert_allclose(
        adjusted[np.array([True, False, True, True])],
        np.array([0.03, 0.04, 0.04]),
        rtol=1e-12,
        atol=0.0,
    )
    assert np.isnan(adjusted[1])


def test_bh_adjustment_preserves_nan_positions() -> None:
    adjusted = benjamini_hochberg(np.array([np.nan, 0.2, np.nan, 0.05]))

    assert np.isnan(adjusted[0])
    assert np.isnan(adjusted[2])
    np.testing.assert_allclose(
        adjusted[np.array([False, True, False, True])],
        np.array([0.2, 0.1]),
        rtol=1e-12,
        atol=0.0,
    )


def test_bh_adjustment_preserves_inf_positions_as_nan() -> None:
    adjusted = benjamini_hochberg(np.array([0.01, np.inf, 0.04, -np.inf, 0.03]))

    assert adjusted.shape == (5,)
    assert np.isnan(adjusted[1])
    assert np.isnan(adjusted[3])
    np.testing.assert_allclose(
        adjusted[np.array([True, False, True, False, True])],
        np.array([0.03, 0.04, 0.04]),
        rtol=1e-12,
        atol=0.0,
    )
