from __future__ import annotations

import numpy as np
import pytest

from phospy.science.statistics.multiple_testing import (
    MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
    MULTIPLE_TESTING_CORRECTION_BENJAMINI_YEKUTIELI,
    MULTIPLE_TESTING_CORRECTION_BONFERRONI,
    MULTIPLE_TESTING_CORRECTION_HOLM,
    MULTIPLE_TESTING_CORRECTION_NONE,
    SUPPORTED_MULTIPLE_TESTING_CORRECTIONS,
    MultipleTestingCorrection,
    adjust_p_values,
    run,
)


@pytest.mark.parametrize(
    ("method", "expected"),
    [
        (
            MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
            (0.025, 0.05, 0.05, 0.01, 0.05),
        ),
        (
            MULTIPLE_TESTING_CORRECTION_BONFERRONI,
            (0.05, 0.2, 0.15, 0.01, 0.25),
        ),
        (
            MULTIPLE_TESTING_CORRECTION_HOLM,
            (0.04, 0.09, 0.09, 0.01, 0.09),
        ),
        (
            MULTIPLE_TESTING_CORRECTION_BENJAMINI_YEKUTIELI,
            (
                0.05708333333333333,
                0.11416666666666667,
                0.11416666666666667,
                0.022833333333333334,
                0.11416666666666667,
            ),
        ),
        (
            MULTIPLE_TESTING_CORRECTION_NONE,
            (0.01, 0.04, 0.03, 0.002, 0.05),
        ),
    ],
)
def test_multiple_testing_correction_known_vectors(
    method: MultipleTestingCorrection,
    expected: tuple[float, ...],
) -> None:
    adjusted = run(
        (0.01, 0.04, 0.03, 0.002, 0.05),
        method=method,
    )

    assert adjusted == pytest.approx(expected)


@pytest.mark.parametrize(
    ("method", "expected"),
    [
        (
            MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
            (0.03, np.nan, 0.04, 0.04),
        ),
        (
            MULTIPLE_TESTING_CORRECTION_BONFERRONI,
            (0.03, np.nan, 0.12, 0.09),
        ),
        (
            MULTIPLE_TESTING_CORRECTION_HOLM,
            (0.03, np.nan, 0.06, 0.06),
        ),
        (
            MULTIPLE_TESTING_CORRECTION_BENJAMINI_YEKUTIELI,
            (0.055, np.nan, 0.07333333333333333, 0.07333333333333333),
        ),
        (
            MULTIPLE_TESTING_CORRECTION_NONE,
            (0.01, np.nan, 0.04, 0.03),
        ),
    ],
)
def test_multiple_testing_correction_preserves_nan_positions(
    method: MultipleTestingCorrection,
    expected: tuple[float, ...],
) -> None:
    adjusted = adjust_p_values(
        np.array([0.01, np.nan, 0.04, 0.03]),
        method=method,
    )

    np.testing.assert_allclose(
        adjusted,
        np.array(expected),
        rtol=1e-12,
        atol=0.0,
        equal_nan=True,
    )


@pytest.mark.parametrize("method", SUPPORTED_MULTIPLE_TESTING_CORRECTIONS)
def test_multiple_testing_correction_preserves_nonfinite_positions_as_missing(
    method: MultipleTestingCorrection,
) -> None:
    adjusted = adjust_p_values(
        np.array([0.01, np.inf, 0.04, -np.inf, 0.03]),
        method=method,
    )

    assert adjusted.shape == (5,)
    assert np.isnan(adjusted[1])
    assert np.isnan(adjusted[3])
    assert np.isfinite(adjusted[np.array([True, False, True, False, True])]).all()

    returned = run((0.01, np.inf, 0.04, -np.inf, 0.03), method=method)
    assert returned[1] is None
    assert returned[3] is None


@pytest.mark.parametrize("method", SUPPORTED_MULTIPLE_TESTING_CORRECTIONS)
def test_multiple_testing_correction_all_nan_vector(
    method: MultipleTestingCorrection,
) -> None:
    adjusted = adjust_p_values((np.nan, None, np.inf), method=method)

    assert adjusted.shape == (3,)
    assert np.isnan(adjusted).all()
    assert run((np.nan, None, np.inf), method=method) == (None, None, None)


@pytest.mark.parametrize("method", SUPPORTED_MULTIPLE_TESTING_CORRECTIONS)
def test_multiple_testing_correction_empty_vector(
    method: MultipleTestingCorrection,
) -> None:
    adjusted = adjust_p_values((), method=method)

    assert adjusted.shape == (0,)
    assert run((), method=method) == ()


def test_multiple_testing_correction_preserves_input_order() -> None:
    adjusted = run(
        (0.04, 0.001, 0.03),
        method=MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
    )

    assert adjusted == pytest.approx((0.04, 0.003, 0.04))


def test_multiple_testing_correction_rejects_storey_q_values() -> None:
    with pytest.raises(
        ValueError,
        match="multiple-testing correction method must be one of",
    ):
        adjust_p_values((0.1, 0.2), method="storey")
