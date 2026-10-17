from __future__ import annotations

import numpy as np
import pytest

from phospy.statistics.multiple_testing import (
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


@pytest.mark.parametrize(
    ("method", "expected_finite"),
    [
        (MULTIPLE_TESTING_CORRECTION_NONE, (0.05, 0.01, 0.01, 0.2, 0.05)),
        (
            MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
            (0.0625, 0.025, 0.025, 0.2, 0.0625),
        ),
        (
            MULTIPLE_TESTING_CORRECTION_BONFERRONI,
            (0.25, 0.05, 0.05, 1.0, 0.25),
        ),
        (
            MULTIPLE_TESTING_CORRECTION_HOLM,
            (0.15, 0.05, 0.05, 0.2, 0.15),
        ),
        (
            MULTIPLE_TESTING_CORRECTION_BENJAMINI_YEKUTIELI,
            (
                0.14270833333333333,
                0.05708333333333333,
                0.05708333333333333,
                0.45666666666666667,
                0.14270833333333333,
            ),
        ),
    ],
)
def test_multiple_testing_ties_and_interleaved_missing_values_are_exact(
    method: MultipleTestingCorrection,
    expected_finite: tuple[float, ...],
) -> None:
    adjusted = adjust_p_values(
        (0.05, None, 0.01, np.nan, 0.01, np.inf, 0.2, -np.inf, 0.05),
        method=method,
    )

    assert adjusted.dtype == np.dtype(np.float64)
    assert adjusted.shape == (9,)
    assert np.isnan(adjusted.take((1, 3, 5, 7))).all()
    np.testing.assert_allclose(adjusted.take((0, 2, 4, 6, 8)), expected_finite)


@pytest.mark.parametrize("method", SUPPORTED_MULTIPLE_TESTING_CORRECTIONS)
@pytest.mark.parametrize("value", (0.0, 1.0))
def test_multiple_testing_singleton_boundaries(
    method: MultipleTestingCorrection,
    value: float,
) -> None:
    adjusted = adjust_p_values((value,), method=method)

    assert adjusted.dtype == np.dtype(np.float64)
    assert adjusted.shape == (1,)
    assert adjusted.tolist() == [value]


def test_multiple_testing_rejects_non_vector_and_out_of_range_inputs_exactly() -> None:
    with pytest.raises(ValueError) as two_dimensional:
        adjust_p_values(np.array([[0.1, 0.2]]), method="none")
    assert str(two_dimensional.value) == "p_values must be one-dimensional"

    for values in ((-0.01,), (1.01,)):
        with pytest.raises(ValueError) as out_of_range:
            adjust_p_values(values, method="none")
        assert str(out_of_range.value) == "p_values must be within [0.0, 1.0]"
