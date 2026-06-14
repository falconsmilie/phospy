"""Multiple-testing correction helpers."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Literal, cast

import numpy as np

MULTIPLE_TESTING_CORRECTION_NONE = "none"
MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG = "benjamini_hochberg"
MultipleTestingCorrection = Literal["none", "benjamini_hochberg"]
SUPPORTED_MULTIPLE_TESTING_CORRECTIONS: tuple[MultipleTestingCorrection, ...] = (
    MULTIPLE_TESTING_CORRECTION_NONE,
    MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
)


def run(
    p_values: Sequence[float | None],
    *,
    method: MultipleTestingCorrection,
) -> tuple[float | None, ...]:
    """Return adjusted p-values for ``p_values`` using ``method``.

    Non-finite adjusted values are returned as ``None``.
    """

    adjusted = adjust_p_values(p_values, method=method)
    return tuple(
        None if not np.isfinite(value) else float(value) for value in adjusted.tolist()
    )


def adjust_p_values(
    p_values: Sequence[float | None] | np.ndarray,
    *,
    method: MultipleTestingCorrection | str,
) -> np.ndarray:
    """Return adjusted p-values as a one-dimensional float array."""

    resolved_method = _require_correction_method(method)
    values = _normalise_p_values(p_values)
    _validate_unit_interval(values)
    if resolved_method == MULTIPLE_TESTING_CORRECTION_NONE:
        return values.copy()
    if resolved_method == MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG:
        return benjamini_hochberg(values)
    raise ValueError(f"unsupported multiple-testing correction: {resolved_method!r}")


def benjamini_hochberg(p_values: Sequence[float | None] | np.ndarray) -> np.ndarray:
    """Return Benjamini-Hochberg adjusted p-values.

    Only finite input p-values are ranked and adjusted. The BH denominator is
    the number of finite p-values passed to this helper, and non-finite input
    positions are preserved as ``NaN`` in the adjusted output.
    """

    values = _normalise_p_values(p_values)
    _validate_unit_interval(values)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    finite_mask = np.isfinite(values)
    if not np.any(finite_mask):
        return adjusted

    finite = values[finite_mask]
    order = np.argsort(finite, kind="mergesort")
    ranked = finite[order]
    n_tests = int(ranked.size)
    ranks = np.arange(1, n_tests + 1, dtype=float)
    raw_adjusted = ranked * float(n_tests) / ranks
    monotone = np.minimum.accumulate(raw_adjusted[::-1])[::-1]
    monotone = np.clip(monotone, 0.0, 1.0)

    restored = np.empty_like(monotone)
    restored[order] = monotone
    adjusted[finite_mask] = restored
    return adjusted


def _normalise_p_values(
    p_values: Sequence[float | None] | np.ndarray,
) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    if values.ndim != 1:
        raise ValueError("p_values must be one-dimensional")
    return values


def _validate_unit_interval(values: np.ndarray) -> None:
    finite_mask = np.isfinite(values)
    invalid_mask = finite_mask & ((values < 0.0) | (values > 1.0))
    if np.any(invalid_mask):
        raise ValueError("p_values must be within [0.0, 1.0]")


def _require_correction_method(value: object) -> MultipleTestingCorrection:
    if (
        not isinstance(value, str)
        or value not in SUPPORTED_MULTIPLE_TESTING_CORRECTIONS
    ):
        allowed = ", ".join(
            repr(method) for method in SUPPORTED_MULTIPLE_TESTING_CORRECTIONS
        )
        raise ValueError(
            f"multiple-testing correction method must be one of: {allowed}"
        )
    return cast(MultipleTestingCorrection, value)


__all__ = [
    "MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG",
    "MULTIPLE_TESTING_CORRECTION_NONE",
    "SUPPORTED_MULTIPLE_TESTING_CORRECTIONS",
    "MultipleTestingCorrection",
    "adjust_p_values",
    "benjamini_hochberg",
    "run",
]
