"""Neutral diagnostics for caller-declared quantitative intensity scales."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.quantitative.models import IntensityScaleKind


def _run_declared_scale_sanity_diagnostics(
    *,
    declared_scale_kind: IntensityScaleKind,
    phospho: pd.DataFrame,
    total: pd.DataFrame | None,
) -> tuple[str, ...]:
    warnings: list[str] = []
    for matrix_label, matrix in _iter_numeric_matrices(phospho=phospho, total=total):
        finite_values = _finite_numeric_values(matrix)
        if finite_values.size == 0:
            continue
        median_value = float(np.median(finite_values))
        min_value = float(np.min(finite_values))
        max_value = float(np.max(finite_values))
        dynamic_range = max_value - min_value
        integer_like_fraction = _fraction_integer_like(finite_values)
        near_zero_signed_fraction = _fraction_near_zero_signed_values(finite_values)
        likely_log2_fraction = _fraction_within_log2_like_range(finite_values)
        if declared_scale_kind is IntensityScaleKind.LOG2:
            large_value_fraction = float(np.mean(finite_values >= 25.0))
            extreme_value_fraction = float(np.mean(finite_values >= 35.0))
            if (
                min_value >= 0.0
                and median_value >= 20.0
                and large_value_fraction >= 0.8
            ):
                warnings.append(
                    f"{matrix_label}: declared log2 scale is suspicious; "
                    f"median={median_value:.3f}, min={min_value:.3f}, "
                    f"fraction_ge_25={large_value_fraction:.3f}"
                )
            if max_value >= 35.0 or extreme_value_fraction >= 0.2:
                warnings.append(
                    f"{matrix_label}: declared log2 scale has very large positive values; "
                    f"max={max_value:.3f}, fraction_ge_35={extreme_value_fraction:.3f}"
                )
            if min_value <= -40.0 or dynamic_range >= 80.0:
                warnings.append(
                    f"{matrix_label}: declared log2 scale has highly suspicious range; "
                    f"min={min_value:.3f}, max={max_value:.3f}, "
                    f"range={dynamic_range:.3f}"
                )
            if median_value >= 1000.0 and integer_like_fraction >= 0.95:
                warnings.append(
                    f"{matrix_label}: declared log2 scale strongly resembles raw "
                    f"linear intensities; median={median_value:.3f}, "
                    f"integer_like_fraction={integer_like_fraction:.3f}"
                )
            continue
        if declared_scale_kind is IntensityScaleKind.LINEAR:
            negative_fraction = float(np.mean(finite_values < 0.0))
            has_positive = bool(np.any(finite_values > 0.0))
            has_negative = bool(np.any(finite_values < 0.0))
            if negative_fraction > 0.0:
                warnings.append(
                    f"{matrix_label}: declared linear scale contains negative values; "
                    f"fraction_negative={negative_fraction:.3f}"
                )
            if near_zero_signed_fraction >= 0.7 and has_positive and has_negative:
                warnings.append(
                    f"{matrix_label}: declared linear scale shows many small signed "
                    f"values consistent with log-ratio style data; "
                    f"fraction_abs_le_5_with_sign={near_zero_signed_fraction:.3f}"
                )
            if (
                likely_log2_fraction >= 0.95
                and integer_like_fraction <= 0.8
                and max_value <= 20.0
                and median_value <= 8.0
                and (has_negative or near_zero_signed_fraction >= 0.7)
            ):
                warnings.append(
                    f"{matrix_label}: declared linear scale strongly resembles "
                    f"already log-transformed data; "
                    f"fraction_in_log2_range={likely_log2_fraction:.3f}, "
                    f"integer_like_fraction={integer_like_fraction:.3f}"
                )
    return tuple(dict.fromkeys(warnings))


def _iter_numeric_matrices(
    *,
    phospho: pd.DataFrame,
    total: pd.DataFrame | None,
) -> tuple[tuple[str, pd.DataFrame], ...]:
    matrices: list[tuple[str, pd.DataFrame]] = [("phospho", phospho)]
    if total is not None:
        matrices.append(("total", total))
    return tuple(matrices)


def _finite_numeric_values(frame: pd.DataFrame) -> np.ndarray:
    numeric = frame.select_dtypes(include="number")
    if numeric.empty:
        return np.array([], dtype=float)
    values = numeric.to_numpy(dtype=float, copy=False).reshape(-1)
    finite_mask = np.isfinite(values)
    if not bool(finite_mask.any()):
        return np.array([], dtype=float)
    return values[finite_mask]


def _fraction_integer_like(values: np.ndarray) -> float:
    return float(np.mean(np.isclose(values, np.round(values), atol=1e-9)))


def _fraction_near_zero_signed_values(values: np.ndarray) -> float:
    near_zero_mask = np.abs(values) <= 5.0
    if not bool(near_zero_mask.any()):
        return 0.0
    signed_mask = values != 0.0
    return float(np.mean(near_zero_mask & signed_mask))


def _fraction_within_log2_like_range(values: np.ndarray) -> float:
    return float(np.mean((values >= -10.0) & (values <= 20.0)))


__all__: list[str] = []
