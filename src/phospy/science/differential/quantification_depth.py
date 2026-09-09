"""Quantification-depth validation for differential empirical-Bayes trends."""

from __future__ import annotations

from typing import Final

import numpy as np
import numpy.typing as npt
import pandas as pd

from phospy.errors.input import PhosPyInputError

QUANTIFICATION_DEPTH_COLUMN: Final[str] = "quantification_depth"
QUANTIFICATION_DEPTH_LOG2_TREND_COVARIATE_NAME: Final[str] = "log2_quantification_depth"
QUANTIFICATION_DEPTH_TREND_TRANSFORMATION: Final[str] = "log2"
QUANTIFICATION_DEPTH_INTEGER_TOLERANCE: Final[float] = 1.0e-8


def validate_quantification_depth_series(
    series: pd.Series,
    *,
    field_name: str,
    expected_index: pd.Index | None = None,
) -> pd.Series:
    """Return validated, integer-normalized quantification depths.

    Depth must be supplied as numeric feature-aligned counts. The returned Series
    preserves the supplied index and uses float values so it can pass through the
    existing numeric differential machinery without extension-dtype ambiguity.
    """

    if not isinstance(series, pd.Series):
        raise PhosPyInputError(f"{field_name} must be a pandas Series")
    if expected_index is not None and not series.index.equals(expected_index):
        raise PhosPyInputError(
            f"{field_name}.index must exactly match differential.matrix.index"
        )

    values_as_object = series.to_numpy(dtype=object, copy=False)
    missing_mask = np.asarray(pd.isna(values_as_object), dtype=bool)
    if bool(missing_mask.any()):
        raise PhosPyInputError(
            f"{field_name} must not contain missing values; "
            f"{_invalid_entry_preview(series=series, invalid_mask=missing_mask)}"
        )

    if pd.api.types.is_bool_dtype(series.dtype):
        raise PhosPyInputError(f"{field_name} must contain numeric count values")
    if not pd.api.types.is_numeric_dtype(series.dtype):
        raise PhosPyInputError(
            f"{field_name} must contain numeric count values; "
            f"observed_dtype={series.dtype}"
        )

    values: npt.NDArray[np.float64] = np.asarray(
        series.to_numpy(dtype=float, copy=False),
        dtype=np.float64,
    )
    finite_mask = np.asarray(np.isfinite(values), dtype=bool)
    if bool((~finite_mask).any()):
        raise PhosPyInputError(
            f"{field_name} must contain finite numeric count values; "
            f"{_invalid_entry_preview(series=series, invalid_mask=~finite_mask)}"
        )

    positive_count_mask = np.asarray(values >= 1.0, dtype=bool)
    if bool((~positive_count_mask).any()):
        raise PhosPyInputError(
            f"{field_name} values must be >= 1; "
            f"{_invalid_entry_preview(series=series, invalid_mask=~positive_count_mask)}"
        )

    integer_mask = np.asarray(
        np.isclose(
            values,
            np.rint(values),
            rtol=0.0,
            atol=QUANTIFICATION_DEPTH_INTEGER_TOLERANCE,
        ),
        dtype=bool,
    )
    if bool((~integer_mask).any()):
        raise PhosPyInputError(
            f"{field_name} count values must be integer-valued within tolerance "
            f"{QUANTIFICATION_DEPTH_INTEGER_TOLERANCE:g}; "
            f"{_invalid_entry_preview(series=series, invalid_mask=~integer_mask)}"
        )

    return pd.Series(
        np.rint(values).astype(float),
        index=series.index.copy(),
        name=QUANTIFICATION_DEPTH_COLUMN,
    )


def log2_quantification_depth_series(
    depth: pd.Series,
    *,
    field_name: str,
    expected_index: pd.Index | None = None,
) -> pd.Series:
    """Return ``log2(quantification_depth)`` after enforcing count validity."""

    validated = validate_quantification_depth_series(
        depth,
        field_name=field_name,
        expected_index=expected_index,
    )
    return pd.Series(
        np.log2(validated.to_numpy(dtype=float, copy=False)),
        index=validated.index.copy(),
        name=QUANTIFICATION_DEPTH_LOG2_TREND_COVARIATE_NAME,
    )


def _invalid_entry_preview(
    *,
    series: pd.Series,
    invalid_mask: npt.NDArray[np.bool_],
) -> str:
    invalid_positions = np.flatnonzero(invalid_mask)
    preview = ", ".join(
        f"({series.index[int(position)]!r}, {series.iloc[int(position)]!r})"
        for position in invalid_positions[:3]
    )
    suffix = (
        ""
        if invalid_positions.size <= 3
        else f", +{int(invalid_positions.size - 3)} more"
    )
    return (
        f"invalid values: {preview}{suffix}; "
        f"invalid_entry_count={int(invalid_positions.size)}"
    )


__all__ = [
    "QUANTIFICATION_DEPTH_COLUMN",
    "QUANTIFICATION_DEPTH_INTEGER_TOLERANCE",
    "QUANTIFICATION_DEPTH_LOG2_TREND_COVARIATE_NAME",
    "QUANTIFICATION_DEPTH_TREND_TRANSFORMATION",
    "log2_quantification_depth_series",
    "validate_quantification_depth_series",
]
