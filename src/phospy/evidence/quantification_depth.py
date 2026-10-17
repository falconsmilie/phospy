"""Base quantification-depth evidence contracts and validation."""

from __future__ import annotations

from typing import Final, Literal

import numpy as np
import numpy.typing as npt
import pandas as pd

from phospy.errors.input import PhosPyInputError

QUANTIFICATION_DEPTH_COLUMN: Final[str] = "quantification_depth"
QUANTIFICATION_DEPTH_INTEGER_TOLERANCE: Final[float] = 1.0e-8
QUANTIFICATION_DEPTH_KIND_PSM_COUNT = "psm_count"
QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT = "peptide_count"
QuantificationDepthKind = Literal["psm_count", "peptide_count"]
SUPPORTED_QUANTIFICATION_DEPTH_KINDS: tuple[str, ...] = (
    QUANTIFICATION_DEPTH_KIND_PSM_COUNT,
    QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT,
)


def validate_quantification_depth_series(
    series: pd.Series,
    *,
    field_name: str,
    expected_index: pd.Index | None = None,
) -> pd.Series:
    """Return validated, integer-normalised quantification depths."""

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
    "QUANTIFICATION_DEPTH_KIND_PEPTIDE_COUNT",
    "QUANTIFICATION_DEPTH_KIND_PSM_COUNT",
    "QuantificationDepthKind",
    "SUPPORTED_QUANTIFICATION_DEPTH_KINDS",
    "validate_quantification_depth_series",
]
