from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np
import pandas as pd

from ...errors import TableSchemaError


def require_dataframe(frame: pd.DataFrame, *, context: str) -> pd.DataFrame:
    """Validate the raw boundary object type without taking ownership yet."""

    if not isinstance(frame, pd.DataFrame):
        msg = f"{context} must be a pandas DataFrame"
        raise TableSchemaError(msg)
    return frame


def require_columns(
    frame: pd.DataFrame,
    *,
    required_columns: Sequence[str],
    context: str,
) -> None:
    """Ensure a frame exposes the required column set."""

    missing_columns = [
        column for column in required_columns if column not in frame.columns
    ]
    if missing_columns:
        joined_columns = ", ".join(missing_columns)
        raise TableSchemaError(
            f"{context} is missing required columns: {joined_columns}"
        )


def require_unique_columns(columns: Iterable[object], *, context: str) -> None:
    """Reject duplicate column names."""

    seen: set[str] = set()
    duplicates: list[str] = []
    for column in columns:
        value = str(column)
        if value in seen and value not in duplicates:
            duplicates.append(value)
        seen.add(value)
    if duplicates:
        duplicates_str = ", ".join(duplicates)
        msg = f"{context} contains duplicate column names: {duplicates_str}"
        raise TableSchemaError(msg)


def require_non_null_column_names(columns: Iterable[object], *, context: str) -> None:
    """Reject null-like column labels."""

    null_like: list[str] = []
    for column in columns:
        if pd.isna(column):
            null_like.append("<null>")
    if null_like:
        msg = f"{context} contains null column names"
        raise TableSchemaError(msg)


def require_non_null_values(
    frame: pd.DataFrame,
    *,
    columns: Sequence[str],
    context: str,
) -> None:
    """Reject null values in required columns."""

    failures: list[str] = []
    for column in columns:
        if frame[column].isna().any():
            failures.append(column)
    if failures:
        failures_str = ", ".join(failures)
        msg = f"{context} contains null values in required columns: {failures_str}"
        raise TableSchemaError(msg)


def require_numeric_series(
    values: pd.Series,
    *,
    column: str,
    context: str,
) -> pd.Series:
    """Coerce a Series to numeric values or raise a table schema error."""

    try:
        numeric_values = pd.to_numeric(values, errors="raise")
    except (TypeError, ValueError) as exc:
        raise TableSchemaError(
            f"{context} requires numeric values in column '{column}'"
        ) from exc

    return pd.Series(numeric_values, index=values.index, copy=False)


def require_numeric_columns(
    frame: pd.DataFrame,
    *,
    columns: Sequence[str],
    context: str,
) -> None:
    """Coerce selected frame columns to numeric values in place."""

    for column in columns:
        frame[column] = require_numeric_series(
            frame[column],
            column=column,
            context=context,
        )


def coerce_numeric_columns(
    frame: pd.DataFrame,
    *,
    columns: Sequence[str],
    context: str,
    copy_frame: bool = True,
) -> pd.DataFrame:
    """Return a numeric frame, copying only when taking external ownership.

    External validation boundaries should keep the default ``copy_frame=True`` so
    callers retain isolation from later mutation. Internal trusted flows may pass
    ``copy_frame=False`` when they already own the frame and want in-place numeric
    coercion instead of another full-frame copy.
    """

    validated = frame.copy(deep=True) if copy_frame else frame
    failures: list[str] = []
    for column in columns:
        converted = pd.to_numeric(validated[column], errors="coerce")
        invalid_mask = validated[column].notna() & converted.isna()
        if invalid_mask.any():
            sample_values = validated.loc[invalid_mask, column].astype(str).unique()[:3]
            sample_preview = ", ".join(str(value) for value in sample_values)
            failures.append(f"{column} ({sample_preview})")
        validated[column] = converted
    if failures:
        failures_str = "; ".join(failures)
        msg = (
            f"{context} contains non-numeric values in numeric columns: {failures_str}"
        )
        raise TableSchemaError(msg)
    return validated


def require_value_range(
    frame: pd.DataFrame,
    *,
    columns: Sequence[str],
    minimum: float,
    maximum: float,
    context: str,
) -> None:
    """Validate inclusive numeric bounds for selected columns."""

    failures: list[str] = []
    for column in columns:
        series = frame[column]
        out_of_range = series.notna() & ((series < minimum) | (series > maximum))
        if out_of_range.any():
            failures.append(column)
    if failures:
        failures_str = ", ".join(failures)
        msg = (
            f"{context} contains values outside the allowed range "
            f"[{minimum}, {maximum}] in columns: {failures_str}"
        )
        raise TableSchemaError(msg)


def require_no_infinite_numeric_values(
    frame: pd.DataFrame,
    *,
    columns: Sequence[str],
    context: str,
) -> None:
    """Reject positive or negative infinity in selected numeric columns."""

    failures: list[str] = []
    for column in columns:
        series = frame[column]
        values = series.to_numpy(dtype=float)
        invalid_mask = np.isinf(values)
        if invalid_mask.any():
            sample_values = series.loc[invalid_mask].astype(str).unique()[:3]
            sample_preview = ", ".join(str(value) for value in sample_values)
            failures.append(f"{column} ({sample_preview})")
    if failures:
        failures_str = "; ".join(failures)
        msg = f"{context} contains infinite values in numeric columns: {failures_str}"
        raise TableSchemaError(msg)


def require_finite_numeric_values(
    frame: pd.DataFrame,
    *,
    columns: Sequence[str],
    context: str,
) -> None:
    """Reject NaN or infinite numeric values in selected columns."""

    failures: list[str] = []
    for column in columns:
        series = frame[column]
        invalid_mask = ~np.isfinite(series.to_numpy(dtype=float))
        if invalid_mask.any():
            sample_values = series.loc[invalid_mask].astype(str).unique()[:3]
            sample_preview = ", ".join(str(value) for value in sample_values)
            failures.append(f"{column} ({sample_preview})")
    if failures:
        failures_str = "; ".join(failures)
        msg = f"{context} contains non-finite values in numeric columns: {failures_str}"
        raise TableSchemaError(msg)


def require_unique_index(frame: pd.DataFrame, *, context: str) -> None:
    """Reject duplicate DataFrame index entries."""

    duplicated = frame.index.duplicated(keep=False)
    if duplicated.any():
        duplicates = frame.index[duplicated].astype(str).unique().tolist()[:5]
        duplicates_str = ", ".join(duplicates)
        msg = f"{context} contains duplicate index entries: {duplicates_str}"
        raise TableSchemaError(msg)


def require_non_null_index(frame: pd.DataFrame, *, context: str) -> None:
    """Reject null index entries."""

    if frame.index.isna().any():
        msg = f"{context} contains null index entries"
        raise TableSchemaError(msg)


__all__ = [
    "coerce_numeric_columns",
    "require_columns",
    "require_dataframe",
    "require_finite_numeric_values",
    "require_no_infinite_numeric_values",
    "require_non_null_column_names",
    "require_non_null_index",
    "require_non_null_values",
    "require_numeric_columns",
    "require_numeric_series",
    "require_unique_columns",
    "require_unique_index",
    "require_value_range",
]
