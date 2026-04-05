from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from numbers import Real

import pandas as pd

from .errors import PhospyValidationError, TableSchemaError


def validate_fraction(value: float, *, name: str) -> float:
    """Validate a finite fraction-like numeric value in the inclusive 0..1 range."""

    if isinstance(value, bool) or not isinstance(value, Real):
        raise PhospyValidationError(
            f"{name} must be a finite numeric value between 0 and 1"
        )

    resolved = float(value)
    if not math.isfinite(resolved) or not 0.0 <= resolved <= 1.0:
        raise PhospyValidationError(
            f"{name} must be a finite numeric value between 0 and 1"
        )
    return resolved


def resolve_required_columns(
    columns: Iterable[str],
    *,
    argument_name: str,
    context: str,
) -> list[str]:
    """Resolve a required list of column names and reject empty collections."""

    resolved_columns = list(columns)
    if not resolved_columns:
        raise PhospyValidationError(
            f"{context} requires at least one column name in '{argument_name}'"
        )
    return resolved_columns


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
    """Coerce the selected frame columns to numeric values in place."""

    for column in columns:
        frame[column] = require_numeric_series(
            frame[column],
            column=column,
            context=context,
        )
