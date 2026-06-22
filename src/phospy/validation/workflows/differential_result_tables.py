"""Validation helpers for differential result table reporting utilities."""

from __future__ import annotations

import math
from typing import Literal, cast

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.validation.common.dataframes import (
    require_columns,
    require_dataframe,
    require_unique_columns,
)


def require_differential_result_columns(
    table: object,
    *,
    columns: tuple[str, ...],
    field_name: str,
) -> pd.DataFrame:
    """Require columns needed by a differential result reporting helper."""

    frame = require_dataframe(
        table,
        field_name=field_name,
        allow_empty=True,
        error_type=PhosPyInputError,
    )
    require_unique_columns(
        frame,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_columns(
        frame,
        field_name=field_name,
        required_columns=columns,
        error_type=PhosPyInputError,
    )
    return frame


def require_numeric_result_column(
    table: pd.DataFrame,
    *,
    column_name: str,
    field_name: str,
) -> pd.Series:
    """Return a numeric result column, allowing missing values."""

    require_differential_result_columns(
        table,
        columns=(column_name,),
        field_name=field_name,
    )
    column = table.loc[:, column_name]
    if not isinstance(column, pd.Series):
        raise PhosPyInputError(f"{field_name}.{column_name} must be one column")
    if pd.api.types.is_bool_dtype(column) or _contains_boolean_value(column):
        raise PhosPyInputError(
            f"{field_name}.{column_name} must contain numeric values, not booleans"
        )
    numeric = pd.to_numeric(column, errors="coerce")
    coerced_missing = numeric.isna() & column.notna()
    if bool(coerced_missing.any()):
        raise PhosPyInputError(
            f"{field_name}.{column_name} must contain numeric values when present"
        )
    values = numeric.to_numpy(dtype="float64", copy=False)
    invalid = ~np.isfinite(values) & ~np.isnan(values)
    if bool(invalid.any()):
        raise PhosPyInputError(
            f"{field_name}.{column_name} must contain finite numeric values when present"
        )
    return cast(pd.Series, numeric)


def require_probability_threshold(
    value: object,
    *,
    field_name: str,
) -> float:
    """Require a p-value style threshold in the closed unit interval."""

    threshold = _require_finite_number(value, field_name=field_name)
    if threshold < 0.0 or threshold > 1.0:
        raise PhosPyInputError(f"{field_name} must be between 0.0 and 1.0")
    return threshold


def require_non_negative_threshold(
    value: object,
    *,
    field_name: str,
) -> float:
    """Require a non-negative numeric threshold."""

    threshold = _require_finite_number(value, field_name=field_name)
    if threshold < 0.0:
        raise PhosPyInputError(f"{field_name} must be >= 0.0")
    return threshold


def require_column_name(value: object, *, field_name: str) -> str:
    """Require a non-empty table column name."""

    if not isinstance(value, str) or value.strip() == "":
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return value


def require_boolean(value: object, *, field_name: str) -> bool:
    """Require a boolean option value."""

    if not isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be a bool")
    return value


def require_na_position(value: object, *, field_name: str) -> Literal["first", "last"]:
    """Require a pandas-compatible missing-value sort placement."""

    if value not in {"first", "last"}:
        raise PhosPyInputError(f"{field_name} must be 'first' or 'last'")
    return cast(Literal["first", "last"], value)


def _require_finite_number(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhosPyInputError(f"{field_name} must be numeric")
    number = float(value)
    if not math.isfinite(number):
        raise PhosPyInputError(f"{field_name} must be finite")
    return number


def _contains_boolean_value(column: pd.Series) -> bool:
    return any(isinstance(value, bool) for value in column.dropna().tolist())


__all__ = [
    "require_boolean",
    "require_column_name",
    "require_differential_result_columns",
    "require_na_position",
    "require_non_negative_threshold",
    "require_numeric_result_column",
    "require_probability_threshold",
]
