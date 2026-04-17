"""Shared DataFrame-level validation helpers."""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

from phospy.errors.validation import PhosPyValidationError

ValidationErrorType = type[PhosPyValidationError]


def require_dataframe(
    value: object,
    *,
    field_name: str,
    allow_empty: bool,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require a pandas DataFrame and optionally reject empty frames."""

    if not isinstance(value, pd.DataFrame):
        raise error_type(f"{field_name} must be a pandas DataFrame")
    if not allow_empty and value.empty:
        raise error_type(f"{field_name} must be non-empty")
    return value


def require_numeric_dataframe(
    value: pd.DataFrame,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require all columns in a DataFrame to be numeric."""

    non_numeric_columns = [
        str(column)
        for column in value.columns
        if not pd.api.types.is_numeric_dtype(value[column])
    ]
    if non_numeric_columns:
        joined_columns = ", ".join(non_numeric_columns)
        raise error_type(
            f"{field_name} must contain only numeric columns; non-numeric columns: "
            f"{joined_columns}"
        )
    return value


def require_unique_index(
    value: pd.DataFrame,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require unique index labels in a DataFrame."""

    if not value.index.is_unique:
        raise error_type(f"{field_name}.index must be unique")
    return value


def require_unique_columns(
    value: pd.DataFrame,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require unique column labels in a DataFrame."""

    if not value.columns.is_unique:
        raise error_type(f"{field_name}.columns must be unique")
    return value


def require_columns(
    value: pd.DataFrame,
    *,
    field_name: str,
    required_columns: Iterable[str],
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require a DataFrame to include the given columns."""

    missing = [column for column in required_columns if column not in value.columns]
    if missing:
        joined = ", ".join(missing)
        raise error_type(f"{field_name} is missing required columns: {joined}")
    return value


def require_exact_index_match(
    *,
    left: pd.Index,
    right: pd.Index,
    left_name: str,
    right_name: str,
    error_type: ValidationErrorType,
) -> None:
    """Require two indexes to be exactly equal (labels and order)."""

    if not left.equals(right):
        raise error_type(f"{left_name} must exactly match {right_name}")


def require_non_empty_string_column(
    value: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require all values in a DataFrame column to be non-empty strings."""

    column_values = value[column_name]
    if column_values.isna().any():
        raise error_type(f"{field_name}.{column_name} must not contain missing values")
    non_string_or_blank = [
        idx
        for idx, raw_value in column_values.items()
        if not isinstance(raw_value, str) or not raw_value.strip()
    ]
    if non_string_or_blank:
        raise error_type(
            f"{field_name}.{column_name} must contain non-empty string values"
        )
    return value
