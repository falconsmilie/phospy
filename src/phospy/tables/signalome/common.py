from __future__ import annotations

import json

import numpy as np
import pandas as pd

from phospy.tables.base import ValidationErrorType


def _column_series(frame: pd.DataFrame, column_name: str) -> pd.Series:
    return frame.loc[:, column_name]


def _numeric_series(frame: pd.DataFrame, column_name: str) -> pd.Series:
    return pd.to_numeric(_column_series(frame, column_name), errors="coerce")


def _require_string_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: ValidationErrorType,
) -> None:
    values = _column_series(frame, column_name)
    if values.isna().any():
        raise error_type(f"{field_name}.{column_name} must not contain missing values")
    if not all(isinstance(value, str) for value in values.tolist()):
        raise error_type(f"{field_name}.{column_name} must contain string values")


def _require_boolean_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: ValidationErrorType,
) -> None:
    values = _column_series(frame, column_name)
    if values.isna().any():
        raise error_type(f"{field_name}.{column_name} must not contain missing values")
    invalid = [
        value for value in values.tolist() if not isinstance(value, (bool, np.bool_))
    ]
    if invalid:
        raise error_type(f"{field_name}.{column_name} must contain boolean values")


def _require_integer_compatible_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: ValidationErrorType,
    allow_missing: bool,
) -> None:
    numeric = _numeric_series(frame, column_name)
    if not allow_missing and numeric.isna().any():
        raise error_type(f"{field_name}.{column_name} must not contain missing values")
    finite_values = numeric.dropna().to_numpy(dtype="float64", copy=False)
    if not np.isfinite(finite_values).all():
        raise error_type(
            f"{field_name}.{column_name} must contain finite integer-compatible values"
        )
    if not np.isclose(finite_values, np.round(finite_values)).all():
        raise error_type(
            f"{field_name}.{column_name} must contain integer-compatible values"
        )


def _require_numeric_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: ValidationErrorType,
    allow_missing: bool,
) -> None:
    numeric = _numeric_series(frame, column_name)
    if not allow_missing and numeric.isna().any():
        raise error_type(f"{field_name}.{column_name} must not contain missing values")
    finite_values = numeric.dropna().to_numpy(dtype="float64", copy=False)
    if not np.isfinite(finite_values).all():
        raise error_type(
            f"{field_name}.{column_name} must contain finite numeric values"
        )


def _require_json_string_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: ValidationErrorType,
) -> None:
    for value in frame.loc[:, column_name].tolist():
        if not isinstance(value, str):
            raise error_type(
                f"{field_name}.{column_name} must contain JSON-encoded strings"
            )
        try:
            json.loads(value)
        except json.JSONDecodeError as exc:
            raise error_type(
                f"{field_name}.{column_name} must contain parseable JSON strings"
            ) from exc


def _require_non_negative_integer_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: ValidationErrorType,
) -> None:
    _require_integer_compatible_column(
        frame,
        field_name=field_name,
        column_name=column_name,
        error_type=error_type,
        allow_missing=False,
    )
    numeric = _numeric_series(frame, column_name)
    values = numeric.to_numpy(dtype="float64", copy=False)
    if (values < 0.0).any():
        raise error_type(
            f"{field_name}.{column_name} must contain non-negative integer values"
        )


def _require_integer_compatible_index(
    index: pd.Index,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> None:
    numeric = pd.to_numeric(index.to_series(index=index), errors="coerce")
    if numeric.isna().any():
        raise error_type(f"{field_name} must contain integer-compatible labels")
    values = numeric.to_numpy(dtype="float64", copy=False)
    if not np.isfinite(values).all():
        raise error_type(f"{field_name} must contain finite integer-compatible labels")
    if not np.isclose(values, np.round(values)).all():
        raise error_type(f"{field_name} must contain integer-compatible labels")


def _require_non_negative_integer_index(
    index: pd.Index,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> None:
    numeric = pd.to_numeric(index.to_series(index=index), errors="coerce")
    values = numeric.to_numpy(dtype="float64", copy=False)
    if (values < 0.0).any():
        raise error_type(f"{field_name} must contain non-negative integer labels")


def _require_numeric_bounds(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: ValidationErrorType,
    minimum: float,
    maximum: float,
    allow_missing: bool,
) -> None:
    numeric = _numeric_series(frame, column_name)
    if not allow_missing and numeric.isna().any():
        raise error_type(f"{field_name}.{column_name} must not contain missing values")
    values = numeric.dropna().to_numpy(dtype="float64", copy=False)
    if ((values < float(minimum)) | (values > float(maximum))).any():
        raise error_type(
            f"{field_name}.{column_name} must be between {minimum:.1f} and {maximum:.1f}"
        )
