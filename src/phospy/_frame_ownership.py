"""Internal pandas ownership helpers.

PhosPy treats DataFrames as owned mutable state internally.

Input DataFrames are copied when accepted into validated dataset/table/result
objects. Internal workflow code may pass owned DataFrames between private
helpers without repeated defensive copies.

Public accessors should make ownership explicit:
- safe/public access returns copied snapshots;
- borrowed access may return the owned DataFrame directly and is unsafe to
  mutate unless intentional owner mutation is desired.

Provenance fingerprints describe the owned internal state at creation time.
"""

from __future__ import annotations

import pandas as pd

ExceptionType = type[Exception]


def own_dataframe(
    value: object,
    *,
    field_name: str,
    error_type: ExceptionType = TypeError,
    assume_owned: bool = False,
) -> pd.DataFrame:
    """Return an owned DataFrame, copying only when ownership is not established."""

    if not isinstance(value, pd.DataFrame):
        raise error_type(f"{field_name} must be a pandas DataFrame")
    if assume_owned:
        return value
    return value.copy(deep=True)


def export_dataframe(value: pd.DataFrame, *, copy: bool = True) -> pd.DataFrame:
    """Return a public DataFrame export.

    When ``copy=True`` (default), callers receive a safe snapshot copy.
    When ``copy=False``, callers receive a borrowed reference to owned internal
    state and mutating it mutates the owning object.
    """

    if copy:
        return value.copy(deep=True)
    return value


def own_optional_dataframe(
    value: object | None,
    *,
    field_name: str,
    error_type: ExceptionType = TypeError,
    assume_owned: bool = False,
) -> pd.DataFrame | None:
    if value is None:
        return None
    return own_dataframe(
        value,
        field_name=field_name,
        error_type=error_type,
        assume_owned=assume_owned,
    )


def export_optional_dataframe(
    value: pd.DataFrame | None, *, copy: bool = True
) -> pd.DataFrame | None:
    if value is None:
        return None
    return export_dataframe(value, copy=copy)


def own_series(
    value: object,
    *,
    field_name: str,
    error_type: ExceptionType = TypeError,
    assume_owned: bool = False,
) -> pd.Series:
    """Return an owned Series, copying only when ownership is not established."""

    if not isinstance(value, pd.Series):
        raise error_type(f"{field_name} must be a pandas Series")
    if assume_owned:
        return value
    return value.copy(deep=True)
