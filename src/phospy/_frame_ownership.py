"""Internal pandas ownership helpers.

PhosPy ownership policy:
- internal objects own mutable DataFrame/Series state;
- public DataFrame exports are defensive snapshots;
- provenance fingerprints describe owned internal state at creation time;
- borrowed DataFrame access is private/internal only.
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


def export_dataframe(value: pd.DataFrame) -> pd.DataFrame:
    """Return a defensive public snapshot of an owned DataFrame."""

    return value.copy(deep=True)


def _borrow_dataframe(value: pd.DataFrame) -> pd.DataFrame:
    """Return internal borrowed DataFrame access without copying."""

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


def export_optional_dataframe(value: pd.DataFrame | None) -> pd.DataFrame | None:
    if value is None:
        return None
    return export_dataframe(value)


def _borrow_optional_dataframe(value: pd.DataFrame | None) -> pd.DataFrame | None:
    """Return internal borrowed optional DataFrame access without copying."""

    return value


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


def export_series(value: pd.Series) -> pd.Series:
    """Return a defensive public snapshot of an owned Series."""

    return value.copy(deep=True)


def _borrow_series(value: pd.Series) -> pd.Series:
    """Return internal borrowed Series access without copying."""

    return value


def own_optional_series(
    value: object | None,
    *,
    field_name: str,
    error_type: ExceptionType = TypeError,
    assume_owned: bool = False,
) -> pd.Series | None:
    if value is None:
        return None
    return own_series(
        value,
        field_name=field_name,
        error_type=error_type,
        assume_owned=assume_owned,
    )


def export_optional_series(value: pd.Series | None) -> pd.Series | None:
    if value is None:
        return None
    return export_series(value)


def _borrow_optional_series(value: pd.Series | None) -> pd.Series | None:
    """Return internal borrowed optional Series access without copying."""

    return value
