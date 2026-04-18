"""Internal pandas ownership helpers.

Policy:
- Public boundary constructors copy caller-provided pandas objects by default.
- Internal assembly paths may transfer already-owned objects without re-copying.
- Internal DTOs should alias owned objects rather than deep-copying again.
"""

from __future__ import annotations

import pandas as pd


def own_dataframe(
    value: object,
    *,
    field_name: str,
    assume_owned: bool = False,
) -> pd.DataFrame:
    """Return an owned DataFrame, copying only when ownership is not established."""

    if not isinstance(value, pd.DataFrame):
        raise TypeError(f"{field_name} must be a pandas DataFrame")
    if assume_owned:
        return value
    return value.copy(deep=True)


def own_optional_dataframe(
    value: object | None,
    *,
    field_name: str,
    assume_owned: bool = False,
) -> pd.DataFrame | None:
    if value is None:
        return None
    return own_dataframe(value, field_name=field_name, assume_owned=assume_owned)


def own_series(
    value: object,
    *,
    field_name: str,
    assume_owned: bool = False,
) -> pd.Series:
    """Return an owned Series, copying only when ownership is not established."""

    if not isinstance(value, pd.Series):
        raise TypeError(f"{field_name} must be a pandas Series")
    if assume_owned:
        return value
    return value.copy(deep=True)
