"""Internal pandas ownership helpers.

PhosPy frame ownership policy:
- `_frame_ownership` and dataset/result models own DataFrame/Series policy.
- public properties/export helpers always return defensive copies.
- internal borrowed access is explicit (`_borrow_*`) and internal-only.
- borrowed frames are mutation-isolated snapshot views for internal read paths.
- mutation is allowed only in owned construction/transformation code paths.
"""

from __future__ import annotations

import pandas as pd

ExceptionType = type[Exception]
_PANDAS_MAJOR_VERSION = int(str(pd.__version__).split(".", maxsplit=1)[0])
_COPY_ON_WRITE_ENSURED = False


def _ensure_copy_on_write_enabled() -> None:
    """Enable pandas copy-on-write where available for safe shallow borrows."""

    global _COPY_ON_WRITE_ENSURED
    if _COPY_ON_WRITE_ENSURED:
        return
    if _PANDAS_MAJOR_VERSION >= 3:
        _COPY_ON_WRITE_ENSURED = True
        return
    mode_options = getattr(pd.options, "mode", None)
    if mode_options is not None and hasattr(mode_options, "copy_on_write"):
        pd.options.mode.copy_on_write = True
    _COPY_ON_WRITE_ENSURED = True


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
    """Return internal borrowed DataFrame access without deep-copy churn.

    Borrowed access is mutation-isolated from the owning frame:
    - pandas>=3: shallow copy uses native copy-on-write semantics.
    - pandas<3: this helper enables copy-on-write mode before borrowing.

    Internal mutation that should affect owned scientific state must happen on
    explicitly owned frames, never through `_borrow_*` accessors.
    """

    if not isinstance(value, pd.DataFrame):
        raise TypeError("borrowed frame access requires a pandas DataFrame")

    _ensure_copy_on_write_enabled()
    return value.copy(deep=False)


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
    """Return internal borrowed optional DataFrame access."""

    if value is not None and not isinstance(value, pd.DataFrame):
        raise TypeError(
            "borrowed optional frame access requires a pandas DataFrame or None"
        )

    if value is None:
        return None
    return _borrow_dataframe(value)


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
    """Return internal borrowed Series access."""

    if not isinstance(value, pd.Series):
        raise TypeError("borrowed series access requires a pandas Series")

    _ensure_copy_on_write_enabled()
    return value.copy(deep=False)


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
    """Return internal borrowed optional Series access."""

    if value is not None and not isinstance(value, pd.Series):
        raise TypeError(
            "borrowed optional series access requires a pandas Series or None"
        )

    if value is None:
        return None
    return _borrow_series(value)


__all__ = [
    "ExceptionType",
    "export_dataframe",
    "export_optional_dataframe",
    "export_optional_series",
    "export_series",
    "own_dataframe",
    "own_optional_dataframe",
    "own_optional_series",
    "own_series",
]
