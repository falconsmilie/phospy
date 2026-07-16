"""Shared pandas DataFrame and Series ownership helpers."""

from __future__ import annotations

from typing import TypeVar, cast

import numpy as np
import pandas as pd

ExceptionType = type[Exception]
_PandasObject = TypeVar("_PandasObject", pd.DataFrame, pd.Series)
_PANDAS_MAJOR_VERSION = int(str(pd.__version__).split(".", maxsplit=1)[0])


def _pandas_has_native_copy_on_write() -> bool:
    """Return whether shallow copies are locally mutation-isolated by pandas."""

    return _PANDAS_MAJOR_VERSION >= 3


def _mark_numpy_blocks_read_only(value: pd.DataFrame | pd.Series) -> bool:
    """Mark a shallow pandas copy's NumPy blocks read-only.

    This intentionally uses pandas' private BlockManager surface in one
    contained helper. If that surface is unavailable, or if a block is backed by
    an extension array that cannot be made read-only this way, callers fall back
    to a deep copy.
    """

    manager = getattr(value, "_mgr", None)
    blocks = getattr(manager, "blocks", None)
    if blocks is None:
        return False

    for block in blocks:
        values = getattr(block, "values", None)
        if not isinstance(values, np.ndarray):
            return False
        read_only_values = values.view()
        read_only_values.flags.writeable = False
        try:
            block.values = read_only_values
        except (AttributeError, TypeError, ValueError):
            return False
    return True


def _borrow_pandas_object(value: _PandasObject) -> _PandasObject:
    """Return a mutation-isolated internal snapshot without global mutation."""

    borrowed = cast(_PandasObject, value.copy(deep=False))
    if _pandas_has_native_copy_on_write():
        return borrowed
    if _mark_numpy_blocks_read_only(borrowed):
        return borrowed
    return cast(_PandasObject, value.copy(deep=True))


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


def borrow_dataframe(value: pd.DataFrame) -> pd.DataFrame:
    """Return internal borrowed DataFrame access without deep-copy churn."""

    return _borrow_dataframe(value)


def _borrow_dataframe(value: pd.DataFrame) -> pd.DataFrame:
    """Return internal borrowed DataFrame access without deep-copy churn.

    Borrowed access is mutation-isolated from the owning frame:
    - pandas>=3: shallow copy uses native copy-on-write semantics.
    - NumPy-backed pandas<3 frames: shallow copy with read-only borrowed blocks.
    - unsupported pandas internals: deep-copy fallback.

    Internal mutation that should affect owned scientific state must happen on
    explicitly owned frames, never through `_borrow_*` accessors. Writes to a
    borrowed object may raise or detach locally; they must not mutate the owner.
    """

    if not isinstance(value, pd.DataFrame):
        raise TypeError("borrowed frame access requires a pandas DataFrame")

    return _borrow_pandas_object(value)


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


def borrow_optional_dataframe(value: pd.DataFrame | None) -> pd.DataFrame | None:
    """Return internal borrowed optional DataFrame access."""

    return _borrow_optional_dataframe(value)


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

    return _borrow_pandas_object(value)


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
    "export_dataframe",
    "export_optional_dataframe",
    "export_optional_series",
    "export_series",
    "own_dataframe",
    "own_optional_dataframe",
    "own_optional_series",
    "own_series",
]
