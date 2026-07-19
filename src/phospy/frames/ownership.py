"""Shared pandas DataFrame and Series ownership helpers.

Ownership policy:

- public constructors copy caller-provided DataFrames and Series before storing
  them, unless an internal caller explicitly passes an already-owned object;
- public exports always return detached snapshots;
- internal borrow helpers also return detached snapshots, not owner-backed
  read-only aliases.

Pandas ``deep=True`` does not recursively copy mutable Python objects held in
object-dtype cells. These helpers therefore isolate supported mutable object
cells explicitly: ``list``, ``dict``, ``set``, ``tuple``/``frozenset``
containers, and ``numpy.ndarray`` values. Unsupported mutable object cells are
rejected with a clear error instead of being stored or exported as aliases.
"""

from __future__ import annotations

from collections.abc import MutableMapping, MutableSequence, MutableSet
from datetime import date, datetime, time, timedelta
from decimal import Decimal
from enum import Enum
from fractions import Fraction
from pathlib import PurePath
from typing import Any, TypeVar, cast
from uuid import UUID

import numpy as np
import pandas as pd

ExceptionType = type[Exception]
_PandasObject = TypeVar("_PandasObject", pd.DataFrame, pd.Series)
_IMMUTABLE_OBJECT_TYPES = (
    bool,
    int,
    float,
    complex,
    str,
    bytes,
    type(None),
    range,
    slice,
    Decimal,
    Fraction,
    date,
    datetime,
    time,
    timedelta,
    PurePath,
    UUID,
    Enum,
    np.generic,
    pd.Timestamp,
    pd.Timedelta,
    pd.Interval,
    pd.Period,
)


def _copy_pandas_object(
    value: _PandasObject,
    *,
    field_name: str,
    error_type: ExceptionType,
) -> _PandasObject:
    copied = cast(_PandasObject, value.copy(deep=True))
    _isolate_object_cells(copied, field_name=field_name, error_type=error_type)
    return copied


def _isolate_object_cells(
    value: pd.DataFrame | pd.Series,
    *,
    field_name: str,
    error_type: ExceptionType,
) -> None:
    memo: dict[int, object] = {}
    active: set[int] = set()
    if isinstance(value, pd.Series):
        if not pd.api.types.is_object_dtype(value.dtype):
            return
        for position in range(len(value.index)):
            location = f"position {position}, index {value.index[position]!r}"
            isolated = _isolate_object_cell_value(
                value.iat[position],
                field_name=field_name,
                location=location,
                error_type=error_type,
                memo=memo,
                active=active,
            )
            cast(Any, value.iat)[position] = isolated
        return

    for column_position, dtype in enumerate(value.dtypes):
        if not pd.api.types.is_object_dtype(dtype):
            continue
        column_label = value.columns[column_position]
        for row_position in range(len(value.index)):
            location = f"row {value.index[row_position]!r}, column {column_label!r}"
            isolated = _isolate_object_cell_value(
                value.iat[row_position, column_position],
                field_name=field_name,
                location=location,
                error_type=error_type,
                memo=memo,
                active=active,
            )
            cast(Any, value.iat)[row_position, column_position] = isolated


def _isolate_object_cell_value(
    value: object,
    *,
    field_name: str,
    location: str,
    error_type: ExceptionType,
    memo: dict[int, object],
    active: set[int],
) -> object:
    value_id = id(value)
    if value_id in memo:
        return memo[value_id]
    if value_id in active:
        raise error_type(
            f"{field_name} contains circular mutable object data at {location}; "
            "object cells must be acyclic to be safely isolated"
        )
    if _is_known_immutable_object(value):
        return value
    if isinstance(value, np.ndarray):
        return _copy_numpy_object_cell(
            value,
            field_name=field_name,
            location=location,
            error_type=error_type,
            memo=memo,
            active=active,
        )
    if isinstance(value, list):
        copied_list: list[object] = []
        memo[value_id] = copied_list
        active.add(value_id)
        try:
            copied_list.extend(
                _isolate_object_cell_value(
                    item,
                    field_name=field_name,
                    location=f"{location}[]",
                    error_type=error_type,
                    memo=memo,
                    active=active,
                )
                for item in value
            )
        finally:
            active.discard(value_id)
        return copied_list
    if isinstance(value, dict):
        copied_dict: dict[object, object] = {}
        memo[value_id] = copied_dict
        active.add(value_id)
        try:
            for key, item in value.items():
                copied_key = _isolate_object_cell_value(
                    key,
                    field_name=field_name,
                    location=f"{location}.<key>",
                    error_type=error_type,
                    memo=memo,
                    active=active,
                )
                copied_item = _isolate_object_cell_value(
                    item,
                    field_name=field_name,
                    location=f"{location}[{key!r}]",
                    error_type=error_type,
                    memo=memo,
                    active=active,
                )
                try:
                    copied_dict[copied_key] = copied_item
                except TypeError as exc:
                    raise error_type(
                        f"{field_name} contains uncopyable dict key at {location}: "
                        f"{type(key).__module__}.{type(key).__qualname__}"
                    ) from exc
        finally:
            active.discard(value_id)
        return copied_dict
    if isinstance(value, set):
        copied_set: set[object] = set()
        memo[value_id] = copied_set
        active.add(value_id)
        try:
            for item in value:
                copied_item = _isolate_object_cell_value(
                    item,
                    field_name=field_name,
                    location=f"{location}.<set-item>",
                    error_type=error_type,
                    memo=memo,
                    active=active,
                )
                try:
                    copied_set.add(copied_item)
                except TypeError as exc:
                    raise error_type(
                        f"{field_name} contains uncopyable set item at {location}: "
                        f"{type(item).__module__}.{type(item).__qualname__}"
                    ) from exc
        finally:
            active.discard(value_id)
        return copied_set
    if isinstance(value, tuple):
        active.add(value_id)
        try:
            copied_tuple = tuple(
                _isolate_object_cell_value(
                    item,
                    field_name=field_name,
                    location=f"{location}[]",
                    error_type=error_type,
                    memo=memo,
                    active=active,
                )
                for item in value
            )
        finally:
            active.discard(value_id)
        memo[value_id] = copied_tuple
        return copied_tuple
    if isinstance(value, frozenset):
        active.add(value_id)
        try:
            copied_items = [
                _isolate_object_cell_value(
                    item,
                    field_name=field_name,
                    location=f"{location}.<frozenset-item>",
                    error_type=error_type,
                    memo=memo,
                    active=active,
                )
                for item in value
            ]
            copied_frozenset = frozenset(copied_items)
        except TypeError as exc:
            raise error_type(
                f"{field_name} contains uncopyable frozenset item at {location}"
            ) from exc
        finally:
            active.discard(value_id)
        memo[value_id] = copied_frozenset
        return copied_frozenset
    if _looks_like_unsupported_mutable_object(value):
        raise error_type(
            f"{field_name} contains unsupported mutable object at {location}: "
            f"{type(value).__module__}.{type(value).__qualname__}; supported "
            "mutable object cells are list, dict, set, tuple/frozenset "
            "containers, and numpy.ndarray"
        )
    return value


def _copy_numpy_object_cell(
    value: np.ndarray,
    *,
    field_name: str,
    location: str,
    error_type: ExceptionType,
    memo: dict[int, object],
    active: set[int],
) -> np.ndarray:
    copied = value.copy()
    memo[id(value)] = copied
    if value.dtype != object:
        return copied

    active.add(id(value))
    try:
        source_flat = value.reshape(-1)
        copied_flat = copied.reshape(-1)
        for position, item in enumerate(source_flat):
            copied_flat[position] = _isolate_object_cell_value(
                item,
                field_name=field_name,
                location=f"{location}.ndarray[{position}]",
                error_type=error_type,
                memo=memo,
                active=active,
            )
    finally:
        active.discard(id(value))
    return copied


def _is_known_immutable_object(value: object) -> bool:
    if value is pd.NA or value is pd.NaT:
        return True
    return isinstance(value, _IMMUTABLE_OBJECT_TYPES)


def _looks_like_unsupported_mutable_object(value: object) -> bool:
    if isinstance(value, MutableMapping | MutableSequence | MutableSet):
        return True
    if isinstance(value, bytearray | memoryview):
        return True
    if hasattr(value, "__dict__"):
        return True
    slots = getattr(type(value), "__slots__", ())
    return bool(slots)


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
    return _copy_pandas_object(
        value,
        field_name=field_name,
        error_type=error_type,
    )


def export_dataframe(value: pd.DataFrame) -> pd.DataFrame:
    """Return a defensive public snapshot of an owned DataFrame."""

    return _copy_pandas_object(
        value,
        field_name="public DataFrame export",
        error_type=TypeError,
    )


def borrow_dataframe(value: pd.DataFrame) -> pd.DataFrame:
    """Return internal borrowed DataFrame access without deep-copy churn."""

    return _borrow_dataframe(value)


def _borrow_dataframe(value: pd.DataFrame) -> pd.DataFrame:
    """Return detached internal DataFrame snapshot access.

    Internal mutation that should affect owned scientific state must happen on an
    explicitly owned frame, never through `_borrow_*` accessors. Restoring
    writeability on arrays exposed by the returned snapshot must not mutate the
    owner.
    """

    if not isinstance(value, pd.DataFrame):
        raise TypeError("borrowed frame access requires a pandas DataFrame")

    return _copy_pandas_object(
        value,
        field_name="borrowed DataFrame snapshot",
        error_type=TypeError,
    )


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
    return _copy_pandas_object(
        value,
        field_name=field_name,
        error_type=error_type,
    )


def export_series(value: pd.Series) -> pd.Series:
    """Return a defensive public snapshot of an owned Series."""

    return _copy_pandas_object(
        value,
        field_name="public Series export",
        error_type=TypeError,
    )


def _borrow_series(value: pd.Series) -> pd.Series:
    """Return detached internal Series snapshot access."""

    if not isinstance(value, pd.Series):
        raise TypeError("borrowed series access requires a pandas Series")

    return _copy_pandas_object(
        value,
        field_name="borrowed Series snapshot",
        error_type=TypeError,
    )


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
