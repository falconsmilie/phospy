"""Deterministic hashing helpers for provenance fingerprints."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from datetime import date, datetime, time
from decimal import Decimal
from typing import Any, cast

import numpy as np
import pandas as pd

from phospy.provenance.models import JsonValue, TableFingerprint

DEFAULT_TABLE_HASH_ALGORITHM = "sha256"
_MISSING_SENTINEL = "<MISSING>"
_PANDAS_MISSING_SCALAR_TYPES = (
    str,
    bytes,
    bool,
    int,
    float,
    complex,
    Decimal,
    date,
    datetime,
    time,
    np.generic,
    np.datetime64,
    np.timedelta64,
)


def hash_table(
    table: pd.DataFrame,
    *,
    name: str,
    algorithm: str = DEFAULT_TABLE_HASH_ALGORITHM,
) -> str:
    """Return a deterministic digest for a table and its full structure."""

    hasher = hashlib.new(algorithm)
    _update(hasher, name)
    _update(hasher, [int(table.shape[0]), int(table.shape[1])])
    _update(hasher, _index_structure(table.index))
    _update(hasher, _index_structure(table.columns))
    _update(hasher, [str(dtype) for dtype in table.dtypes.tolist()])
    values = table.to_numpy(dtype=object, copy=False)
    for row in values:
        for value in row:
            _update(hasher, _normalize_value(value))
    return hasher.hexdigest()


def fingerprint_table(
    table: pd.DataFrame,
    *,
    name: str,
    algorithm: str = DEFAULT_TABLE_HASH_ALGORITHM,
) -> TableFingerprint:
    """Build a typed deterministic table fingerprint."""

    return TableFingerprint(
        name=name,
        rows=int(table.shape[0]),
        columns=int(table.shape[1]),
        index_name=None if table.index.name is None else str(table.index.name),
        column_names=tuple(str(label) for label in table.columns.tolist()),
        dtypes=tuple(str(dtype) for dtype in table.dtypes.tolist()),
        hash_algorithm=algorithm,
        hash_value=hash_table(table, name=name, algorithm=algorithm),
        index_structure=_index_structure(table.index),
        column_index_structure=_index_structure(table.columns),
    )


def fingerprint_optional_table(
    table: pd.DataFrame | None,
    *,
    name: str,
    algorithm: str = DEFAULT_TABLE_HASH_ALGORITHM,
) -> TableFingerprint | None:
    """Build an optional table fingerprint when a table exists."""

    if table is None:
        return None
    return fingerprint_table(table, name=name, algorithm=algorithm)


def _update(hasher: hashlib._Hash, payload: Any) -> None:  # type: ignore[attr-defined]
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    hasher.update(encoded)
    hasher.update(b"\n")


def _normalize_value(value: object) -> JsonValue:
    if value is None or _is_missing_scalar(value):
        return {"kind": "missing", "value": _MISSING_SENTINEL}
    if isinstance(value, bool):
        return {"kind": "bool", "value": bool(value)}
    if isinstance(value, int):
        return {"kind": "int", "value": int(value)}
    if isinstance(value, np.integer):
        return {"kind": "int", "value": int(cast(int, value))}
    if isinstance(value, float):
        numeric = float(value)
        if math.isnan(numeric):
            return {"kind": "missing", "value": _MISSING_SENTINEL}
        if math.isinf(numeric):
            return {
                "kind": "float",
                "value": "Infinity" if numeric > 0.0 else "-Infinity",
            }
        return {"kind": "float", "value": _normalize_float_string(numeric)}
    if isinstance(value, np.floating):
        numeric = float(cast(float, value))
        if math.isnan(numeric):
            return {"kind": "missing", "value": _MISSING_SENTINEL}
        if math.isinf(numeric):
            return {
                "kind": "float",
                "value": "Infinity" if numeric > 0.0 else "-Infinity",
            }
        return {"kind": "float", "value": _normalize_float_string(numeric)}
    if isinstance(value, Decimal):
        return {"kind": "decimal", "value": format(value, "f")}
    if isinstance(value, str):
        return {"kind": "str", "value": value}
    if isinstance(value, bytes):
        return {"kind": "bytes", "value": value.decode("utf-8", errors="replace")}
    if isinstance(value, (datetime, date, time)):
        return {"kind": "datetime", "value": value.isoformat()}
    if isinstance(value, np.datetime64):
        return {"kind": "datetime64", "value": pd.Timestamp(value).isoformat()}
    if isinstance(value, np.timedelta64):
        return {"kind": "timedelta64", "value": str(pd.Timedelta(value))}
    if isinstance(value, Mapping):
        normalized_pairs: list[list[JsonValue]] = [
            [str(key), _normalize_value(item)]
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        ]
        return cast(
            JsonValue,
            {"kind": "mapping", "value": normalized_pairs},
        )
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return {"kind": "sequence", "value": [_normalize_value(item) for item in value]}
    return {"kind": "repr", "value": repr(value)}


def _normalize_axis_name(value: object) -> JsonValue:
    if value is None:
        return {"kind": "none", "value": None}
    return _normalize_value(value)


def _index_structure(index: pd.Index) -> dict[str, JsonValue]:
    if isinstance(index, pd.MultiIndex):
        return {
            "type": "multi_index",
            "index_class": type(index).__name__,
            "nlevels": int(index.nlevels),
            "names": [_normalize_axis_name(name) for name in index.names],
            "level_dtypes": [str(level.dtype) for level in index.levels],
            "values": [
                [_normalize_value(component) for component in tuple(value)]
                for value in index.tolist()
            ],
        }
    if isinstance(index, pd.RangeIndex):
        return {
            "type": "range_index",
            "index_class": type(index).__name__,
            "name": _normalize_axis_name(index.name),
            "start": int(index.start),
            "stop": int(index.stop),
            "step": int(index.step),
            "dtype": str(index.dtype),
        }
    return {
        "type": "index",
        "index_class": type(index).__name__,
        "name": _normalize_axis_name(index.name),
        "dtype": str(index.dtype),
        "values": [_normalize_value(label) for label in index.tolist()],
    }


def _is_missing_scalar(value: object) -> bool:
    if not isinstance(value, _PANDAS_MISSING_SCALAR_TYPES):
        return False
    try:
        missing = pd.isna(value)
    except Exception:
        return False
    return isinstance(missing, (bool, np.bool_)) and bool(missing)


def _normalize_float_string(value: float) -> str:
    # Hex form is an exact IEEE-754 representation and remains stable across
    # Python versions/platforms for the same float payload.
    if value == 0.0:
        return "0"
    return value.hex()


__all__ = [
    "DEFAULT_TABLE_HASH_ALGORITHM",
    "fingerprint_optional_table",
    "fingerprint_table",
    "hash_table",
]
