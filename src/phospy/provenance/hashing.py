"""Deterministic hashing helpers for provenance fingerprints."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from datetime import date, datetime, time
from decimal import Decimal
from typing import Any

import numpy as np
import pandas as pd

from phospy.provenance.models import TableFingerprint

DEFAULT_TABLE_HASH_ALGORITHM = "sha256"
_MISSING_SENTINEL = "<MISSING>"


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
    _update(hasher, None if table.index.name is None else str(table.index.name))
    _update(hasher, [str(label) for label in table.columns.tolist()])
    _update(hasher, [str(dtype) for dtype in table.dtypes.tolist()])
    for label in table.index.tolist():
        _update(hasher, _normalize_value(label))
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


def _normalize_value(value: object) -> object:
    if value is None or _is_missing_scalar(value):
        return {"kind": "missing", "value": _MISSING_SENTINEL}
    if isinstance(value, bool):
        return {"kind": "bool", "value": bool(value)}
    if isinstance(value, (np.integer, int)):
        return {"kind": "int", "value": int(value)}
    if isinstance(value, (np.floating, float)):
        numeric = float(value)
        if math.isnan(numeric):
            return {"kind": "missing", "value": _MISSING_SENTINEL}
        if math.isinf(numeric):
            return {
                "kind": "float",
                "value": "Infinity" if numeric > 0.0 else "-Infinity",
            }
        return {"kind": "float", "value": format(numeric, ".17g")}
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
        normalized_pairs = [
            (str(key), _normalize_value(item))
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        ]
        return {"kind": "mapping", "value": normalized_pairs}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return {"kind": "sequence", "value": [_normalize_value(item) for item in value]}
    return {"kind": "repr", "value": repr(value)}


def _is_missing_scalar(value: object) -> bool:
    missing = pd.isna(value)
    return isinstance(missing, (bool, np.bool_)) and bool(missing)


__all__ = [
    "DEFAULT_TABLE_HASH_ALGORITHM",
    "fingerprint_optional_table",
    "fingerprint_table",
    "hash_table",
]
