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

from phospy.provenance.immutability import thaw_json_value
from phospy.provenance.models import JsonValue, TableFingerprint

DEFAULT_TABLE_HASH_ALGORITHM = "sha256"
DEFAULT_STABLE_JSON_HASH_ALGORITHM = "sha256-stable-json-v1"
DEFAULT_EXACT_TABLE_HASH_ALGORITHM = "sha256-stable-json-v1"
DEFAULT_TOLERANCE_TABLE_HASH_ALGORITHM = "sha256-float-round-8dp-v1"
_MISSING_SENTINEL = "<MISSING>"
_FLOAT_HASH_DECIMAL_PLACES = 8
_SITE_METADATA_PROVENANCE_IGNORED_COLUMNS = ("display_id", "site_key")
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


def hash_table_exact(
    table: pd.DataFrame,
    *,
    name: str,
    algorithm: str = DEFAULT_TABLE_HASH_ALGORITHM,
) -> str:
    """Return an exact deterministic digest for a table and full structure."""

    return _hash_table(table, name=name, algorithm=algorithm, round_floats=False)


def hash_table_tolerance(
    table: pd.DataFrame,
    *,
    name: str,
    algorithm: str = DEFAULT_TABLE_HASH_ALGORITHM,
) -> str:
    """Return a tolerance-oriented deterministic digest for a table."""

    return _hash_table(table, name=name, algorithm=algorithm, round_floats=True)


def _hash_table(
    table: pd.DataFrame,
    *,
    name: str,
    algorithm: str,
    round_floats: bool,
) -> str:
    """Return a deterministic digest for a table and its full structure."""

    hasher = hashlib.new(algorithm)
    _update(hasher, name)
    _update(hasher, [int(table.shape[0]), int(table.shape[1])])
    _update(hasher, _index_structure(table.index, round_floats=round_floats))
    _update(hasher, _index_structure(table.columns, round_floats=round_floats))
    _update(hasher, [str(dtype) for dtype in table.dtypes.tolist()])
    values = table.to_numpy(dtype=object, copy=False)
    for row in values:
        for value in row:
            _update(hasher, _normalize_value(value, round_floats=round_floats))
    return hasher.hexdigest()


def fingerprint_table(
    table: pd.DataFrame,
    *,
    name: str,
    algorithm: str = DEFAULT_TABLE_HASH_ALGORITHM,
) -> TableFingerprint:
    """Build a typed deterministic table fingerprint."""

    normalized_table = _normalize_provenance_table_view(table, name=name)
    exact_hash_algorithm = _exact_hash_algorithm_name(algorithm)
    tolerance_hash_algorithm = _tolerance_hash_algorithm_name(algorithm)
    exact_hash = hash_table_exact(normalized_table, name=name, algorithm=algorithm)
    tolerance_hash = hash_table_tolerance(
        normalized_table, name=name, algorithm=algorithm
    )

    return TableFingerprint(
        name=name,
        rows=int(normalized_table.shape[0]),
        columns=int(normalized_table.shape[1]),
        index_name=(
            None
            if normalized_table.index.name is None
            else str(normalized_table.index.name)
        ),
        column_names=tuple(str(label) for label in normalized_table.columns.tolist()),
        dtypes=tuple(str(dtype) for dtype in normalized_table.dtypes.tolist()),
        exact_hash_algorithm=exact_hash_algorithm,
        exact_hash_value=exact_hash,
        tolerance_hash_algorithm=tolerance_hash_algorithm,
        tolerance_hash_value=tolerance_hash,
        index_structure=_index_structure(normalized_table.index, round_floats=False),
        column_index_structure=_index_structure(
            normalized_table.columns, round_floats=False
        ),
    )


def fingerprint_matrix(
    matrix: pd.DataFrame,
    *,
    name: str,
    algorithm: str = DEFAULT_TABLE_HASH_ALGORITHM,
) -> TableFingerprint:
    """Build an order-sensitive deterministic matrix fingerprint."""

    return fingerprint_table(matrix, name=name, algorithm=algorithm)


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


def fingerprint_optional_matrix(
    matrix: pd.DataFrame | None,
    *,
    name: str,
    algorithm: str = DEFAULT_TABLE_HASH_ALGORITHM,
) -> TableFingerprint | None:
    """Build an optional matrix fingerprint when a matrix exists."""

    if matrix is None:
        return None
    return fingerprint_matrix(matrix, name=name, algorithm=algorithm)


def hash_json_payload(
    payload: JsonValue,
    *,
    algorithm: str = DEFAULT_TABLE_HASH_ALGORITHM,
) -> str:
    """Return a deterministic digest for a JSON-compatible payload."""

    hasher = hashlib.new(algorithm)
    _update(hasher, payload)
    return hasher.hexdigest()


def _fingerprint_optional_table_with_normalized_axes(
    table: pd.DataFrame | None,
    *,
    name: str,
    algorithm: str = DEFAULT_TABLE_HASH_ALGORITHM,
) -> TableFingerprint | None:
    """Build an optional table fingerprint after deterministic axis sorting."""

    if table is None:
        return None
    normalized_table = _normalize_table_axes_for_fingerprint(table)
    return fingerprint_table(normalized_table, name=name, algorithm=algorithm)


def _update(hasher: hashlib._Hash, payload: Any) -> None:  # type: ignore[attr-defined]
    safe_payload = thaw_json_value(payload, field_name="provenance_hash_payload")
    encoded = json.dumps(
        safe_payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    hasher.update(encoded)
    hasher.update(b"\n")


def _normalize_value(value: object, *, round_floats: bool) -> JsonValue:
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
        return {
            "kind": "float",
            "value": (
                _normalize_float_string(numeric)
                if round_floats
                else _normalize_exact_float_string(numeric)
            ),
        }
    if isinstance(value, np.floating):
        numeric = float(cast(float, value))
        if math.isnan(numeric):
            return {"kind": "missing", "value": _MISSING_SENTINEL}
        if math.isinf(numeric):
            return {
                "kind": "float",
                "value": "Infinity" if numeric > 0.0 else "-Infinity",
            }
        return {
            "kind": "float",
            "value": (
                _normalize_float_string(numeric)
                if round_floats
                else _normalize_exact_float_string(numeric)
            ),
        }
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
            [str(key), _normalize_value(item, round_floats=round_floats)]
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        ]
        return cast(
            JsonValue,
            {"kind": "mapping", "value": normalized_pairs},
        )
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return {
            "kind": "sequence",
            "value": [
                _normalize_value(item, round_floats=round_floats) for item in value
            ],
        }
    return {"kind": "repr", "value": repr(value)}


def _normalize_axis_name(value: object, *, round_floats: bool) -> JsonValue:
    if value is None:
        return {"kind": "none", "value": None}
    return _normalize_value(value, round_floats=round_floats)


def _index_structure(index: pd.Index, *, round_floats: bool) -> dict[str, JsonValue]:
    if isinstance(index, pd.MultiIndex):
        return {
            "type": "multi_index",
            "index_class": type(index).__name__,
            "nlevels": int(index.nlevels),
            "names": [
                _normalize_axis_name(name, round_floats=round_floats)
                for name in index.names
            ],
            "level_dtypes": [str(level.dtype) for level in index.levels],
            "values": [
                [
                    _normalize_value(component, round_floats=round_floats)
                    for component in tuple(value)
                ]
                for value in index.tolist()
            ],
        }
    if isinstance(index, pd.RangeIndex):
        return {
            "type": "range_index",
            "index_class": type(index).__name__,
            "name": _normalize_axis_name(index.name, round_floats=round_floats),
            "start": int(index.start),
            "stop": int(index.stop),
            "step": int(index.step),
            "dtype": str(index.dtype),
        }
    return {
        "type": "index",
        "index_class": type(index).__name__,
        "name": _normalize_axis_name(index.name, round_floats=round_floats),
        "dtype": str(index.dtype),
        "values": [
            _normalize_value(label, round_floats=round_floats)
            for label in index.tolist()
        ],
    }


def _is_missing_scalar(value: object) -> bool:
    if not isinstance(value, _PANDAS_MISSING_SCALAR_TYPES):
        return False
    if isinstance(value, float):
        return math.isnan(value)
    if isinstance(value, np.floating):
        return bool(np.isnan(value))
    if value is pd.NA or value is pd.NaT:
        return True
    if isinstance(value, np.datetime64 | np.timedelta64):
        return bool(np.isnat(value))
    return False


def _normalize_float_string(value: float) -> str:
    # Normalize floats to a fixed decimal precision so tiny parser/ULP
    # differences across environments do not churn provenance fingerprints.
    if value == 0.0:
        return "0"
    rounded = round(value, _FLOAT_HASH_DECIMAL_PLACES)
    if rounded == 0.0:
        return "0"
    return format(rounded, f".{_FLOAT_HASH_DECIMAL_PLACES}f").rstrip("0").rstrip(".")


def _normalize_exact_float_string(value: float) -> str:
    return value.hex()


def _exact_hash_algorithm_name(algorithm: str) -> str:
    if algorithm == DEFAULT_TABLE_HASH_ALGORITHM:
        return DEFAULT_EXACT_TABLE_HASH_ALGORITHM
    return f"{algorithm}-stable-json-v1"


def _tolerance_hash_algorithm_name(algorithm: str) -> str:
    if algorithm == DEFAULT_TABLE_HASH_ALGORITHM:
        return DEFAULT_TOLERANCE_TABLE_HASH_ALGORITHM
    return f"{algorithm}-float-round-{_FLOAT_HASH_DECIMAL_PLACES}dp-v1"


def _normalize_provenance_table_view(
    table: pd.DataFrame,
    *,
    name: str,
) -> pd.DataFrame:
    if name != "dataset.site_metadata":
        return table
    removable = [
        column
        for column in _SITE_METADATA_PROVENANCE_IGNORED_COLUMNS
        if column in table.columns
    ]
    if not removable:
        return table
    # Keep provenance contracts stable while builder-owned identity fields are
    # introduced. A later provenance schema ticket can make this explicit.
    return table.drop(columns=removable)


def _normalize_table_axes_for_fingerprint(table: pd.DataFrame) -> pd.DataFrame:
    normalized = table
    try:
        normalized = normalized.sort_index(axis=0, kind="mergesort")
    except Exception:
        pass
    try:
        normalized = normalized.sort_index(axis=1, kind="mergesort")
    except Exception:
        pass
    return normalized


__all__ = [
    "DEFAULT_EXACT_TABLE_HASH_ALGORITHM",
    "DEFAULT_STABLE_JSON_HASH_ALGORITHM",
    "DEFAULT_TABLE_HASH_ALGORITHM",
    "DEFAULT_TOLERANCE_TABLE_HASH_ALGORITHM",
    "fingerprint_matrix",
    "fingerprint_optional_matrix",
    "fingerprint_optional_table",
    "fingerprint_table",
    "hash_json_payload",
    "hash_table_exact",
    "hash_table_tolerance",
]
