"""Primitive payload coercion for bundle manifest/config decoding."""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import TypeAlias

from phospy.errors.input import PhosPyInputError

JsonPrimitive: TypeAlias = None | str | bool | int | float
JsonValue: TypeAlias = JsonPrimitive | list["JsonValue"] | dict[str, "JsonValue"]


def require_mapping(value: object, *, field_name: str) -> Mapping[str, object]:
    """Require a mapping payload field."""

    if isinstance(value, Mapping):
        return value
    raise PhosPyInputError(f"{field_name} must be an object")


def require_str(value: object, *, field_name: str) -> str:
    """Require a non-empty string payload field."""

    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a string")
    normalized = value.strip()
    if not normalized:
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return normalized


def require_bool(value: object, *, field_name: str) -> bool:
    """Require a bool payload field."""

    if not isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be a bool")
    return value


def require_int(value: object, *, field_name: str) -> int:
    """Require an int payload field (excluding bool)."""

    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(f"{field_name} must be an int")
    return value


def require_float(value: object, *, field_name: str) -> float:
    """Require a float payload field (accepting ints, excluding bool)."""

    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise PhosPyInputError(f"{field_name} must be a float")
    return float(value)


def validate_json_safe_mapping(
    value: object,
    *,
    field_name: str,
) -> dict[str, JsonValue]:
    """Validate and normalize a mapping as stable JSON-safe metadata."""

    payload = require_mapping(value, field_name=field_name)
    validated = _validate_json_safe_value(payload, path=field_name)
    if not isinstance(validated, dict):
        raise PhosPyInputError(f"{field_name} must be a JSON object")
    return validated


def _validate_json_safe_value(value: object, *, path: str) -> JsonValue:
    if value is None:
        return None
    if isinstance(value, str):
        return value
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise PhosPyInputError(f"{path} must contain only finite float values")
        return value
    if isinstance(value, list):
        return [
            _validate_json_safe_value(item, path=f"{path}[{index}]")
            for index, item in enumerate(value)
        ]
    if isinstance(value, tuple):
        return [
            _validate_json_safe_value(item, path=f"{path}[{index}]")
            for index, item in enumerate(value)
        ]
    if isinstance(value, Mapping):
        normalized: dict[str, JsonValue] = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise PhosPyInputError(
                    f"{path} must contain only string keys; got key "
                    f"{key!r} ({type(key).__name__})"
                )
            normalized[key] = _validate_json_safe_value(item, path=f"{path}.{key}")
        return normalized
    raise PhosPyInputError(
        f"{path} contains unsupported value type "
        f"{type(value).__module__}.{type(value).__name__}; expected JSON-safe "
        "scalars, arrays, or objects"
    )
