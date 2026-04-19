"""Primitive payload coercion for bundle manifest/config decoding."""

from __future__ import annotations

from collections.abc import Mapping

from phospy.errors.input import PhosPyInputError


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
