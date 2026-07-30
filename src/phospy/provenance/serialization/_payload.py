"""Shared provenance payload coercion helpers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import cast

from phospy.errors.input import PhosPyInputError
from phospy.provenance.immutability import thaw_json_value
from phospy.provenance.models import JsonValue

_LEGACY_PROVENANCE_SCHEMA_ERROR = (
    "Legacy provenance schemas are no longer supported. Regenerate the result "
    "with the current PhosPy version."
)


def to_json_safe(value: object) -> object:
    return thaw_json_value(value, field_name="provenance")


def to_json_value(value: object) -> JsonValue:
    return cast(JsonValue, to_json_safe(value))


def require_mapping(value: object, *, field_name: str) -> Mapping[str, object]:
    if isinstance(value, Mapping):
        mapping = cast(Mapping[object, object], value)
        result: dict[str, object] = {}
        for key, item in mapping.items():
            if not isinstance(key, str):
                raise PhosPyInputError(
                    f"{field_name} JSON object keys must be strings; "
                    f"got {type(key).__name__}"
                )
            if key in result:
                raise PhosPyInputError(
                    f"{field_name} contains duplicate JSON object key {key!r}"
                )
            result[key] = item
        return result
    raise PhosPyInputError(f"{field_name} must be an object")


def require_sequence(value: object, *, field_name: str) -> list[object]:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return list(cast(Sequence[object], value))
    raise PhosPyInputError(f"{field_name} must be an array")


def require_str(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a string")
    normalized = value.strip()
    if not normalized:
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return normalized


def optional_str(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    return require_str(value, field_name=field_name)


def require_raw_str(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a string")
    return value


def optional_raw_str(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    return require_raw_str(value, field_name=field_name)


def optional_mapping(value: object, *, field_name: str) -> Mapping[str, object] | None:
    if value is None:
        return None
    return require_mapping(value, field_name=field_name)


def require_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(f"{field_name} must be an int")
    return int(value)


def optional_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    return require_int(value, field_name=field_name)


def reject_legacy_provenance_fields(
    payload: Mapping[str, object],
    *,
    field_name: str,
    legacy_fields: frozenset[str],
) -> None:
    present = sorted(key for key in legacy_fields if key in payload)
    if present:
        raise_legacy_provenance_schema()


def raise_legacy_provenance_schema() -> None:
    raise PhosPyInputError(_LEGACY_PROVENANCE_SCHEMA_ERROR)


def require_shape(value: object, *, field_name: str) -> tuple[int, int]:
    sequence = require_sequence(value, field_name=field_name)
    if len(sequence) != 2:
        raise PhosPyInputError(f"{field_name} must contain exactly two integers")
    return (
        require_int(sequence[0], field_name=f"{field_name}[0]"),
        require_int(sequence[1], field_name=f"{field_name}[1]"),
    )
