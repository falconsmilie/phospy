"""Shared provenance model primitives and scalar guards."""

from __future__ import annotations

import math
from collections.abc import Sequence
from numbers import Integral
from typing import TypeAlias

from phospy.errors.input import PhosPyInputError
from phospy.provenance.immutability import FrozenJsonValue

JsonPrimitive = str | int | float | bool | None

JsonValue: TypeAlias = (
    FrozenJsonValue
    | tuple["JsonValue", ...]
    | list["JsonValue"]
    | dict[str, "JsonValue"]
)


def _empty_json_mapping() -> dict[str, JsonValue]:
    return {}


def _required_provenance_text(value: object, *, field_name: str) -> str:
    text = str(value).strip()
    if not text:
        raise PhosPyInputError(f"{field_name} must be non-empty")
    return text


def _required_provenance_bool(value: object, *, field_name: str) -> bool:
    if not isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be a bool")
    return bool(value)


def _optional_provenance_text(value: object | None) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _optional_provenance_float(
    value: object | None, *, field_name: str
) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise PhosPyInputError(f"{field_name} must be a finite number")
    numeric = float(value)
    if not math.isfinite(numeric):
        raise PhosPyInputError(f"{field_name} must be a finite number")
    return numeric


def _required_non_negative_row_count(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise PhosPyInputError(f"{field_name} must be a non-negative integer")
    count = int(value)
    if count < 0:
        raise PhosPyInputError(f"{field_name} must be non-negative")
    return count


def _required_provenance_text_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(
        values,
        Sequence,
    ):
        raise PhosPyInputError(f"{field_name} must be a sequence of strings")
    return tuple(
        _required_provenance_text(value, field_name=f"{field_name}[]")
        for value in values
    )


def _provenance_string_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(
        values,
        Sequence,
    ):
        raise PhosPyInputError(f"{field_name} must be a sequence of strings")
    return tuple(str(value) for value in values)


def _required_shape(value: object, *, field_name: str) -> tuple[int, int]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise PhosPyInputError(f"{field_name} must contain exactly two integers")
    sequence = tuple(value)
    if len(sequence) != 2:
        raise PhosPyInputError(f"{field_name} must contain exactly two integers")
    return (
        _required_non_negative_row_count(sequence[0], field_name=f"{field_name}[0]"),
        _required_non_negative_row_count(sequence[1], field_name=f"{field_name}[1]"),
    )
