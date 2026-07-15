"""Canonical immutable JSON-like provenance containers."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import NoReturn, SupportsIndex, TypeAlias, cast

from phospy.errors.input import PhosPyInputError

JsonPrimitive: TypeAlias = str | int | float | bool | None


class FrozenJsonMapping(dict[str, object]):
    """Dict-compatible recursively frozen provenance mapping."""

    __slots__ = ()

    def __readonly(self) -> NoReturn:
        raise TypeError("frozen provenance mapping is immutable")

    def __setitem__(self, key: str, value: object) -> NoReturn:
        self.__readonly()

    def __delitem__(self, key: str) -> NoReturn:
        self.__readonly()

    def clear(self) -> NoReturn:
        self.__readonly()

    def pop(self, key: str, default: object = None) -> NoReturn:
        self.__readonly()

    def popitem(self) -> NoReturn:
        self.__readonly()

    def setdefault(self, key: str, default: object = None) -> NoReturn:
        self.__readonly()

    def update(self, *args: object, **kwargs: object) -> NoReturn:
        self.__readonly()

    def __ior__(self, other: object) -> FrozenJsonMapping:
        self.__readonly()
        return self


class FrozenJsonSequence(list[object]):
    """List-compatible recursively frozen provenance sequence."""

    __slots__ = ()

    def __readonly(self) -> NoReturn:
        raise TypeError("frozen provenance sequence is immutable")

    def __setitem__(self, index: object, value: object) -> NoReturn:
        self.__readonly()

    def __delitem__(self, index: object) -> NoReturn:
        self.__readonly()

    def append(self, value: object) -> NoReturn:
        self.__readonly()

    def clear(self) -> NoReturn:
        self.__readonly()

    def extend(self, values: object) -> NoReturn:
        self.__readonly()

    def insert(self, index: SupportsIndex, value: object) -> NoReturn:
        self.__readonly()

    def pop(self, index: SupportsIndex = -1) -> NoReturn:
        self.__readonly()

    def remove(self, value: object) -> NoReturn:
        self.__readonly()

    def reverse(self) -> NoReturn:
        self.__readonly()

    def sort(
        self,
        *,
        key: object = None,
        reverse: bool = False,
    ) -> NoReturn:
        self.__readonly()

    def __iadd__(self, other: object) -> FrozenJsonSequence:
        self.__readonly()
        return self

    def __imul__(self, other: object) -> FrozenJsonSequence:
        self.__readonly()
        return self

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Sequence) and not isinstance(
            other,
            (str, bytes, bytearray),
        ):
            return list(self) == list(other)
        return super().__eq__(other)


FrozenJsonValue: TypeAlias = JsonPrimitive | FrozenJsonMapping | FrozenJsonSequence


def freeze_json_value(value: object, *, field_name: str) -> FrozenJsonValue:
    """Return a recursively immutable JSON-like value."""

    if value is None:
        return None
    if isinstance(value, bool):
        return bool(value)
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise PhosPyInputError(f"{field_name} must be a finite JSON number")
        return float(value)
    if isinstance(value, Mapping):
        return freeze_json_mapping(value, field_name=field_name)
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return FrozenJsonSequence(
            [
                freeze_json_value(item, field_name=f"{field_name}[{position}]")
                for position, item in enumerate(value)
            ]
        )
    raise PhosPyInputError(
        f"{field_name} must contain only JSON-compatible scalars, mappings, "
        f"or sequences; got {type(value).__name__}"
    )


def freeze_json_mapping(
    value: object,
    *,
    field_name: str,
) -> FrozenJsonMapping:
    """Return a recursively immutable JSON object mapping."""

    if not isinstance(value, Mapping):
        raise PhosPyInputError(f"{field_name} must be a mapping")
    frozen = {
        str(key): freeze_json_value(item, field_name=f"{field_name}.{str(key)}")
        for key, item in value.items()
    }
    return FrozenJsonMapping(frozen)


def freeze_optional_json_mapping(
    value: object | None,
    *,
    field_name: str,
) -> FrozenJsonMapping | None:
    """Return a frozen JSON object mapping or None."""

    if value is None:
        return None
    return freeze_json_mapping(value, field_name=field_name)


def thaw_json_value(value: object, *, field_name: str) -> object:
    """Return a fresh mutable JSON-safe payload value from frozen provenance."""

    if value is None:
        return None
    if isinstance(value, bool):
        return bool(value)
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise PhosPyInputError(f"{field_name} must be a finite JSON number")
        return float(value)
    if isinstance(value, Mapping):
        mapping = cast(Mapping[object, object], value)
        return {
            str(key): thaw_json_value(item, field_name=f"{field_name}.{str(key)}")
            for key, item in mapping.items()
        }
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        sequence = cast(Sequence[object], value)
        return [
            thaw_json_value(item, field_name=f"{field_name}[{position}]")
            for position, item in enumerate(sequence)
        ]
    raise PhosPyInputError(
        f"{field_name} must contain only JSON-compatible scalars, mappings, "
        f"or sequences; got {type(value).__name__}"
    )


def thaw_json_mapping(value: object, *, field_name: str) -> dict[str, object]:
    """Return a fresh mutable JSON object mapping from frozen provenance."""

    if not isinstance(value, Mapping):
        raise PhosPyInputError(f"{field_name} must be a mapping")
    return cast(dict[str, object], thaw_json_value(value, field_name=field_name))


__all__ = [
    "FrozenJsonValue",
    "FrozenJsonMapping",
    "FrozenJsonSequence",
    "JsonPrimitive",
    "freeze_json_mapping",
    "freeze_json_value",
    "freeze_optional_json_mapping",
    "thaw_json_mapping",
    "thaw_json_value",
]
