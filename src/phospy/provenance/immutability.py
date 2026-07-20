"""Canonical immutable JSON-like provenance containers."""

from __future__ import annotations

import math
from collections.abc import Iterable, Iterator, Mapping, Sequence
from typing import TypeAlias, TypeVar, cast

from phospy.errors.input import PhosPyInputError

JsonPrimitive: TypeAlias = str | int | float | bool | None
_JsonErrorT = TypeVar("_JsonErrorT", bound=Exception)


class FrozenJsonMapping(Mapping[str, object]):
    """Recursively immutable provenance mapping backed by private tuple storage."""

    __slots__ = ("__items",)

    def __init__(
        self,
        value: Mapping[object, object] | Iterable[tuple[object, object]] = (),
        *,
        field_name: str = "frozen_json_mapping",
    ) -> None:
        raw_items = value.items() if isinstance(value, Mapping) else value
        items: list[tuple[str, FrozenJsonValue]] = []
        seen: set[str] = set()
        for key, item in raw_items:
            json_key = _require_json_object_key(key, field_name=field_name)
            if json_key in seen:
                raise PhosPyInputError(
                    f"{field_name} contains duplicate JSON object key {json_key!r}"
                )
            seen.add(json_key)
            items.append(
                (
                    json_key,
                    freeze_json_value(
                        item,
                        field_name=f"{field_name}.{_format_json_object_key(json_key)}",
                    ),
                )
            )
        object.__setattr__(self, "_FrozenJsonMapping__items", tuple(items))

    def __getitem__(self, key: str) -> object:
        for stored_key, value in self.__items:
            if stored_key == key:
                return value
        raise KeyError(key)

    def __iter__(self) -> Iterator[str]:
        return (key for key, _ in self.__items)

    def __len__(self) -> int:
        return len(self.__items)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({dict(self.items())!r})"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Mapping):
            try:
                return thaw_json_mapping(
                    self,
                    field_name="frozen_json_mapping",
                ) == thaw_json_mapping(other, field_name="comparison_mapping")
            except PhosPyInputError:
                return False
        return False

    def __setattr__(self, name: str, value: object) -> None:
        raise AttributeError("frozen provenance mapping is immutable")

    def __delattr__(self, name: str) -> None:
        raise AttributeError("frozen provenance mapping is immutable")

    def copy(self) -> dict[str, object]:
        """Return a fresh mutable JSON-safe copy."""

        return thaw_json_mapping(self, field_name="frozen_json_mapping")

    def __deepcopy__(self, memo: dict[int, object]) -> dict[str, object]:
        """Return a fresh mutable JSON-safe copy for dataclass helpers."""

        return thaw_json_mapping(self, field_name="frozen_json_mapping")


FrozenJsonSequence = tuple
FrozenJsonValue: TypeAlias = (
    JsonPrimitive | FrozenJsonMapping | tuple["FrozenJsonValue", ...]
)


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
        return tuple(
            freeze_json_value(item, field_name=f"{field_name}[{position}]")
            for position, item in enumerate(value)
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
    return FrozenJsonMapping(
        cast(Mapping[object, object], value), field_name=field_name
    )


def freeze_optional_json_mapping(
    value: object | None,
    *,
    field_name: str,
) -> FrozenJsonMapping | None:
    """Return a frozen JSON object mapping or None."""

    if value is None:
        return None
    return freeze_json_mapping(value, field_name=field_name)


def freeze_json_mapping_with_error_type(
    value: object,
    *,
    field_name: str,
    error_type: type[_JsonErrorT],
) -> FrozenJsonMapping:
    """Return a frozen JSON mapping while preserving a caller-owned boundary error."""

    try:
        return freeze_json_mapping(value, field_name=field_name)
    except PhosPyInputError as exc:
        if error_type is PhosPyInputError:
            raise
        raise error_type(str(exc)) from exc


def thaw_json_value(value: object, *, field_name: str) -> object:
    """Return a fresh mutable JSON-safe payload value from provenance."""

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
        thawed: dict[str, object] = {}
        for key, item in mapping.items():
            json_key = _require_json_object_key(key, field_name=field_name)
            if json_key in thawed:
                raise PhosPyInputError(
                    f"{field_name} contains duplicate JSON object key {json_key!r}"
                )
            thawed[json_key] = thaw_json_value(
                item,
                field_name=f"{field_name}.{_format_json_object_key(json_key)}",
            )
        return thawed
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
    """Return a fresh mutable JSON object mapping from provenance."""

    if not isinstance(value, Mapping):
        raise PhosPyInputError(f"{field_name} must be a mapping")
    return cast(dict[str, object], thaw_json_value(value, field_name=field_name))


def _require_json_object_key(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise PhosPyInputError(
            f"{field_name} JSON object keys must be strings; got {type(value).__name__}"
        )
    return value


def _format_json_object_key(key: str) -> str:
    return repr(key)


__all__ = [
    "FrozenJsonValue",
    "FrozenJsonMapping",
    "FrozenJsonSequence",
    "JsonPrimitive",
    "freeze_json_mapping",
    "freeze_json_mapping_with_error_type",
    "freeze_json_value",
    "freeze_optional_json_mapping",
    "thaw_json_mapping",
    "thaw_json_value",
]
