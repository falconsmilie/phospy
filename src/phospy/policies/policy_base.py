from __future__ import annotations

from collections.abc import Iterable
from enum import Enum
from typing import TypeVar

from phospy.errors.input import PhosPyInputError

_PolicyEnumT = TypeVar("_PolicyEnumT", bound="PolicyEnum")
_EnumT = TypeVar("_EnumT", bound=Enum)
ValidationErrorType = type[Exception]


def format_supported_values(
    supported_values: Iterable[str],
    *,
    sort_values: bool = True,
) -> str:
    """Return one stable, comma-separated supported-values string."""

    values = tuple(str(value) for value in supported_values)
    if sort_values:
        values = tuple(sorted(values))
    return ", ".join(values)


def coerce_policy_enum(
    enum_type: type[_EnumT],
    value: object,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> _EnumT:
    """Coerce one value into an enum-backed policy with strict error messaging."""

    if isinstance(value, enum_type):
        return value
    if isinstance(value, str):
        normalized = value.strip()
        try:
            return enum_type(normalized)
        except ValueError:
            pass
    supported = format_supported_values(
        (str(member.value) for member in enum_type),
        sort_values=False,
    )
    raise error_type(f"{field_name} must be one of: {supported}; got {value!r}")


class PolicyEnum(str, Enum):
    """Base class for stable policy enums with strict parsing helpers."""

    def __str__(self) -> str:
        return self.value

    @classmethod
    def parse(
        cls: type[_PolicyEnumT],
        value: object,
        *,
        field_name: str,
    ) -> _PolicyEnumT:
        return coerce_policy_enum(
            cls,
            value,
            field_name=field_name,
            error_type=PhosPyInputError,
        )


__all__ = ["PolicyEnum", "coerce_policy_enum", "format_supported_values"]
