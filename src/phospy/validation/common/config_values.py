"""Shared scalar/config-value validation helpers."""

from __future__ import annotations

from collections.abc import Collection, Iterable
from enum import Enum
from typing import TypeVar

ValidationErrorType = type[Exception]
_EnumT = TypeVar("_EnumT", bound=Enum)
_ValueT = TypeVar("_ValueT")


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


def require_supported_literal(
    value: object,
    *,
    field_name: str,
    supported_values: Collection[str],
    error_type: ValidationErrorType,
    sort_supported_values: bool = True,
) -> str:
    """Require one string-like literal to be in the supported values set."""

    if value in supported_values:
        return str(value)
    supported = format_supported_values(
        supported_values,
        sort_values=sort_supported_values,
    )
    raise error_type(f"{field_name} must be one of: {supported}")


def require_instance(
    value: object,
    *,
    expected_type: type[_ValueT] | tuple[type[_ValueT], ...],
    field_name: str,
    error_type: ValidationErrorType,
) -> _ValueT:
    """Require one object to be an instance of the expected type."""

    if isinstance(value, expected_type):
        return value
    raise error_type(f"{field_name} must be a {_expected_type_label(expected_type)}")


def _expected_type_label(
    expected_type: type[object] | tuple[type[object], ...],
) -> str:
    if isinstance(expected_type, tuple):
        return " or ".join(item.__name__ for item in expected_type)
    return expected_type.__name__


def require_non_empty_string(
    value: object,
    *,
    field_name: str,
    error_type: ValidationErrorType,
    when_provided: bool = False,
) -> str:
    """Require one non-empty string value."""

    if isinstance(value, str) and value.strip():
        return value
    if when_provided:
        raise error_type(f"{field_name} must be a non-empty string when provided")
    raise error_type(f"{field_name} must be a non-empty string")


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
