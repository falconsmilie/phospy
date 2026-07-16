"""Contract-owned scalar validation helpers for config DTOs."""

from __future__ import annotations

from collections.abc import Collection, Iterable
from enum import Enum
from typing import TypeVar

from phospy.errors.validation import ContractValidationError

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
    error_type: ValidationErrorType = ContractValidationError,
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
    error_type: ValidationErrorType = ContractValidationError,
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
    error_type: ValidationErrorType = ContractValidationError,
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
    error_type: ValidationErrorType = ContractValidationError,
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


def require_int_at_least(
    value: object,
    *,
    field_name: str,
    minimum: int,
    error_type: ValidationErrorType = ContractValidationError,
) -> int:
    """Require an integer value greater than or equal to ``minimum``."""

    if isinstance(value, bool) or not isinstance(value, int):
        raise error_type(f"{field_name} must be an int")
    if value < minimum:
        raise error_type(f"{field_name} must be greater than or equal to {minimum}")
    return value


def require_real_between(
    value: object,
    *,
    field_name: str,
    minimum: float,
    maximum: float,
    error_type: ValidationErrorType = ContractValidationError,
) -> float:
    """Require a real numeric value constrained to an inclusive range."""

    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise error_type(
            f"{field_name} must be a float between {float(minimum):.1f} and "
            f"{float(maximum):.1f}"
        )
    numeric_value = float(value)
    if not minimum <= numeric_value <= maximum:
        raise error_type(
            f"{field_name} must be between {float(minimum):.1f} and "
            f"{float(maximum):.1f}"
        )
    return numeric_value


def require_optional_real_between(
    value: object | None,
    *,
    field_name: str,
    minimum: float,
    maximum: float,
    error_type: ValidationErrorType = ContractValidationError,
) -> float | None:
    """Require an optional real numeric value constrained to an inclusive range."""

    if value is None:
        return None
    return require_real_between(
        value,
        field_name=field_name,
        minimum=minimum,
        maximum=maximum,
        error_type=error_type,
    )


def require_int_between(
    value: object,
    *,
    field_name: str,
    minimum: int,
    maximum: int,
    error_type: ValidationErrorType = ContractValidationError,
) -> int:
    """Require an integer value constrained to an inclusive range."""

    validated = require_int_at_least(
        value,
        field_name=field_name,
        minimum=minimum,
        error_type=error_type,
    )
    if validated > maximum:
        raise error_type(f"{field_name} must be less than or equal to {maximum}")
    return validated


def require_optional_int_at_least(
    value: object | None,
    *,
    field_name: str,
    minimum: int,
    error_type: ValidationErrorType = ContractValidationError,
) -> int | None:
    """Require an optional integer value greater than or equal to ``minimum``."""

    if value is None:
        return None
    return require_int_at_least(
        value,
        field_name=field_name,
        minimum=minimum,
        error_type=error_type,
    )


def require_local_filesystem_path(
    value: object,
    *,
    field_name: str,
    error_type: ValidationErrorType = ContractValidationError,
    when_provided: bool = False,
) -> str:
    """Require one local filesystem path string and reject remote URLs."""

    path = require_non_empty_string(
        value,
        field_name=field_name,
        error_type=error_type,
        when_provided=when_provided,
    )
    if "://" in path.lower():
        raise error_type(
            f"{field_name} must be a local filesystem path; remote URLs are not "
            "supported"
        )
    return path


__all__ = [
    "ValidationErrorType",
    "coerce_policy_enum",
    "format_supported_values",
    "require_instance",
    "require_int_at_least",
    "require_int_between",
    "require_local_filesystem_path",
    "require_non_empty_string",
    "require_optional_int_at_least",
    "require_optional_real_between",
    "require_real_between",
    "require_supported_literal",
]
