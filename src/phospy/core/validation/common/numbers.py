"""Shared numeric validation helpers."""

from __future__ import annotations

ValidationErrorType = type[Exception]


def require_int_at_least(
    value: object,
    *,
    field_name: str,
    minimum: int,
    error_type: ValidationErrorType,
) -> int:
    """Require an integer value that is greater than or equal to ``minimum``."""

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
    error_type: ValidationErrorType,
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
    error_type: ValidationErrorType,
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
    error_type: ValidationErrorType,
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
    error_type: ValidationErrorType,
) -> int | None:
    """Require an optional integer value that is greater than or equal to ``minimum``."""

    if value is None:
        return None
    return require_int_at_least(
        value,
        field_name=field_name,
        minimum=minimum,
        error_type=error_type,
    )
