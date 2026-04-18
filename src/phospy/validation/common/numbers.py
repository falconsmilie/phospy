"""Shared numeric validation helpers."""

from __future__ import annotations

from phospy.errors.validation import PhosPyValidationError

ValidationErrorType = type[PhosPyValidationError]


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
