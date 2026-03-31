from __future__ import annotations

from .errors import PhospyValidationError


def validate_non_negative_int(value: int, name: str) -> int:
    """Validate that an integer configuration value is zero or greater."""

    if value < 0:
        raise PhospyValidationError(f"{name} must be at least 0")
    return value


def validate_positive_int(value: int, name: str) -> int:
    """Validate that an integer configuration value is one or greater."""

    if value < 1:
        raise PhospyValidationError(f"{name} must be at least 1")
    return value


__all__ = ["validate_non_negative_int", "validate_positive_int"]
