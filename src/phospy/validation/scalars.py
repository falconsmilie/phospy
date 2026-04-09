from __future__ import annotations

import math
from numbers import Real

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


def validate_fraction(value: float, *, name: str) -> float:
    """Validate a finite fraction-like numeric value in the inclusive 0..1 range."""

    if isinstance(value, bool) or not isinstance(value, Real):
        raise PhospyValidationError(
            f"{name} must be a finite numeric value between 0 and 1"
        )

    resolved = float(value)
    if not math.isfinite(resolved) or not 0.0 <= resolved <= 1.0:
        raise PhospyValidationError(
            f"{name} must be a finite numeric value between 0 and 1"
        )
    return resolved


__all__ = [
    "validate_fraction",
    "validate_non_negative_int",
    "validate_positive_int",
]
