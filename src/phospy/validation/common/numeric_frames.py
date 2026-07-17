"""Compatibility imports for neutral numeric DataFrame primitives."""

from __future__ import annotations

from phospy.frames.numeric import (
    ValidationErrorType,
    require_numeric_matrix,
    require_numeric_unit_interval,
)

__all__ = [
    "ValidationErrorType",
    "require_numeric_matrix",
    "require_numeric_unit_interval",
]
