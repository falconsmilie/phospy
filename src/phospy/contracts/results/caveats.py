"""Compatibility re-exports for common result caveats."""

from __future__ import annotations

from phospy.contracts.result_caveats import (
    ResultCaveat,
    ResultCaveatSeverity,
    validate_result_caveats,
)

__all__ = [
    "ResultCaveat",
    "ResultCaveatSeverity",
    "validate_result_caveats",
]
