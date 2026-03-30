from __future__ import annotations

from ..types import DegenerateProbabilityPolicy, PredictionSvmMode
from ..validation.errors import PhospyValidationError


def validate_positive_int(value: int, name: str) -> None:
    if value < 1:
        raise PhospyValidationError(f"{name} must be at least 1")


def validate_degenerate_probability_policy(
    value: DegenerateProbabilityPolicy,
) -> DegenerateProbabilityPolicy:
    if value not in {"uniform", "error"}:
        msg = "degenerate_probability_policy must be one of: 'uniform', 'error'"
        raise PhospyValidationError(msg)
    return value


def validate_svm_mode(value: PredictionSvmMode) -> PredictionSvmMode:
    if value not in {"default", "r_parity"}:
        msg = "svm_mode must be one of: 'default', 'r_parity'"
        raise PhospyValidationError(msg)
    return value


__all__ = [
    "validate_degenerate_probability_policy",
    "validate_positive_int",
    "validate_svm_mode",
]
