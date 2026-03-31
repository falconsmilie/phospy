from __future__ import annotations

import math

from ..types import PredictionSvmMode, PredictionTraceFormat, PredictionTraceLevel
from ..validation.errors import PhospyValidationError


def validate_positive_int(value: int, name: str) -> None:
    if value < 1:
        raise PhospyValidationError(f"{name} must be at least 1")


def validate_probability_threshold(
    value: float,
    name: str = "score_threshold",
) -> float:
    if not math.isfinite(value):
        raise PhospyValidationError(f"{name} must be a finite number")
    if not 0.0 <= value <= 1.0:
        raise PhospyValidationError(f"{name} must be between 0.0 and 1.0")
    return value


def validate_svm_mode(value: PredictionSvmMode) -> PredictionSvmMode:
    if value not in {"default", "r_parity"}:
        msg = "svm_mode must be one of: 'default', 'r_parity'"
        raise PhospyValidationError(msg)
    return value


def validate_trace_level(value: PredictionTraceLevel) -> PredictionTraceLevel:
    if value not in {"none", "summary", "full"}:
        msg = "trace_level must be one of: 'none', 'summary', 'full'"
        raise PhospyValidationError(msg)
    return value


def validate_trace_format(value: PredictionTraceFormat) -> PredictionTraceFormat:
    if value not in {"csv", "parquet"}:
        msg = "trace_sink_format must be one of: 'csv', 'parquet'"
        raise PhospyValidationError(msg)
    return value


__all__ = [
    "validate_positive_int",
    "validate_probability_threshold",
    "validate_svm_mode",
    "validate_trace_format",
    "validate_trace_level",
]
