from __future__ import annotations

from ..types import PredictionSvmMode, PredictionTraceFormat, PredictionTraceLevel
from ..validation.errors import PhospyValidationError
from ..validation.primitives import validate_positive_int


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
    "validate_svm_mode",
    "validate_trace_format",
    "validate_trace_level",
]
