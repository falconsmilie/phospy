"""Shared application error hierarchy.

This package owns the cross-domain error hierarchy used by validation, I/O, and
other application layers. Domain packages should import shared error types from
this package rather than defining parallel hierarchies.
"""

from .base import PhospyError
from .validation import (
    InputCompatibilityError,
    NoCandidateKinasesError,
    PhospyValidationError,
    PredictionConfigurationError,
    RequestValidationError,
    TableSchemaError,
    TraceError,
    format_empty_prediction_matrix_message,
    format_no_candidate_kinases_message,
    format_overlap_failure_message,
)

__all__ = [
    "format_empty_prediction_matrix_message",
    "format_no_candidate_kinases_message",
    "format_overlap_failure_message",
    "InputCompatibilityError",
    "NoCandidateKinasesError",
    "PhospyError",
    "PhospyValidationError",
    "PredictionConfigurationError",
    "RequestValidationError",
    "TableSchemaError",
    "TraceError",
]
