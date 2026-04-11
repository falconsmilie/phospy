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
)

__all__ = [
    "InputCompatibilityError",
    "NoCandidateKinasesError",
    "PhospyError",
    "PhospyValidationError",
    "PredictionConfigurationError",
    "RequestValidationError",
    "TableSchemaError",
    "TraceError",
]
