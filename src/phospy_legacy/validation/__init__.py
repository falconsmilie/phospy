"""Validation domain.

This package owns validation concerns organised by validation type, including
request validation, schema validation, and compatibility checks. Shared
validation error classes are imported from ``phospy.errors``.

The supported import surface here is intentionally narrower than the owning
submodules. Use the request models and primary validator entry points below.
Trusted input bundles and low-level helper utilities stay in their owning
modules.
"""

from ..errors import (
    InputCompatibilityError,
    NoCandidateKinasesError,
    PhospyError,
    PhospyValidationError,
    PredictionConfigurationError,
    RequestValidationError,
    TableSchemaError,
    TraceError,
)
from .compatibility import (
    DEFAULT_MIN_PRED_MAT_OVERLAP,
    DEFAULT_MIN_PRED_MAT_OVERLAP_FRACTION,
    ProteinCorrectionMatchSummary,
    validate_core_column_alignment,
    validate_pred_mat_overlap,
    validate_protein_correction_inputs,
    validate_signalome_alignment,
    validate_workflow_matrix_inputs,
)
from .schema.tables import (
    PhosphoInputSchema,
    PredictionScoreMatrixSchema,
    PredMatSchema,
    SiteMatrixSchema,
    SiteMatrixSourceSchema,
    TotalInputSchema,
)
from .values.enums import validate_svm_mode, validate_trace_format, validate_trace_level
from .values.numeric import (
    validate_fraction,
    validate_non_negative_int,
    validate_positive_int,
)

__all__ = [
    "InputCompatibilityError",
    "NoCandidateKinasesError",
    "PhosphoInputSchema",
    "PhospyError",
    "PhospyValidationError",
    "PredMatSchema",
    "DEFAULT_MIN_PRED_MAT_OVERLAP",
    "DEFAULT_MIN_PRED_MAT_OVERLAP_FRACTION",
    "PredictionConfigurationError",
    "PredictionScoreMatrixSchema",
    "ProteinCorrectionMatchSummary",
    "RequestValidationError",
    "SiteMatrixSchema",
    "SiteMatrixSourceSchema",
    "TableSchemaError",
    "TotalInputSchema",
    "TraceError",
    "validate_core_column_alignment",
    "validate_fraction",
    "validate_non_negative_int",
    "validate_positive_int",
    "validate_pred_mat_overlap",
    "validate_protein_correction_inputs",
    "validate_signalome_alignment",
    "validate_svm_mode",
    "validate_trace_format",
    "validate_trace_level",
    "validate_workflow_matrix_inputs",
]
