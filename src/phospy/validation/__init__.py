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
    ProteinCorrectionMatchSummary,
    validate_core_column_alignment,
    validate_pred_mat_overlap,
    validate_protein_correction_inputs,
    validate_signalome_alignment,
    validate_workflow_matrix_inputs,
)
from .requests import (
    CorePipelineRequest,
    KinaseActivityRequest,
    KinaseWorkflowRequest,
    PredictionRequest,
    SignalomeRequest,
    validate_analysis_request,
    validate_pipeline_construction_request,
    validate_pipeline_runtime_compatibility,
    validate_signalome_request,
    validate_workflow_request,
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
    "CorePipelineRequest",
    "InputCompatibilityError",
    "KinaseActivityRequest",
    "KinaseWorkflowRequest",
    "NoCandidateKinasesError",
    "PhosphoInputSchema",
    "PhospyError",
    "PhospyValidationError",
    "PredMatSchema",
    "PredictionConfigurationError",
    "PredictionRequest",
    "PredictionScoreMatrixSchema",
    "ProteinCorrectionMatchSummary",
    "RequestValidationError",
    "SignalomeRequest",
    "SiteMatrixSchema",
    "SiteMatrixSourceSchema",
    "TableSchemaError",
    "TotalInputSchema",
    "TraceError",
    "validate_analysis_request",
    "validate_core_column_alignment",
    "validate_fraction",
    "validate_non_negative_int",
    "validate_pipeline_construction_request",
    "validate_pipeline_runtime_compatibility",
    "validate_positive_int",
    "validate_pred_mat_overlap",
    "validate_protein_correction_inputs",
    "validate_signalome_alignment",
    "validate_signalome_request",
    "validate_svm_mode",
    "validate_trace_format",
    "validate_trace_level",
    "validate_workflow_matrix_inputs",
    "validate_workflow_request",
]
