"""Public exception taxonomy exports."""

from phospy.errors.base import PhosPyError
from phospy.errors.build import DatasetBuildError, PhosPyBuildError
from phospy.errors.input import PhosPyInputError, UnsupportedInputFormatError
from phospy.errors.references import (
    PhosPyReferenceError,
    ReferenceCompatibilityError,
    ReferenceResolutionError,
    UnsupportedOrganismError,
)
from phospy.errors.transformations import (
    InvalidTransformationStateError,
    PhosPyTransformationError,
    TransformationStateEstablishmentError,
    TransformerExecutionError,
)
from phospy.errors.validation import (
    DatasetValidationError,
    PhosPyValidationError,
    ReferenceValidationError,
    TransformationValidationError,
    WorkflowValidationError,
)
from phospy.errors.workflows import (
    PhosPyWorkflowError,
    SignalomeModuleCountValidationError,
    SignalomeScaleError,
    WorkflowBoundaryError,
    WorkflowStageError,
)

__all__ = [
    "DatasetBuildError",
    "DatasetValidationError",
    "InvalidTransformationStateError",
    "PhosPyBuildError",
    "PhosPyError",
    "PhosPyInputError",
    "PhosPyReferenceError",
    "PhosPyTransformationError",
    "PhosPyValidationError",
    "PhosPyWorkflowError",
    "ReferenceCompatibilityError",
    "ReferenceResolutionError",
    "ReferenceValidationError",
    "SignalomeModuleCountValidationError",
    "SignalomeScaleError",
    "TransformationStateEstablishmentError",
    "TransformationValidationError",
    "TransformerExecutionError",
    "UnsupportedInputFormatError",
    "UnsupportedOrganismError",
    "WorkflowBoundaryError",
    "WorkflowStageError",
    "WorkflowValidationError",
]
