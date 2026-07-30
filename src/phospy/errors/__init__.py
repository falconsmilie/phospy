"""Public exception taxonomy exports."""

from phospy.errors.base import PhosPyError
from phospy.errors.build import DatasetBuildError, PhosPyBuildError
from phospy.errors.input import (
    DatasetProcessingStateError,
    PhosPyInputError,
    UnsupportedInputFormatError,
)
from phospy.errors.provenance import (
    PhosPyProvenanceError,
    ProvenanceFingerprintError,
)
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
    ContractValidationError,
    DatasetValidationError,
    PhosPyValidationError,
    ReferenceIdentifierNormalisationValidationError,
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
    "DatasetProcessingStateError",
    "DatasetValidationError",
    "ContractValidationError",
    "InvalidTransformationStateError",
    "PhosPyBuildError",
    "PhosPyError",
    "PhosPyInputError",
    "PhosPyProvenanceError",
    "PhosPyReferenceError",
    "PhosPyTransformationError",
    "PhosPyValidationError",
    "PhosPyWorkflowError",
    "ProvenanceFingerprintError",
    "ReferenceCompatibilityError",
    "ReferenceIdentifierNormalisationValidationError",
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
