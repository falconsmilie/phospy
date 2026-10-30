"""Public exception taxonomy exports."""

from phospy.core.errors.base import PhosPyError
from phospy.core.errors.build import DatasetBuildError, PhosPyBuildError
from phospy.core.errors.input import (
    DatasetProcessingStateError,
    PhosPyInputError,
    UnsupportedInputFormatError,
)
from phospy.core.errors.provenance import (
    PhosPyProvenanceError,
    ProvenanceFingerprintError,
)
from phospy.core.errors.references import (
    PhosPyReferenceError,
    ReferenceCompatibilityError,
    ReferenceResolutionError,
    UnsupportedOrganismError,
)
from phospy.core.errors.transformations import (
    InvalidTransformationStateError,
    PhosPyTransformationError,
    TransformationStateEstablishmentError,
    TransformerExecutionError,
)
from phospy.core.errors.validation import (
    ContractValidationError,
    DatasetValidationError,
    PhosPyValidationError,
    ReferenceIdentifierNormalisationValidationError,
    ReferenceValidationError,
    TransformationValidationError,
    WorkflowValidationError,
)
from phospy.core.errors.workflows import (
    PhosPyWorkflowError,
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
    "TransformationStateEstablishmentError",
    "TransformationValidationError",
    "TransformerExecutionError",
    "UnsupportedInputFormatError",
    "UnsupportedOrganismError",
    "WorkflowBoundaryError",
    "WorkflowStageError",
    "WorkflowValidationError",
]
