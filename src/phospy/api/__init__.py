"""Authoritative namespace for the supported PhosPy public contract.

`phospy.api` owns the supported package contract: requests, configs, results,
workflows, dataset/reference entrypoints, and public exception types.
"""

from phospy.api.builders import AnalysisReadyDatasetBuilder
from phospy.api.configs import (
    DatasetComparisonBuildingConfig,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DatasetTotalProteinCorrectionConfig,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    SignalomeConfig,
)
from phospy.api.datasets import AnalysisReadyPhosphoDataset
from phospy.api.enums import Organism, ReferencePreset
from phospy.api.requests import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    SignalomeWorkflowRequest,
)
from phospy.api.results import (
    KinaseActivityResult,
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
    SignalomeWorkflowResult,
)
from phospy.api.workflows import KinaseWorkflow, SignalomeWorkflow
from phospy.errors import (
    DatasetBuildError,
    DatasetValidationError,
    InvalidTransformationStateError,
    PhosPyBuildError,
    PhosPyError,
    PhosPyInputError,
    PhosPyReferenceError,
    PhosPyTransformationError,
    PhosPyValidationError,
    PhosPyWorkflowError,
    ReferenceCompatibilityError,
    ReferenceResolutionError,
    ReferenceValidationError,
    TransformationStateEstablishmentError,
    TransformationValidationError,
    TransformerExecutionError,
    UnsupportedInputFormatError,
    UnsupportedOrganismError,
    WorkflowBoundaryError,
    WorkflowStageError,
    WorkflowValidationError,
)
from phospy.references.models import ReferenceBundle

__all__ = [
    "AnalysisReadyDatasetBuilder",
    "AnalysisReadyPhosphoDataset",
    "DatasetComparisonBuildingConfig",
    "DatasetBuildRequest",
    "DatasetMissingDataConfig",
    "DatasetPreprocessingConfig",
    "DatasetSiteMatrixConfig",
    "DatasetTotalProteinCorrectionConfig",
    "KinaseActivityConfig",
    "KinaseActivityResult",
    "KinasePredictionConfig",
    "KinasePredictionResult",
    "KinaseScoringConfig",
    "KinaseScoringResult",
    "KinaseWorkflowRequest",
    "KinaseWorkflowResult",
    "Organism",
    "ReferenceBundle",
    "ReferencePreset",
    "SignalomeConfig",
    "SignalomeWorkflow",
    "SignalomeWorkflowRequest",
    "SignalomeWorkflowResult",
    "KinaseWorkflow",
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
    "TransformationStateEstablishmentError",
    "TransformationValidationError",
    "TransformerExecutionError",
    "UnsupportedInputFormatError",
    "UnsupportedOrganismError",
    "WorkflowBoundaryError",
    "WorkflowStageError",
    "WorkflowValidationError",
]
