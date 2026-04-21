"""Primary user-facing facade for the PhosPy public API.

Import guidance:
- User code should import from top-level `phospy`.
- `phospy.api` is the canonical authored API-definition namespace.
- Top-level exports are a curated facade, not a full mirror of all public modules.
- `phospy.errors` exposes the full error taxonomy; top-level keeps a focused
  subset of user-handleable exceptions.
"""

from phospy.api import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DatasetBuildRequest,
    DatasetComparisonBuildingConfig,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DatasetTotalProteinCorrectionConfig,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    KinaseWorkflowResult,
    Organism,
    ReferenceBundle,
    ReferencePreset,
    SignalomeConfig,
    SignalomeWorkflow,
    SignalomeWorkflowRequest,
    SignalomeWorkflowResult,
)
from phospy.errors import (
    PhosPyBuildError,
    PhosPyError,
    PhosPyInputError,
    PhosPyReferenceError,
    PhosPyTransformationError,
    PhosPyValidationError,
    PhosPyWorkflowError,
    UnsupportedInputFormatError,
    UnsupportedOrganismError,
    WorkflowBoundaryError,
)

__all__ = [
    "AnalysisReadyDatasetBuilder",
    "AnalysisReadyPhosphoDataset",
    "DatasetBuildRequest",
    "DatasetComparisonBuildingConfig",
    "DatasetMissingDataConfig",
    "DatasetPreprocessingConfig",
    "DatasetSiteMatrixConfig",
    "DatasetTotalProteinCorrectionConfig",
    "KinaseActivityConfig",
    "KinasePredictionConfig",
    "KinaseScoringConfig",
    "KinaseWorkflow",
    "KinaseWorkflowRequest",
    "KinaseWorkflowResult",
    "Organism",
    "ReferenceBundle",
    "ReferencePreset",
    "SignalomeConfig",
    "SignalomeWorkflow",
    "SignalomeWorkflowRequest",
    "SignalomeWorkflowResult",
    "PhosPyBuildError",
    "PhosPyError",
    "PhosPyInputError",
    "PhosPyReferenceError",
    "PhosPyTransformationError",
    "PhosPyValidationError",
    "PhosPyWorkflowError",
    "UnsupportedInputFormatError",
    "UnsupportedOrganismError",
    "WorkflowBoundaryError",
]
