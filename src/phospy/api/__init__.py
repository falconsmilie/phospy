"""Canonical namespace where the PhosPy public API is defined.

`phospy.api` owns the package-authored public surface (requests, configs,
results, workflows, dataset and enum entry points).

Normal user code should import from top-level `phospy`, which provides the
primary stable facade. Top-level `phospy` re-exports this namespace
intentionally.
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
]
