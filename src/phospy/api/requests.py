"""Public compatibility wrapper for internal contract ownership."""

from phospy.contracts import requests as _request_contracts
from phospy.contracts.requests import (
    Contrast,
    ContrastMatrix,
    DatasetBuildRequest,
    DesignMatrix,
    DifferentialAnalysisRequest,
    EmpiricalBayesConfig,
    ExperimentalDesign,
    KinaseWorkflowRequest,
    SampleDesignRecord,
    SignalomeWorkflowRequest,
)

# Compatibility aliases intentionally re-exported at module scope.
# Keep these as explicit assignments (not import-only) so static "unused import"
# cleanups do not remove them.
DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING = (
    _request_contracts.DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING
)
DATASET_MULTI_SITE_POLICY_KEEP_JOINT = (
    _request_contracts.DATASET_MULTI_SITE_POLICY_KEEP_JOINT
)
DATASET_MULTI_SITE_POLICY_REJECT = _request_contracts.DATASET_MULTI_SITE_POLICY_REJECT
DATASET_MULTI_SITE_POLICY_SPLIT = _request_contracts.DATASET_MULTI_SITE_POLICY_SPLIT
DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE = (
    _request_contracts.DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE
)
DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED = (
    _request_contracts.DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED
)

__all__ = [
    "Contrast",
    "ContrastMatrix",
    "DesignMatrix",
    "DatasetBuildRequest",
    "DifferentialAnalysisRequest",
    "EmpiricalBayesConfig",
    "ExperimentalDesign",
    "KinaseWorkflowRequest",
    "SampleDesignRecord",
    "SignalomeWorkflowRequest",
]
