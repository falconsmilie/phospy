"""Public compatibility wrapper for internal contract ownership."""

from phospy.contracts import requests as _request_contracts
from phospy.contracts.configs import SignalomeConfig
from phospy.contracts.requests import (
    BatchCovariate,
    CategoricalCovariate,
    ContinuousCovariate,
    Contrast,
    ContrastMatrix,
    DatasetBuildRequest,
    DesignMatrix,
    DifferentialAnalysisRequest,
    EmpiricalBayesConfig,
    ExperimentalDesign,
    FixedEffectCovariate,
    KinaseWorkflowRequest,
    PhosphositeImporter,
    PhosphositeImportRequest,
    SampleDesignRecord,
    SignalomeWorkflowRequest,
)

# SignalomeConfig is an implementation dependency here; configs remain owned by
# phospy.api.configs, so this private helper is intentionally not exported.
_SIGNALOME_CONFIG_TYPE_HINT_ALIAS = SignalomeConfig

# Compatibility constants intentionally re-exported at module scope.
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
    "BatchCovariate",
    "Contrast",
    "ContrastMatrix",
    "CategoricalCovariate",
    "ContinuousCovariate",
    "DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING",
    "DATASET_MULTI_SITE_POLICY_KEEP_JOINT",
    "DATASET_MULTI_SITE_POLICY_REJECT",
    "DATASET_MULTI_SITE_POLICY_SPLIT",
    "DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE",
    "DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED",
    "DesignMatrix",
    "DatasetBuildRequest",
    "DifferentialAnalysisRequest",
    "EmpiricalBayesConfig",
    "ExperimentalDesign",
    "FixedEffectCovariate",
    "KinaseWorkflowRequest",
    "PhosphositeImporter",
    "PhosphositeImportRequest",
    "SampleDesignRecord",
    "SignalomeWorkflowRequest",
]
