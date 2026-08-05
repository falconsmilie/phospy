# pyright: reportUnsupportedDunderAll=false
# ruff: noqa: F401
"""Stable public request models and request compatibility constants."""

from __future__ import annotations

from phospy._api_inventory import STABLE_REQUEST_API
from phospy.api._compat import deprecated_request_export
from phospy.contracts import requests as _request_contracts
from phospy.contracts.requests import (
    BatchCovariate,
    CategoricalCovariate,
    ContinuousCovariate,
    Contrast,
    ContrastMatrix,
    DatasetBuildRequest,
    DesignMatrix,
    DifferentialAnalysisRequest,
    EnrichmentIdentifierSetProvenance,
    EnrichmentIdentifierSetSourceType,
    EnrichmentSet,
    EnrichmentSetCollection,
    EnrichmentWorkflowRequest,
    ExperimentalDesign,
    FixedEffectCovariate,
    GeneSetCollection,
    KinaseWorkflowRequest,
    PhosphositeImporter,
    PhosphositeImportRequest,
    PtmSetCollection,
    SampleDesignRecord,
    SignalomeWorkflowRequest,
    all_pairwise_contrasts,
    contrasts_vs_control,
)

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

__all__ = STABLE_REQUEST_API


def __getattr__(name: str) -> object:
    return deprecated_request_export(name, old_module=__name__)
