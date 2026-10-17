from __future__ import annotations

import inspect
from typing import get_args, get_origin, get_type_hints

import pytest

import phospy
import phospy.api as request_models
import phospy.api as workflow_models
from phospy import (
    KinaseWorkflow,
    SignalomeWorkflow,
)
from phospy.advanced import (
    DifferentialAnalysisConfig,
    KinaseScoringConfig,
    MultipleTestingConfig,
    SignalomeConfig,
    SignalomeOutputConfig,
)
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DifferentialAnalysisRequest,
    DifferentialAnalysisResult,
    DifferentialAnalysisWorkflow,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    EnrichmentWorkflowResult,
    ExperimentalDesign,
    KinaseWorkflowRequest,
    KinaseWorkflowResult,
    SignalomeWorkflowRequest,
    SignalomeWorkflowResult,
    all_pairwise_contrasts,
    contrasts_vs_control,
)
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors import ContractValidationError

REMOVED_REQUEST_COMPATIBILITY_CONSTANTS = {
    "DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING",
    "DATASET_MULTI_SITE_POLICY_REJECT",
    "DATASET_MULTI_SITE_POLICY_SPLIT",
    "DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE",
    "DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED",
}

EXPECTED_REQUEST_EXPORTS = {
    "BatchCovariate",
    "CategoricalCovariate",
    "Contrast",
    "ContinuousCovariate",
    "DatasetBuildRequest",
    "DifferentialAnalysisRequest",
    "EnrichmentDerivedQuantitativeSetProvenance",
    "EnrichmentDerivedSetMissingValueRule",
    "EnrichmentDerivedSetSourceResultKind",
    "EnrichmentDerivedSetThresholdDirection",
    "EnrichmentDerivedSetValueMeaning",
    "EnrichmentDerivedSetValueScale",
    "EnrichmentIdentifierSetProvenance",
    "EnrichmentIdentifierSetSourceType",
    "EnrichmentSet",
    "EnrichmentSetCollection",
    "EnrichmentWorkflowRequest",
    "InputIntensityScaleEvidence",
    "ExperimentalDesign",
    "FixedEffectCovariate",
    "GeneSetCollection",
    "KinaseWorkflowRequest",
    "PhosphositeImporter",
    "PhosphositeImportRequest",
    "PtmSetCollection",
    "SampleDesignRecord",
    "SignalomeWorkflowRequest",
    "all_pairwise_contrasts",
    "contrasts_vs_control",
}


def _public_methods(cls: type[object]) -> set[str]:
    return {
        name
        for name, value in inspect.getmembers(cls)
        if callable(value) and not name.startswith("_")
    }


def test_public_workflow_and_request_exports_match_contract() -> None:
    assert EXPECTED_REQUEST_EXPORTS <= set(request_models.__all__)
    assert {
        "DifferentialAnalysisWorkflow",
        "EnrichmentWorkflow",
        "KinaseWorkflow",
        "SignalomeWorkflow",
    } <= set(workflow_models.__all__)
    assert {
        "KinaseWorkflow",
        "SignalomeWorkflow",
    }.issubset(set(phospy.__all__))
    assert "DifferentialAnalysisWorkflow" not in phospy.__all__
    assert "EnrichmentWorkflow" not in phospy.__all__
    assert "KinaseWorkflowRequest" not in phospy.__all__
    assert "SignalomeWorkflowRequest" not in phospy.__all__
    assert "KinaseWorkflowResult" not in phospy.__all__
    assert "SignalomeWorkflowResult" not in phospy.__all__
    assert "DifferentialAnalysisRequest" not in phospy.__all__
    assert "DifferentialAnalysisResult" not in phospy.__all__
    assert callable(all_pairwise_contrasts)
    assert callable(contrasts_vs_control)


def test_removed_request_compatibility_constants_are_unavailable() -> None:
    assert REMOVED_REQUEST_COMPATIBILITY_CONSTANTS.isdisjoint(request_models.__all__)
    for name in REMOVED_REQUEST_COMPATIBILITY_CONSTANTS:
        assert not hasattr(request_models, name)
    assert "DATASET_MULTI_SITE_POLICY_KEEP_JOINT" not in request_models.__all__
    assert not hasattr(request_models, "DATASET_MULTI_SITE_POLICY_KEEP_JOINT")


def test_request_star_import_exposes_public_contract_without_internals() -> None:
    namespace: dict[str, object] = {}

    exec("from phospy.api import *", namespace)

    exported_names = {name for name in namespace if name != "__builtins__"}
    assert exported_names == set(request_models.__all__)
    for name in EXPECTED_REQUEST_EXPORTS:
        assert namespace[name] is getattr(request_models, name)
    assert "_request_contracts" not in namespace
    assert "_SIGNALOME_CONFIG_TYPE_HINT_ALIAS" not in namespace
    assert "SignalomeConfig" not in namespace


def test_public_workflows_expose_run_only() -> None:
    assert _public_methods(KinaseWorkflow) == {"run"}
    assert _public_methods(SignalomeWorkflow) == {"run"}
    assert _public_methods(DifferentialAnalysisWorkflow) == {"run"}
    assert _public_methods(EnrichmentWorkflow) == {"run"}
    assert not hasattr(KinaseWorkflow, "execute")
    assert not hasattr(SignalomeWorkflow, "execute")
    assert not hasattr(DifferentialAnalysisWorkflow, "execute")
    assert not hasattr(EnrichmentWorkflow, "execute")
    assert not hasattr(KinaseWorkflow, "run_from_analysis_ready")
    assert not hasattr(SignalomeWorkflow, "run_from_analysis_ready")
    assert not hasattr(DifferentialAnalysisWorkflow, "run_from_analysis_ready")
    assert not hasattr(EnrichmentWorkflow, "run_from_analysis_ready")


def test_workflow_run_type_contracts_are_request_to_result() -> None:
    differential_top_level_hints = get_type_hints(DifferentialAnalysisWorkflow.run)
    differential_hints = get_type_hints(DifferentialAnalysisWorkflow.run)
    enrichment_hints = get_type_hints(EnrichmentWorkflow.run)
    kinase_hints = get_type_hints(KinaseWorkflow.run)
    signalome_hints = get_type_hints(SignalomeWorkflow.run)
    assert differential_top_level_hints["request"] is DifferentialAnalysisRequest
    assert differential_top_level_hints["return"] is DifferentialAnalysisResult
    assert differential_hints["request"] is DifferentialAnalysisRequest
    assert differential_hints["return"] is DifferentialAnalysisResult
    assert enrichment_hints["request"] is EnrichmentWorkflowRequest
    assert enrichment_hints["return"] is EnrichmentWorkflowResult
    assert kinase_hints["request"] is KinaseWorkflowRequest
    assert kinase_hints["return"] is KinaseWorkflowResult
    assert signalome_hints["request"] is SignalomeWorkflowRequest
    assert signalome_hints["return"] is SignalomeWorkflowResult


def test_workflow_requests_keep_ingestion_outside_workflows() -> None:
    differential_request_hints = get_type_hints(DifferentialAnalysisRequest)
    kinase_request_hints = get_type_hints(KinaseWorkflowRequest)
    signalome_request_hints = get_type_hints(SignalomeWorkflowRequest)
    assert differential_request_hints["design"] is ExperimentalDesign
    contrasts_hint = differential_request_hints["contrasts"]
    assert get_origin(contrasts_hint) is tuple
    assert get_args(contrasts_hint) == (Contrast, Ellipsis)
    assert differential_request_hints["config"] is DifferentialAnalysisConfig
    assert kinase_request_hints["dataset"] is AnalysisReadyPhosphoDataset
    assert signalome_request_hints["kinase_result"] is KinaseWorkflowResult
    assert signalome_request_hints["config"] is SignalomeConfig
    assert kinase_request_hints["dataset"] is not DatasetBuildRequest
    assert MultipleTestingConfig().method == "benjamini_hochberg"


def test_workflow_configs_self_validate_local_policy_constraints() -> None:
    with pytest.raises(ContractValidationError, match="scoring_config.min_substrates"):
        KinaseScoringConfig(reliability_profile="custom", min_substrates=1)
    with pytest.raises(
        ContractValidationError,
        match="signalome workflow request config.output.network_policy",
    ):
        SignalomeConfig(
            output=SignalomeOutputConfig(network_policy="invalid")  # type: ignore[arg-type]
        )
