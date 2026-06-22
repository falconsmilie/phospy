from __future__ import annotations

import inspect
from typing import get_args, get_origin, get_type_hints

import pytest

import phospy
import phospy.api.requests as request_models
import phospy.api.workflows as workflow_models
from phospy import (
    DifferentialAnalysisWorkflow,
    KinaseWorkflow,
    SignalomeWorkflow,
)
from phospy.api.configs import (
    DifferentialAnalysisConfig,
    KinaseScoringConfig,
    MultipleTestingConfig,
    SignalomeConfig,
    SignalomeOutputConfig,
)
from phospy.api.requests import (
    Contrast,
    DatasetBuildRequest,
    DifferentialAnalysisRequest,
    EnrichmentWorkflowRequest,
    ExperimentalDesign,
    KinaseWorkflowRequest,
    SignalomeWorkflowRequest,
    all_pairwise_contrasts,
    contrasts_vs_control,
)
from phospy.api.results import (
    DifferentialAnalysisResult,
    EnrichmentWorkflowResult,
    KinaseWorkflowResult,
    SignalomeWorkflowResult,
)
from phospy.api.workflows import EnrichmentWorkflow
from phospy.errors import WorkflowValidationError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset

INTENTIONAL_REQUEST_COMPATIBILITY_CONSTANTS = {
    "DATASET_MULTI_SITE_POLICY_EXCLUDE_FROM_SEQUENCE_SCORING",
    "DATASET_MULTI_SITE_POLICY_KEEP_JOINT",
    "DATASET_MULTI_SITE_POLICY_REJECT",
    "DATASET_MULTI_SITE_POLICY_SPLIT",
    "DATASET_SITE_RESOLUTION_MODE_PEPTIDE_EVIDENCE",
    "DATASET_SITE_RESOLUTION_MODE_SITE_LEVEL_RESOLVED",
}

EXPECTED_REQUEST_EXPORTS = {
    "BatchCovariate",
    "CategoricalCovariate",
    "Contrast",
    "ContrastMatrix",
    "ContinuousCovariate",
    "DesignMatrix",
    "DatasetBuildRequest",
    "DifferentialAnalysisRequest",
    "EmpiricalBayesConfig",
    "EnrichmentIdentifierKind",
    "EnrichmentSet",
    "EnrichmentSetCollection",
    "EnrichmentWorkflowRequest",
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
} | INTENTIONAL_REQUEST_COMPATIBILITY_CONSTANTS


def _public_methods(cls: type[object]) -> set[str]:
    return {
        name
        for name, value in inspect.getmembers(cls)
        if callable(value) and not name.startswith("_")
    }


def test_public_workflow_and_request_exports_match_contract() -> None:
    assert set(request_models.__all__) == EXPECTED_REQUEST_EXPORTS
    assert set(workflow_models.__all__) == {
        "DifferentialAnalysisWorkflow",
        "EnrichmentWorkflow",
        "KinaseWorkflow",
        "SignalomeWorkflow",
    }
    assert {
        "DifferentialAnalysisWorkflow",
        "KinaseWorkflow",
        "SignalomeWorkflow",
    }.issubset(set(phospy.__all__))
    assert "EnrichmentWorkflow" not in phospy.__all__
    assert "KinaseWorkflowRequest" not in phospy.__all__
    assert "SignalomeWorkflowRequest" not in phospy.__all__
    assert "KinaseWorkflowResult" not in phospy.__all__
    assert "SignalomeWorkflowResult" not in phospy.__all__
    assert "DifferentialAnalysisRequest" not in phospy.__all__
    assert "DifferentialAnalysisResult" not in phospy.__all__
    assert callable(all_pairwise_contrasts)
    assert callable(contrasts_vs_control)


def test_request_compatibility_constants_are_public_exports() -> None:
    assert INTENTIONAL_REQUEST_COMPATIBILITY_CONSTANTS <= set(request_models.__all__)


def test_request_star_import_exposes_public_contract_without_internals() -> None:
    namespace: dict[str, object] = {}

    exec("from phospy.api.requests import *", namespace)

    exported_names = {name for name in namespace if name != "__builtins__"}
    assert exported_names == EXPECTED_REQUEST_EXPORTS
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
    signalome_request_hints = get_type_hints(
        SignalomeWorkflowRequest,
        globalns={
            **request_models.__dict__,
            "KinaseWorkflowResult": KinaseWorkflowResult,
        },
    )
    assert differential_request_hints["design"] is ExperimentalDesign
    contrasts_hint = differential_request_hints["contrasts"]
    assert get_origin(contrasts_hint) is tuple
    assert get_args(contrasts_hint) == (Contrast, Ellipsis)
    assert differential_request_hints["config"] is DifferentialAnalysisConfig
    assert kinase_request_hints["dataset"] is AnalysisReadyPhosphoDataset
    assert signalome_request_hints["kinase_result"] is KinaseWorkflowResult
    assert kinase_request_hints["dataset"] is not DatasetBuildRequest
    assert MultipleTestingConfig().method == "benjamini_hochberg"


def test_workflow_configs_self_validate_local_policy_constraints() -> None:
    with pytest.raises(WorkflowValidationError, match="scoring_config.min_substrates"):
        KinaseScoringConfig(min_substrates=1)
    with pytest.raises(
        WorkflowValidationError,
        match="signalome workflow request config.output.network_policy",
    ):
        SignalomeConfig(
            output=SignalomeOutputConfig(network_policy="invalid")  # type: ignore[arg-type]
        )
