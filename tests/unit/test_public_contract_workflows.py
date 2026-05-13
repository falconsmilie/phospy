from __future__ import annotations

import inspect
from typing import get_args, get_origin, get_type_hints

import pytest

import phospy
import phospy.api.requests as request_models
import phospy.api.workflows as workflow_models
from phospy import DifferentialAnalysisWorkflow, KinaseWorkflow, SignalomeWorkflow
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
    ExperimentalDesign,
    KinaseWorkflowRequest,
    SignalomeWorkflowRequest,
)
from phospy.api.results import (
    DifferentialAnalysisResult,
    KinaseWorkflowResult,
    SignalomeWorkflowResult,
)
from phospy.errors import WorkflowValidationError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset


def _public_methods(cls: type[object]) -> set[str]:
    return {
        name
        for name, value in inspect.getmembers(cls)
        if callable(value) and not name.startswith("_")
    }


def test_public_workflow_and_request_exports_match_contract() -> None:
    assert set(request_models.__all__) == {
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
    }
    assert set(workflow_models.__all__) == {
        "DifferentialAnalysisWorkflow",
        "KinaseWorkflow",
        "SignalomeWorkflow",
    }
    assert {
        "DifferentialAnalysisWorkflow",
        "KinaseWorkflow",
        "SignalomeWorkflow",
    }.issubset(set(phospy.__all__))
    assert "KinaseWorkflowRequest" not in phospy.__all__
    assert "SignalomeWorkflowRequest" not in phospy.__all__
    assert "KinaseWorkflowResult" not in phospy.__all__
    assert "SignalomeWorkflowResult" not in phospy.__all__
    assert "DifferentialAnalysisRequest" not in phospy.__all__
    assert "DifferentialAnalysisResult" not in phospy.__all__


def test_public_workflows_expose_run_only() -> None:
    assert _public_methods(KinaseWorkflow) == {"run"}
    assert _public_methods(SignalomeWorkflow) == {"run"}
    assert _public_methods(DifferentialAnalysisWorkflow) == {"run"}
    assert not hasattr(KinaseWorkflow, "execute")
    assert not hasattr(SignalomeWorkflow, "execute")
    assert not hasattr(DifferentialAnalysisWorkflow, "execute")
    assert not hasattr(KinaseWorkflow, "run_from_analysis_ready")
    assert not hasattr(SignalomeWorkflow, "run_from_analysis_ready")
    assert not hasattr(DifferentialAnalysisWorkflow, "run_from_analysis_ready")


def test_workflow_run_type_contracts_are_request_to_result() -> None:
    differential_top_level_hints = get_type_hints(DifferentialAnalysisWorkflow.run)
    differential_hints = get_type_hints(DifferentialAnalysisWorkflow.run)
    kinase_hints = get_type_hints(KinaseWorkflow.run)
    signalome_hints = get_type_hints(SignalomeWorkflow.run)
    assert differential_top_level_hints["request"] is DifferentialAnalysisRequest
    assert differential_top_level_hints["return"] is DifferentialAnalysisResult
    assert differential_hints["request"] is DifferentialAnalysisRequest
    assert differential_hints["return"] is DifferentialAnalysisResult
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
