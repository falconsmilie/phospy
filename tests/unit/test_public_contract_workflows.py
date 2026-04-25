from __future__ import annotations

import inspect
from typing import get_type_hints

import pytest

import phospy
import phospy.api.requests as request_models
import phospy.api.workflows as workflow_models
from phospy import KinaseWorkflow, SignalomeWorkflow
from phospy.api.configs import KinaseScoringConfig, SignalomeConfig
from phospy.api.requests import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    SignalomeWorkflowRequest,
)
from phospy.api.results import KinaseWorkflowResult, SignalomeWorkflowResult
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors import WorkflowValidationError


def _public_methods(cls: type[object]) -> set[str]:
    return {
        name
        for name, value in inspect.getmembers(cls)
        if callable(value) and not name.startswith("_")
    }


def test_public_workflow_and_request_exports_match_contract() -> None:
    assert set(request_models.__all__) == {
        "DatasetBuildRequest",
        "KinaseWorkflowRequest",
        "SignalomeWorkflowRequest",
    }
    assert set(workflow_models.__all__) == {"KinaseWorkflow", "SignalomeWorkflow"}
    assert {"KinaseWorkflow", "SignalomeWorkflow"}.issubset(set(phospy.__all__))
    assert "KinaseWorkflowRequest" not in phospy.__all__
    assert "SignalomeWorkflowRequest" not in phospy.__all__
    assert "KinaseWorkflowResult" not in phospy.__all__
    assert "SignalomeWorkflowResult" not in phospy.__all__


def test_public_workflows_expose_run_only() -> None:
    assert _public_methods(KinaseWorkflow) == {"run"}
    assert _public_methods(SignalomeWorkflow) == {"run"}
    assert not hasattr(KinaseWorkflow, "execute")
    assert not hasattr(SignalomeWorkflow, "execute")
    assert not hasattr(KinaseWorkflow, "run_from_analysis_ready")
    assert not hasattr(SignalomeWorkflow, "run_from_analysis_ready")


def test_workflow_run_type_contracts_are_request_to_result() -> None:
    kinase_hints = get_type_hints(KinaseWorkflow.run)
    signalome_hints = get_type_hints(SignalomeWorkflow.run)
    assert kinase_hints["request"] is KinaseWorkflowRequest
    assert kinase_hints["return"] is KinaseWorkflowResult
    assert signalome_hints["request"] is SignalomeWorkflowRequest
    assert signalome_hints["return"] is SignalomeWorkflowResult


def test_workflow_requests_keep_ingestion_outside_workflows() -> None:
    kinase_request_hints = get_type_hints(KinaseWorkflowRequest)
    signalome_request_hints = get_type_hints(
        SignalomeWorkflowRequest,
        globalns={
            **request_models.__dict__,
            "KinaseWorkflowResult": KinaseWorkflowResult,
        },
    )
    assert kinase_request_hints["dataset"] is AnalysisReadyPhosphoDataset
    assert signalome_request_hints["kinase_result"] is KinaseWorkflowResult
    assert kinase_request_hints["dataset"] is not DatasetBuildRequest


def test_workflow_configs_self_validate_local_policy_constraints() -> None:
    with pytest.raises(WorkflowValidationError, match="scoring_config.min_substrates"):
        KinaseScoringConfig(min_substrates=1)
    with pytest.raises(
        WorkflowValidationError,
        match="signalome workflow request config.network_policy",
    ):
        SignalomeConfig(network_policy="invalid")  # type: ignore[arg-type]
