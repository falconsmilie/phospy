from __future__ import annotations

from dataclasses import fields, is_dataclass
from pathlib import Path
from typing import get_type_hints

from phospy.api import (
    DatasetLoadOptions,
    KinaseActivityConfig,
    PredictionRunConfig,
    SignalomeRunConfig,
    SignalomeWorkflow,
    SimpleKinaseWorkflow,
)
from phospy.api.workflow_results import SimpleKinaseWorkflowResult
from phospy.datasets import AnalysisReadyPhosphoDataset


def test_supported_public_workflows_live_in_api_modules() -> None:
    assert SimpleKinaseWorkflow.__module__ == "phospy.api.simple_workflows"
    assert SignalomeWorkflow.__module__ == "phospy.api.signalome_workflows"
    assert SimpleKinaseWorkflowResult.__module__ == "phospy.api.workflow_results"


def test_simple_workflow_result_public_shape_is_explicit() -> None:
    assert is_dataclass(SimpleKinaseWorkflowResult)
    assert [field.name for field in fields(SimpleKinaseWorkflowResult)] == [
        "analysis_ready_dataset",
        "reference_bundle",
        "scoring_result",
        "prediction_result",
        "kinase_activity_result",
    ]
    assert "pred_mat_result" not in SimpleKinaseWorkflowResult.__dataclass_fields__


def test_api_package_exports_only_supported_surface() -> None:
    import phospy.api as api

    assert set(api.__all__) == {
        "DatasetLoadOptions",
        "KinaseActivityConfig",
        "PredictionRunConfig",
        "SimpleKinaseWorkflowConfigSnapshot",
        "SignalomeRunConfig",
        "SignalomeWorkflow",
        "SimpleKinaseWorkflow",
    }
    assert hasattr(api, "SimpleKinaseWorkflow")
    assert hasattr(api, "SignalomeWorkflow")
    assert not hasattr(api, "KinaseWorkflow")
    assert not hasattr(api, "PredMatWorkflow")

    assert DatasetLoadOptions.__module__ == "phospy.api.contracts"
    assert PredictionRunConfig.__module__ == "phospy.api.contracts"
    assert KinaseActivityConfig.__module__ == "phospy.api.contracts"
    assert SignalomeRunConfig.__module__ == "phospy.api.contracts"


def test_removed_public_workflow_modules_are_not_present() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    assert not (repo_root / "src" / "phospy" / "api" / "kinase_workflows.py").exists()
    assert not (repo_root / "src" / "phospy" / "pipeline.py").exists()


def test_signalomes_package_does_not_export_trusted_execution_helper() -> None:
    import phospy.signalomes as signalomes

    assert not hasattr(signalomes, "execute_signalome_inputs")


def test_signalome_workflow_analysis_ready_signature_is_explicit() -> None:
    hints = get_type_hints(SignalomeWorkflow.run_from_analysis_ready)
    assert hints["dataset"] is AnalysisReadyPhosphoDataset
