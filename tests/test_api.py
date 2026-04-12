from __future__ import annotations

from pathlib import Path

from phospy import (
    KinaseWorkflow,
    PredMatWorkflow,
    SignalomeWorkflow,
    SimpleKinaseWorkflow,
)
from phospy.api.workflows import (
    KinaseWorkflowResult,
    PredMatWorkflowResult,
    SimpleKinaseWorkflowResult,
)


def test_public_workflows_are_defined_under_api_package() -> None:
    assert KinaseWorkflow.__module__ == "phospy.api.workflows"
    assert PredMatWorkflow.__module__ == "phospy.api.workflows"
    assert SignalomeWorkflow.__module__ == "phospy.api.workflows"
    assert SimpleKinaseWorkflow.__module__ == "phospy.api.workflows"
    assert KinaseWorkflowResult.__module__ == "phospy.api.workflows"
    assert PredMatWorkflowResult.__module__ == "phospy.api.workflows"
    assert SimpleKinaseWorkflowResult.__module__ == "phospy.api.workflows"


def test_legacy_workflow_module_has_been_removed() -> None:
    workflow_module = (
        Path(__file__).resolve().parents[1] / "src" / "phospy" / "workflow.py"
    )
    assert not workflow_module.exists()


def test_api_package_exports_only_supported_workflows() -> None:
    import phospy.api as api

    assert set(api.__all__) == {
        "KinaseWorkflow",
        "PredMatWorkflow",
        "SignalomeWorkflow",
        "SimpleKinaseWorkflow",
    }
    assert not hasattr(api, "KinaseWorkflowResult")
    assert not hasattr(api, "PredMatWorkflowResult")
    assert not hasattr(api, "SimpleKinaseWorkflowResult")
