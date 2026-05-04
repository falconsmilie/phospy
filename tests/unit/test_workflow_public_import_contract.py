from __future__ import annotations

import importlib
from pathlib import Path

import pytest

import phospy.workflows as workflows_package
import phospy.workflows.kinase as kinase_package
import phospy.workflows.signalome as signalome_package
from phospy.api.workflows import KinaseWorkflow as ApiKinaseWorkflow
from phospy.api.workflows import SignalomeWorkflow as ApiSignalomeWorkflow

ROOT = Path(__file__).resolve().parents[2]

REMOVED_COMPATIBILITY_IMPORTS = (
    "phospy.workflows.kinase.components",
    "phospy.workflows.signalome.components",
)


def test_supported_public_kinase_workflow_imports() -> None:
    from phospy.workflows import KinaseWorkflow
    from phospy.workflows.kinase import KinaseWorkflow as PackageKinaseWorkflow

    assert KinaseWorkflow is ApiKinaseWorkflow
    assert PackageKinaseWorkflow is ApiKinaseWorkflow


def test_supported_public_signalome_workflow_imports() -> None:
    from phospy.workflows import SignalomeWorkflow
    from phospy.workflows.signalome import SignalomeWorkflow as PackageSignalomeWorkflow

    assert SignalomeWorkflow is ApiSignalomeWorkflow
    assert PackageSignalomeWorkflow is ApiSignalomeWorkflow


@pytest.mark.parametrize("import_path", REMOVED_COMPATIBILITY_IMPORTS)
def test_removed_workflow_compatibility_imports_fail(import_path: str) -> None:
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module(import_path)


def test_workflow_package_exports_only_intended_public_symbols() -> None:
    assert set(workflows_package.__all__) == {"KinaseWorkflow", "SignalomeWorkflow"}
    assert set(kinase_package.__all__) == {"KinaseWorkflow"}
    assert set(signalome_package.__all__) == {"SignalomeWorkflow"}


def test_docs_and_examples_do_not_reference_removed_workflow_compatibility_paths() -> (
    None
):
    targets = [
        ROOT / "README.md",
        *sorted((ROOT / "docs").rglob("*.md")),
        *sorted((ROOT / "examples").rglob("*.py")),
    ]
    combined_text = "\n".join(path.read_text(encoding="utf-8") for path in targets)
    for import_path in REMOVED_COMPATIBILITY_IMPORTS:
        assert import_path not in combined_text
