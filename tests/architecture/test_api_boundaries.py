from __future__ import annotations

import ast
from pathlib import Path

import pytest

import phospy.api as public_api

PROJECT_ROOT = Path(__file__).resolve().parents[2]
API_INIT = PROJECT_ROOT / "src" / "phospy" / "api" / "__init__.py"

FORBIDDEN_VALIDATOR_NAMES = frozenset(
    {
        "AnalysisReadyDatasetModelBoundaryValidator",
        "BatchCorrectionAdequacyValidator",
        "DatasetBuildRequestValidator",
        "DatasetInputSourceValidator",
        "DatasetPreprocessingConfigValidator",
        "KinaseWorkflowValidator",
        "ReferenceBundleBuildRequestValidator",
        "SignalomeWorkflowValidator",
    }
)

FORBIDDEN_LOW_LEVEL_WORKFLOW_NAMES = frozenset(
    {
        "DifferentialAnalysisInterpreter",
        "DifferentialAnalysisExecutor",
        "EnrichmentWorkflowInterpreter",
        "EnrichmentWorkflowExecutor",
        "KinaseWorkflowInterpreter",
        "KinaseWorkflowExecutor",
        "SignalomeWorkflowInterpreter",
        "SignalomeWorkflowExecutor",
    }
)


def _assert_from_import_fails(symbol_name: str) -> None:
    with pytest.raises(ImportError):
        exec(f"from phospy.api import {symbol_name}", {})


def _imported_modules(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    modules: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            modules.extend(alias.name for alias in node.names)
        if isinstance(node, ast.ImportFrom) and node.module is not None:
            modules.append(node.module)
    return modules


def test_public_api_root_does_not_import_validation_or_low_level_workflow_modules() -> (
    None
):
    imported = _imported_modules(API_INIT)

    offenders = sorted(
        module
        for module in imported
        if module.startswith("phospy.validation")
        or module.endswith(".validator")
        or module.endswith(".interpreter")
        or module.endswith(".executor")
    )

    assert offenders == []


def test_validators_are_not_public_api_exports() -> None:
    for symbol_name in FORBIDDEN_VALIDATOR_NAMES:
        assert symbol_name not in public_api.__all__
        assert not hasattr(public_api, symbol_name)
        _assert_from_import_fails(symbol_name)


def test_low_level_workflow_internals_are_not_public_api_exports() -> None:
    for symbol_name in FORBIDDEN_LOW_LEVEL_WORKFLOW_NAMES:
        assert symbol_name not in public_api.__all__
        assert not hasattr(public_api, symbol_name)
        _assert_from_import_fails(symbol_name)


def test_future_validator_export_by_name_fails_architecture_boundary() -> None:
    leaked = sorted(name for name in public_api.__all__ if name.endswith("Validator"))

    assert leaked == []
