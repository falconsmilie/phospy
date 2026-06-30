from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

import phospy.api as public_api

ROOT = Path(__file__).resolve().parents[3]
README = ROOT / "README.md"

PRE_REDUCTION_EXPORT_COUNT = 202
MAX_CURATED_EXPORT_COUNT = 120

FORBIDDEN_AGGREGATE_EXPORTS = frozenset(
    {
        "AnalysisReadyDatasetModelBoundaryValidator",
        "BatchCorrectionAdequacyValidator",
        "DatasetBuildRequestValidator",
        "DatasetInputSourceValidator",
        "DifferentialAnalysisInterpreter",
        "DifferentialAnalysisExecutor",
        "KinaseWorkflowInterpreter",
        "KinaseWorkflowExecutor",
        "SignalomeWorkflowInterpreter",
        "SignalomeWorkflowExecutor",
        "ReferenceBundleBuildRequestValidator",
        "ReferenceBundleValidationReport",
        "DatasetProcessingState",
        "BatchCorrectionReport",
        "IMPORTER_QUALITY_STATUS_REPORTED",
    }
)


def _assert_from_import_fails(symbol_name: str) -> None:
    with pytest.raises(ImportError):
        exec(f"from phospy.api import {symbol_name}", {})


def _iter_python_code_blocks(source: str) -> tuple[str, ...]:
    return tuple(
        match.group("code").strip()
        for match in re.finditer(
            r"```python(?:[^\n]*)\n(?P<code>.*?)\n```",
            source,
            flags=re.DOTALL,
        )
    )


def _api_imports_from_markdown(source: str) -> set[str]:
    imported: set[str] = set()
    for block in _iter_python_code_blocks(source):
        try:
            tree = ast.parse(block)
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module == "phospy.api":
                imported.update(alias.name for alias in node.names)
    return imported


def test_public_api_surface_is_curated_smaller_than_previous_dump() -> None:
    assert len(public_api.__all__) < PRE_REDUCTION_EXPORT_COUNT
    assert len(public_api.__all__) <= MAX_CURATED_EXPORT_COUNT
    assert len(public_api.__all__) == len(set(public_api.__all__))


def test_stable_public_imports_work() -> None:
    namespace: dict[str, object] = {}

    for symbol_name in public_api._STABLE_PUBLIC_API:
        exec(f"from phospy.api import {symbol_name}", namespace)
        assert namespace[symbol_name] is getattr(public_api, symbol_name)


def test_readme_example_imports_use_stable_api_names() -> None:
    readme_imports = _api_imports_from_markdown(README.read_text(encoding="utf-8"))

    assert readme_imports
    assert readme_imports <= set(public_api._STABLE_PUBLIC_API)


def test_forbidden_internal_names_are_not_aggregate_exports() -> None:
    for symbol_name in FORBIDDEN_AGGREGATE_EXPORTS:
        assert symbol_name not in public_api.__all__
        assert not hasattr(public_api, symbol_name)
        _assert_from_import_fails(symbol_name)
