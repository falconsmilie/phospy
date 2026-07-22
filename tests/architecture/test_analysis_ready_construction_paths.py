from __future__ import annotations

import ast
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = PROJECT_ROOT / "src" / "phospy"
_BUILDER_OUTPUT_ALLOWED_PATHS = {
    Path("science/datasets/builders/executor.py"),
    Path("science/datasets/models.py"),
}


def _production_python_files() -> tuple[Path, ...]:
    return tuple(sorted(SOURCE_ROOT.rglob("*.py")))


def _call_name(node: ast.Call) -> str | None:
    func = node.func
    if isinstance(func, ast.Name):
        return func.id
    if isinstance(func, ast.Attribute):
        return func.attr
    return None


def test_no_production_code_calls_analysis_ready_direct_constructor() -> None:
    offenders: list[str] = []
    for path in _production_python_files():
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            if _call_name(node) == "AnalysisReadyPhosphoDataset":
                offenders.append(f"{path.relative_to(SOURCE_ROOT)}:{node.lineno}")

    assert offenders == []


def test_builder_output_factory_is_only_used_by_dataset_builder() -> None:
    offenders: list[str] = []
    for path in _production_python_files():
        relative = path.relative_to(SOURCE_ROOT)
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            if _call_name(node) != "_from_builder_output":
                continue
            if relative not in _BUILDER_OUTPUT_ALLOWED_PATHS:
                offenders.append(f"{relative}:{node.lineno}")

    assert offenders == []
