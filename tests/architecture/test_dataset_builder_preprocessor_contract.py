from __future__ import annotations

import ast
from pathlib import Path

import phospy.api as public_api

PROJECT_ROOT = Path(__file__).resolve().parents[2]
BUILDER_ROOT = PROJECT_ROOT / "src" / "phospy" / "science" / "datasets" / "builders"
EXECUTOR = BUILDER_ROOT / "executor.py"

_PREPROCESSOR_RUN_KEYWORDS = {
    "phospho",
    "site_metadata",
    "sample_metadata",
    "total",
    "plan",
    "corrected_preprocessing_output",
    "initial_quantitative_scale_kind",
    "initial_quantitative_meaning",
}
_PREPROCESSOR_PREFLIGHT_KEYWORDS = {
    "plan",
    "initial_quantitative_scale_kind",
    "initial_quantitative_meaning",
}


def _builder_python_files() -> tuple[Path, ...]:
    return tuple(sorted(BUILDER_ROOT.rglob("*.py")))


def _dotted_name(node: ast.AST) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        parent = _dotted_name(node.value)
        if parent is None:
            return node.attr
        return f"{parent}.{node.attr}"
    return None


def _node_mentions_preprocessor(node: ast.AST) -> bool:
    return "preprocessor" in ast.unparse(node)


def _keyword_names(call: ast.Call) -> set[str | None]:
    return {keyword.arg for keyword in call.keywords}


def test_builder_preprocessor_collaboration_has_no_reflective_negotiation() -> None:
    offenders: list[str] = []

    for path in _builder_python_files():
        relative = path.relative_to(BUILDER_ROOT)
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source, filename=str(path))
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Call)
                and _dotted_name(node.func) == "inspect.signature"
            ):
                offenders.append(f"{relative}:{node.lineno}: inspect.signature")
            if (
                isinstance(node, ast.Call)
                and _dotted_name(node.func) == "getattr"
                and len(node.args) >= 2
                and isinstance(node.args[1], ast.Constant)
                and node.args[1].value == "validate_quantitative_contracts"
            ):
                offenders.append(
                    f"{relative}:{node.lineno}: getattr(validate_quantitative_contracts)"
                )
            if (
                isinstance(node, ast.Call)
                and _dotted_name(node.func) == "isinstance"
                and any(_node_mentions_preprocessor(arg) for arg in node.args)
            ):
                offenders.append(f"{relative}:{node.lineno}: preprocessor isinstance")
            if (
                isinstance(node, ast.Call)
                and _dotted_name(node.func) == "type"
                and node.args
                and _node_mentions_preprocessor(node.args[0])
            ):
                offenders.append(f"{relative}:{node.lineno}: preprocessor type branch")

    executor_source = EXECUTOR.read_text(encoding="utf-8")
    for forbidden in (
        "Any",
        "preprocessor_kwargs",
        "self._preprocessor.run(**",
        "_preprocessor_accepts_quantitative_contract_seed",
        "_validate_quantitative_operation_contracts_before_preprocessing",
        "get_preprocessing_stage_metadata",
        "initial_quantitative_contract_state",
        "validate_and_transition",
    ):
        if forbidden in executor_source:
            offenders.append(f"executor.py: {forbidden}")

    assert offenders == []


def test_executor_calls_preprocessor_with_one_fixed_signature() -> None:
    tree = ast.parse(EXECUTOR.read_text(encoding="utf-8"), filename=str(EXECUTOR))
    run_calls: list[ast.Call] = []
    preflight_calls: list[ast.Call] = []

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        called = _dotted_name(node.func)
        if called == "self._preprocessor.run":
            run_calls.append(node)
        if called == "self._preprocessor.validate_quantitative_contracts":
            preflight_calls.append(node)

    assert len(run_calls) == 1
    assert len(preflight_calls) == 1
    assert _keyword_names(run_calls[0]) == _PREPROCESSOR_RUN_KEYWORDS
    assert _keyword_names(preflight_calls[0]) == _PREPROCESSOR_PREFLIGHT_KEYWORDS


def test_internal_preprocessor_contract_is_not_public_api() -> None:
    assert "DatasetPreprocessorContract" not in public_api.__all__
    assert not hasattr(public_api, "DatasetPreprocessorContract")
