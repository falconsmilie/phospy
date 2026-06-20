from __future__ import annotations

import ast
from pathlib import Path

_TESTS_ROOT = Path(__file__).resolve().parents[1]
_APPROVED_UNSAFE_HELPER = _TESTS_ROOT / "support" / "unsafe_dataset_states.py"
_PRIVATE_DATASET_STATE_ATTRS = {
    "_allow_opaque_site_values",
    "_comparisons",
    "_imputation_observation_metadata",
    "_phospho",
    "_sample_metadata",
    "_site_metadata",
    "_total",
}


def test_dataset_private_state_corruption_stays_in_unsafe_helper() -> None:
    violations: list[str] = []
    for path in sorted(_TESTS_ROOT.rglob("*.py")):
        if path == _APPROVED_UNSAFE_HELPER:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    _record_private_state_assignment(violations, path, target)
            elif isinstance(node, (ast.AnnAssign, ast.AugAssign)):
                _record_private_state_assignment(violations, path, node.target)
            elif isinstance(node, ast.Call) and _is_dataset_object_setattr(node):
                violations.append(_format_location(path, node.lineno))

    assert not violations, (
        "Direct AnalysisReadyPhosphoDataset state mutation is only allowed in "
        f"{_APPROVED_UNSAFE_HELPER.relative_to(_TESTS_ROOT)}; found "
        + ", ".join(violations)
    )


def _record_private_state_assignment(
    violations: list[str],
    path: Path,
    target: ast.AST,
) -> None:
    for child in ast.walk(target):
        if (
            isinstance(child, ast.Attribute)
            and child.attr in _PRIVATE_DATASET_STATE_ATTRS
        ):
            violations.append(_format_location(path, child.lineno))
            return


def _is_dataset_object_setattr(node: ast.Call) -> bool:
    if not _is_object_setattr(node.func):
        return False
    if len(node.args) < 1:
        return False
    first_arg = node.args[0]
    if _looks_like_dataset_expr(first_arg):
        return True
    if len(node.args) >= 2 and _is_private_dataset_state_name(node.args[1]):
        return True
    return False


def _is_object_setattr(node: ast.AST) -> bool:
    return (
        isinstance(node, ast.Attribute)
        and node.attr == "__setattr__"
        and isinstance(node.value, ast.Name)
        and node.value.id == "object"
    )


def _looks_like_dataset_expr(node: ast.AST) -> bool:
    if isinstance(node, ast.Name):
        return "dataset" in node.id.lower()
    if isinstance(node, ast.Attribute):
        return node.attr == "dataset" or _looks_like_dataset_expr(node.value)
    if isinstance(node, ast.Subscript):
        return _looks_like_dataset_expr(node.value)
    return False


def _is_private_dataset_state_name(node: ast.AST) -> bool:
    return (
        isinstance(node, ast.Constant)
        and isinstance(node.value, str)
        and node.value in _PRIVATE_DATASET_STATE_ATTRS
    )


def _format_location(path: Path, line_number: int) -> str:
    return f"{path.relative_to(_TESTS_ROOT)}:{line_number}"
