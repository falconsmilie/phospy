from __future__ import annotations

import ast
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
PACKAGE_ROOT = PROJECT_ROOT / "src" / "phospy"
CENTRAL_DEPRECATION_HELPER = PACKAGE_ROOT / "_deprecations.py"


def test_source_modules_do_not_emit_raw_deprecation_warnings() -> None:
    violations: list[str] = []
    for path in sorted(PACKAGE_ROOT.rglob("*.py")):
        if path == CENTRAL_DEPRECATION_HELPER:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            if not _is_warnings_warn_call(node):
                continue
            if not _uses_deprecation_warning_category(node):
                continue
            relative = path.relative_to(PROJECT_ROOT)
            violations.append(f"{relative}:{node.lineno}: use phospy._deprecations")

    assert violations == []


def _is_warnings_warn_call(node: ast.Call) -> bool:
    if isinstance(node.func, ast.Attribute):
        return node.func.attr == "warn"
    return isinstance(node.func, ast.Name) and node.func.id == "warn"


def _uses_deprecation_warning_category(node: ast.Call) -> bool:
    positional_categories = node.args[1:]
    keyword_categories = [
        keyword.value
        for keyword in node.keywords
        if keyword.arg in {"category", "category_"}
    ]
    return any(
        _is_deprecation_warning_name(candidate)
        for candidate in (*positional_categories, *keyword_categories)
    )


def _is_deprecation_warning_name(node: ast.AST) -> bool:
    if isinstance(node, ast.Name):
        return node.id in {"DeprecationWarning", "PhosPyDeprecationWarning"}
    if isinstance(node, ast.Attribute):
        return node.attr in {"DeprecationWarning", "PhosPyDeprecationWarning"}
    return False
