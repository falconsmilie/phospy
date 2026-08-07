from __future__ import annotations

import ast
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CONTRACT_TEST_ROOT = ROOT / "tests" / "contract"

UNSUPPORTED_PHOSPY_PREFIXES = (
    "phospy.science",
    "phospy.validation",
)
UNSUPPORTED_WORKFLOW_IMPLEMENTATION_MODULES = frozenset(
    {
        "contracts",
        "executor",
        "interpreter",
        "resolved_validator",
        "science",
        "validator",
        "workflow",
    }
)


@dataclass(frozen=True)
class ImportRoute:
    path: Path
    line_number: int
    route: str

    def as_message(self) -> str:
        relative = self.path.relative_to(ROOT).as_posix()
        return f"{relative}:{self.line_number}: {self.route}"


def test_contract_tests_use_only_supported_external_consumer_imports() -> None:
    violations = [
        imported.as_message()
        for path in sorted(CONTRACT_TEST_ROOT.rglob("test_*.py"))
        for imported in _import_routes(path)
        if _is_unsupported_contract_import(imported.route)
    ]

    assert violations == []


def _import_routes(path: Path) -> Iterator[ImportRoute]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                yield ImportRoute(path=path, line_number=node.lineno, route=alias.name)
            continue
        if isinstance(node, ast.ImportFrom):
            if node.module is None:
                continue
            yield ImportRoute(path=path, line_number=node.lineno, route=node.module)
            for alias in node.names:
                if alias.name == "*":
                    continue
                yield ImportRoute(
                    path=path,
                    line_number=node.lineno,
                    route=f"{node.module}.{alias.name}",
                )


def _is_unsupported_contract_import(route: str) -> bool:
    normalized = route.strip()
    return (
        _imports_project_test_helper(normalized)
        or _imports_unsupported_phospy_prefix(normalized)
        or _imports_private_phospy_module(normalized)
        or _imports_workflow_implementation_module(normalized)
    )


def _imports_project_test_helper(route: str) -> bool:
    return route == "tests" or route.startswith("tests.") or route == "conftest"


def _imports_unsupported_phospy_prefix(route: str) -> bool:
    return any(
        route == prefix or route.startswith(f"{prefix}.")
        for prefix in UNSUPPORTED_PHOSPY_PREFIXES
    )


def _imports_private_phospy_module(route: str) -> bool:
    parts = route.split(".")
    return parts[:1] == ["phospy"] and any(part.startswith("_") for part in parts[1:])


def _imports_workflow_implementation_module(route: str) -> bool:
    parts = route.split(".")
    return (
        len(parts) >= 4
        and parts[:2] == ["phospy", "workflows"]
        and any(
            part in UNSUPPORTED_WORKFLOW_IMPLEMENTATION_MODULES for part in parts[3:]
        )
    )
