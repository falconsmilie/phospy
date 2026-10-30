from __future__ import annotations

import ast
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import pytest

import phospy.api as ordinary_api
import phospy.api.advanced as specialist_api

ROOT = Path(__file__).resolve().parents[2]
CONTRACT_TEST_ROOT = ROOT / "tests" / "contract"

SUPPORTED_EXTERNAL_MODULES = {
    "phospy": frozenset(),
    "phospy.api": frozenset(ordinary_api.__all__),
    "phospy.api.advanced": frozenset(specialist_api.__all__),
}
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


@pytest.mark.parametrize(
    "route",
    (
        "phospy.datasets.models",
        "phospy.kinase.public",
        "phospy.errors",
        "phospy.frames",
        "phospy.policies",
        "phospy.validation",
        "phospy.core.errors",
        "phospy.core.frames",
        "phospy.core.policies",
        "phospy.core.validation",
        "phospy.core.quantitative.models",
        "phospy.api.builders",
        "phospy.advanced",
    ),
)
def test_external_consumer_policy_rejects_owner_and_retired_routes(route: str) -> None:
    assert _imports_unsupported_phospy_prefix(route)


def test_external_consumer_policy_accepts_only_governed_aggregate_exports() -> None:
    assert not _imports_unsupported_phospy_prefix("phospy")
    assert not _imports_unsupported_phospy_prefix("phospy.api")
    assert not _imports_unsupported_phospy_prefix(
        "phospy.api.AnalysisReadyDatasetBuilder"
    )
    assert not _imports_unsupported_phospy_prefix("phospy.api.advanced")
    assert not _imports_unsupported_phospy_prefix(
        "phospy.api.advanced.KinaseScoringConfig"
    )
    assert _imports_unsupported_phospy_prefix("phospy.AnalysisReadyDatasetBuilder")


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
    if route == "phospy":
        return False
    for module_name in sorted(SUPPORTED_EXTERNAL_MODULES, key=len, reverse=True):
        if route == module_name:
            return False
        prefix = f"{module_name}."
        if not route.startswith(prefix):
            continue
        relative = route.removeprefix(prefix)
        return (
            "." in relative or relative not in SUPPORTED_EXTERNAL_MODULES[module_name]
        )
    return route.startswith("phospy.")


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
