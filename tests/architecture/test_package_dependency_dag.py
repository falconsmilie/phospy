from __future__ import annotations

import ast
import os
import subprocess
import sys
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
PACKAGE_ROOT = SRC_ROOT / "phospy"

ALLOWED_PACKAGE_EDGES = frozenset(
    {
        ("phospy._deprecations", "phospy._api_inventory"),
        ("phospy.advanced", "phospy._api_inventory"),
        ("phospy.advanced", "phospy.contracts"),
        ("phospy.advanced", "phospy.io"),
        ("phospy.advanced", "phospy.science"),
        ("phospy.api", "phospy._deprecations"),
        ("phospy.api", "phospy._api_inventory"),
        ("phospy.api", "phospy.advanced"),
        ("phospy.api", "phospy.contracts"),
        ("phospy.api", "phospy.errors"),
        ("phospy.api", "phospy.io"),
        ("phospy.api", "phospy.science"),
        ("phospy.api", "phospy.validation"),
        ("phospy.api", "phospy.workflows"),
        ("phospy.contracts", "phospy._deprecations"),
        ("phospy.contracts", "phospy.errors"),
        ("phospy.contracts", "phospy.frames"),
        ("phospy.contracts", "phospy.policies"),
        ("phospy.contracts", "phospy.provenance"),
        ("phospy.contracts", "phospy.science"),
        ("phospy.frames", "phospy.errors"),
        ("phospy.io", "phospy._deprecations"),
        ("phospy.io", "phospy.contracts"),
        ("phospy.io", "phospy.errors"),
        ("phospy.io", "phospy.provenance"),
        ("phospy.io", "phospy.science"),
        ("phospy.io", "phospy.validation"),
        ("phospy.policies", "phospy.errors"),
        ("phospy.provenance", "phospy.errors"),
        ("phospy.science", "phospy._deprecations"),
        ("phospy.science", "phospy.errors"),
        ("phospy.science", "phospy.frames"),
        ("phospy.science", "phospy.policies"),
        ("phospy.science", "phospy.provenance"),
        ("phospy.tables", "phospy.frames"),
        ("phospy.tables", "phospy.science"),
        ("phospy.validation", "phospy.contracts"),
        ("phospy.validation", "phospy.errors"),
        ("phospy.validation", "phospy.frames"),
        ("phospy.validation", "phospy.provenance"),
        ("phospy.validation", "phospy.science"),
        ("phospy.workflows", "phospy._deprecations"),
        ("phospy.workflows", "phospy.contracts"),
        ("phospy.workflows", "phospy.errors"),
        ("phospy.workflows", "phospy.provenance"),
        ("phospy.workflows", "phospy.science"),
        ("phospy.workflows", "phospy.validation"),
    }
)

FORBIDDEN_PACKAGE_EDGES = frozenset(
    {
        ("phospy.errors", "phospy.api"),
        ("phospy.errors", "phospy.contracts"),
        ("phospy.errors", "phospy.io"),
        ("phospy.errors", "phospy.science"),
        ("phospy.errors", "phospy.tables"),
        ("phospy.errors", "phospy.validation"),
        ("phospy.errors", "phospy.workflows"),
        ("phospy.io", "phospy.api"),
        ("phospy.science", "phospy.contracts"),
        ("phospy.science", "phospy.io"),
        ("phospy.science", "phospy.tables"),
        ("phospy.science", "phospy.validation"),
        ("phospy.science", "phospy.workflows"),
        ("phospy.validation", "phospy.workflows"),
    }
)

CONTRACTS_PUBLIC_SCIENCE_PREFIXES = frozenset(
    {
        "phospy.science.configs",
    }
)

CONTRACTS_SCIENCE_FACADE_ROLE_MARKER = "__phospy_contracts_facade_role__"

CONTRACTS_PUBLIC_SCIENCE_FACADE_ROLES = frozenset(
    {
        "science_owned_public_constant",
        "science_owned_public_enum",
        "science_owned_public_helper",
        "science_owned_public_model",
        "science_owned_public_table_contract",
    }
)

CONTRACTS_FORBIDDEN_FACADE_DECLARATION_SUFFIXES = (
    "Validator",
    "Loader",
    "Builder",
    "ConstructionService",
    "Stage",
    "Executor",
    "Runner",
    "Orchestrator",
    "Resolver",
)

CONTRACTS_FORBIDDEN_FACADE_DECLARATION_FRAGMENTS = (
    "WorkflowExecutor",
    "WorkflowRunner",
    "WorkflowOrchestrator",
    "InternalDatasetView",
    "InternalView",
)

CONTRACTS_FORBIDDEN_FACADE_FUNCTION_PREFIXES = (
    "build_",
    "construct_",
    "execute_",
    "load_",
    "orchestrate_",
    "resolve_",
    "run_",
    "validate_",
)

CONTRACTS_FORBIDDEN_SCIENCE_STRUCTURE_SEGMENTS = frozenset(
    {
        "builder",
        "builders",
        "construction",
        "executor",
        "executors",
        "execution",
    }
)

CONTRACTS_FORBIDDEN_SCIENCE_WORKFLOW_LEAVES = frozenset(
    {
        "assembler",
        "interpreter",
        "orchestrator",
        "runner",
        "workflow",
        "workflows",
    }
)

PRIVATE_SCIENCE_MODULES = frozenset(
    {
        "phospy.science.differential.internal_view",
        "phospy.science.datasets.internal_view",
        "phospy.science.prediction.internal_view",
    }
)


@dataclass(frozen=True, slots=True)
class ImportRecord:
    source_module: str
    source_path: Path
    target: str
    line: int


@dataclass(frozen=True, slots=True)
class ImportGraph:
    modules: frozenset[str]
    module_paths: Mapping[str, Path]
    module_edges: dict[str, frozenset[str]]
    package_edges: frozenset[tuple[str, str]]
    records: tuple[ImportRecord, ...]


@dataclass(frozen=True, slots=True)
class ContractsScienceFacadeRoleMarker:
    module_name: str
    source_path: Path
    role: str | None
    line: int


def test_package_dependency_graph_has_no_static_cycles() -> None:
    graph = _build_import_graph()
    cycles = [
        tuple(sorted(component))
        for component in _strongly_connected_components(
            _package_adjacency(graph.package_edges)
        )
        if len(component) > 1
    ]

    assert cycles == [], _format_cycles(cycles)


def test_package_dependency_graph_matches_allowed_edge_table() -> None:
    graph = _build_import_graph()

    unexpected = sorted(graph.package_edges - ALLOWED_PACKAGE_EDGES)
    stale_allowed = sorted(ALLOWED_PACKAGE_EDGES - graph.package_edges)

    assert unexpected == [], _format_package_edges(unexpected, graph.records)
    assert stale_allowed == [], "\n".join(f"{s} -> {t}" for s, t in stale_allowed)


def test_forbidden_package_edges_are_absent() -> None:
    graph = _build_import_graph()
    offenders = sorted(graph.package_edges & FORBIDDEN_PACKAGE_EDGES)

    assert offenders == [], _format_package_edges(offenders, graph.records)


def test_contracts_do_not_import_private_science_modules() -> None:
    graph = _build_import_graph()
    offenders = sorted(
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if record.source_module.startswith("phospy.contracts")
        if _is_private_science_module(record.target)
    )

    assert offenders == []


def test_contracts_do_not_import_science_implementation_modules() -> None:
    graph = _build_import_graph()

    offenders = sorted(
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if record.source_module.startswith("phospy.contracts")
        and record.target.startswith("phospy.science")
        if _is_forbidden_contracts_science_implementation_module(record.target)
    )

    assert offenders == []


def test_contracts_to_science_imports_match_public_facade_rules() -> None:
    graph = _build_import_graph()
    offenders = sorted(
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if record.source_module.startswith("phospy.contracts")
        and record.target.startswith("phospy.science")
        if not _is_public_contracts_science_facade_import(
            record.target,
            module_paths=graph.module_paths,
        )
    )

    assert offenders == []


def test_contracts_science_facade_role_markers_are_used_and_valid() -> None:
    graph = _build_import_graph()
    observed = {
        record.target
        for record in graph.records
        if record.source_module.startswith("phospy.contracts")
        and record.target.startswith("phospy.science")
    }
    problems = _contracts_science_facade_role_marker_problems(
        module_paths=graph.module_paths,
        observed_contracts_science_targets=observed,
    )

    assert problems == []


def test_science_package_does_not_import_contracts() -> None:
    graph = _build_import_graph()
    offenders = sorted(
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if record.source_module.startswith("phospy.science")
        and record.target.startswith("phospy.contracts")
    )

    assert offenders == []


def test_contracts_science_facade_rule_allows_config_prefix_without_marker() -> None:
    assert _is_public_contracts_science_facade_import(
        "phospy.science.configs.preprocessing.validation",
        module_paths={},
    )


def test_contracts_science_facade_rule_rejects_unmarked_science_module(
    tmp_path: Path,
) -> None:
    module_name = "phospy.science.public_unmarked.models"
    module_paths = _facade_rule_fixture_module_paths(
        tmp_path,
        module_name=module_name,
        source='"""Public-looking but unmarked science module."""\n',
    )

    assert not _is_public_contracts_science_facade_import(
        module_name,
        module_paths=module_paths,
    )


def test_contracts_science_facade_rule_accepts_marked_public_owner(
    tmp_path: Path,
) -> None:
    module_name = "phospy.science.public_marked.models"
    module_paths = _facade_rule_fixture_module_paths(
        tmp_path,
        module_name=module_name,
        source=(
            '"""Marked public science model owner."""\n'
            f"{CONTRACTS_SCIENCE_FACADE_ROLE_MARKER} = "
            '"science_owned_public_model"\n'
        ),
    )

    assert _is_public_contracts_science_facade_import(
        module_name,
        module_paths=module_paths,
    )


@pytest.mark.parametrize(
    ("symbol_source", "expected_symbol"),
    (
        ("class WorkflowExecutor:\n    pass\n", "WorkflowExecutor"),
        ("class PublicRequestValidator:\n    pass\n", "PublicRequestValidator"),
        (
            "class KinaseLibraryResourceLoader:\n    pass\n",
            "KinaseLibraryResourceLoader",
        ),
        ("class ReferenceResourceBuilder:\n    pass\n", "ReferenceResourceBuilder"),
        (
            "class ProteinAwarePreparationStage:\n    pass\n",
            "ProteinAwarePreparationStage",
        ),
    ),
)
def test_contracts_science_facade_role_rejects_implementation_declarations(
    tmp_path: Path,
    symbol_source: str,
    expected_symbol: str,
) -> None:
    module_name = "phospy.science.public_marked.models"
    module_paths = _facade_rule_fixture_module_paths(
        tmp_path,
        module_name=module_name,
        source=(
            '"""Marked module with an implementation declaration."""\n'
            f"{CONTRACTS_SCIENCE_FACADE_ROLE_MARKER} = "
            '"science_owned_public_model"\n\n'
            f"{symbol_source}"
        ),
    )

    assert not _is_public_contracts_science_facade_import(
        module_name,
        module_paths=module_paths,
    )
    problems = _contracts_science_facade_role_marker_problems(
        module_paths=module_paths,
        observed_contracts_science_targets={module_name},
    )

    assert len(problems) == 1
    assert module_name in problems[0]
    assert "science_owned_public_model" in problems[0]
    assert expected_symbol in problems[0]
    assert ":4:" in problems[0]


def test_contracts_science_facade_role_rejects_loader_function(
    tmp_path: Path,
) -> None:
    module_name = "phospy.science.public_marked.models"
    module_paths = _facade_rule_fixture_module_paths(
        tmp_path,
        module_name=module_name,
        source=(
            '"""Marked module with an executable loader function."""\n'
            f"{CONTRACTS_SCIENCE_FACADE_ROLE_MARKER} = "
            '"science_owned_public_model"\n\n'
            "def load_reference_resource() -> object:\n"
            "    return object()\n"
        ),
    )

    assert not _is_public_contracts_science_facade_import(
        module_name,
        module_paths=module_paths,
    )
    problems = _contracts_science_facade_role_marker_problems(
        module_paths=module_paths,
        observed_contracts_science_targets={module_name},
    )

    assert problems == [
        f"{module_name}:4: {CONTRACTS_SCIENCE_FACADE_ROLE_MARKER} role "
        "'science_owned_public_model' is not role-pure: offending symbol "
        "'load_reference_resource'"
    ]


def test_contracts_science_facade_role_rejects_implementation_import(
    tmp_path: Path,
) -> None:
    module_name = "phospy.science.public_marked.models"
    module_paths = _facade_rule_fixture_module_paths(
        tmp_path,
        module_name=module_name,
        source=(
            '"""Marked module importing validation implementation."""\n'
            f"{CONTRACTS_SCIENCE_FACADE_ROLE_MARKER} = "
            '"science_owned_public_model"\n\n'
            "from phospy.science.references.validation.bundle import "
            "ReferenceBundleValidator\n"
        ),
    )

    assert not _is_public_contracts_science_facade_import(
        module_name,
        module_paths=module_paths,
    )
    problems = _contracts_science_facade_role_marker_problems(
        module_paths=module_paths,
        observed_contracts_science_targets={module_name},
    )

    assert problems == [
        f"{module_name}:4: {CONTRACTS_SCIENCE_FACADE_ROLE_MARKER} role "
        "'science_owned_public_model' is not role-pure: offending import symbol "
        "'ReferenceBundleValidator'",
    ]


@pytest.mark.parametrize(
    "module_name",
    (
        "phospy.science.public_marked._private_model",
        "phospy.science.datasets.builders.contracts",
        "phospy.science.datasets.construction.service",
        "phospy.science.public_marked.executor",
        "phospy.science.public_marked.validation",
        "phospy.science.references.builder",
        "phospy.science.references.validation.manifest_schema",
        "phospy.science.public_marked.interpreter",
    ),
)
def test_contracts_science_facade_rule_rejects_forbidden_modules_even_when_marked(
    tmp_path: Path,
    module_name: str,
) -> None:
    module_paths = _facade_rule_fixture_module_paths(
        tmp_path,
        module_name=module_name,
        source=(
            '"""Forbidden implementation module with a stale marker."""\n'
            f"{CONTRACTS_SCIENCE_FACADE_ROLE_MARKER} = "
            '"science_owned_public_model"\n'
        ),
    )

    assert not _is_public_contracts_science_facade_import(
        module_name,
        module_paths=module_paths,
    )


def test_contracts_science_facade_role_marker_audit_detects_stale_marker(
    tmp_path: Path,
) -> None:
    module_name = "phospy.science.stale_public.models"
    module_paths = _facade_rule_fixture_module_paths(
        tmp_path,
        module_name=module_name,
        source=(
            '"""No contracts module imports this science facade."""\n'
            f"{CONTRACTS_SCIENCE_FACADE_ROLE_MARKER} = "
            '"science_owned_public_model"\n'
        ),
    )

    problems = _contracts_science_facade_role_marker_problems(
        module_paths=module_paths,
        observed_contracts_science_targets=frozenset(),
    )

    assert problems == [f"{module_name}: unused {CONTRACTS_SCIENCE_FACADE_ROLE_MARKER}"]


def test_contracts_science_facade_role_marker_audit_detects_invalid_role(
    tmp_path: Path,
) -> None:
    module_name = "phospy.science.invalid_public.models"
    module_paths = _facade_rule_fixture_module_paths(
        tmp_path,
        module_name=module_name,
        source=(
            '"""Marked with an unsupported role."""\n'
            f"{CONTRACTS_SCIENCE_FACADE_ROLE_MARKER} = "
            '"public_model"\n'
        ),
    )

    assert not _is_public_contracts_science_facade_import(
        module_name,
        module_paths=module_paths,
    )
    problems = _contracts_science_facade_role_marker_problems(
        module_paths=module_paths,
        observed_contracts_science_targets={module_name},
    )

    assert problems == [
        f"{module_name}:2: invalid {CONTRACTS_SCIENCE_FACADE_ROLE_MARKER} "
        "'public_model'"
    ]


def test_contracts_science_facade_rule_accepts_real_role_pure_science_import() -> None:
    graph = _build_import_graph()

    assert _is_public_contracts_science_facade_import(
        "phospy.science.evidence.dataset_resolution.models",
        module_paths=graph.module_paths,
    )


def test_workflows_do_not_import_unresolved_peptide_evidence_models() -> None:
    graph = _build_import_graph()
    offenders = sorted(
        f"{record.source_module}:{record.line} -> {record.target}"
        for record in graph.records
        if record.source_module.startswith("phospy.workflows")
        if record.target.startswith("phospy.science.evidence")
    )

    assert offenders == []


def test_import_extractor_includes_all_static_import_forms() -> None:
    source = """
from typing import TYPE_CHECKING
import phospy.science.datasets.models
from phospy.validation.datasets import preprocessing
from .science import datasets
import importlib
importlib.import_module("phospy.io.bundles")
__import__("phospy.workflows.kinase")
if TYPE_CHECKING:
    from phospy.science.references.models import ReferenceBundle
def _local_imports():
    import phospy.science.configs.dataset
    from .science.signalomes import models
    importlib.import_module("phospy.io.bundles.signalome")
"""
    path = PACKAGE_ROOT / "__init__.py"
    imports = tuple(
        imported
        for imported, _line in _imported_modules("phospy", path, ast.parse(source))
    )

    assert "phospy.science.datasets.models" in imports
    assert "phospy.validation.datasets" in imports
    assert "phospy.validation.datasets.preprocessing" in imports
    assert "phospy.science" in imports
    assert "phospy.science.datasets" in imports
    assert "phospy.science.references.models" in imports
    assert "phospy.science.references.models.ReferenceBundle" in imports
    assert "phospy.science.configs.dataset" in imports
    assert "phospy.science.signalomes" in imports
    assert "phospy.science.signalomes.models" in imports
    assert "phospy.io.bundles" in imports
    assert "phospy.io.bundles.signalome" in imports
    assert "phospy.workflows.kinase" in imports


def test_import_extractor_records_ast_line_numbers() -> None:
    source = """
import phospy.science.datasets.models
from phospy.validation.datasets import preprocessing
import importlib
importlib.import_module("phospy.io.bundles")
"""
    path = PACKAGE_ROOT / "__init__.py"
    imports = tuple(_imported_modules("phospy", path, ast.parse(source)))

    assert imports == (
        ("phospy.science.datasets.models", 2),
        ("phospy.validation.datasets", 3),
        ("phospy.validation.datasets.preprocessing", 3),
        ("importlib", 4),
        ("phospy.io.bundles", 5),
    )


def test_public_import_routes_load_in_clean_interpreters() -> None:
    code = "\n".join(
        [
            "from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow",
            "from phospy.api import ReferenceBundleBuilder",
            "from phospy.contracts.requests import DatasetBuildRequest",
            "from phospy.errors import PhosPyError",
            "from phospy.science.references.models import ReferenceBundle",
            "from phospy.science.datasets.builders.public import "
            "AnalysisReadyDatasetBuilder as ScienceDatasetBuilder",
            "from phospy.validation.references.bundle import ReferenceBundleValidator",
            "assert AnalysisReadyDatasetBuilder",
            "assert ReferenceBundleBuilder",
            "assert DatasetBuildRequest",
            "assert PhosPyError",
            "assert ReferenceBundle",
            "assert ScienceDatasetBuilder",
            "assert ReferenceBundleValidator",
            "assert KinaseWorkflow",
        ]
    )
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join(
        [
            str(SRC_ROOT),
            environment.get("PYTHONPATH", ""),
        ]
    )

    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=PROJECT_ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr


def _build_import_graph() -> ImportGraph:
    paths = {_module_name(path): path for path in PACKAGE_ROOT.rglob("*.py")}
    modules = frozenset(paths)
    records: list[ImportRecord] = []
    module_edges: dict[str, frozenset[str]] = {}
    for module_name, path in paths.items():
        resolved_edges: set[str] = set()
        tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
        for imported, line_number in _imported_modules(module_name, path, tree):
            if not imported.startswith("phospy"):
                continue
            resolved = _resolve_module(imported, modules)
            if resolved is None or resolved == module_name:
                continue
            resolved_edges.add(resolved)
            records.append(
                ImportRecord(
                    source_module=module_name,
                    source_path=path,
                    target=resolved,
                    line=line_number,
                )
            )
        module_edges[module_name] = frozenset(resolved_edges)
    package_edges = frozenset(_package_edges(module_edges))
    return ImportGraph(
        modules=modules,
        module_paths=paths,
        module_edges=module_edges,
        package_edges=package_edges,
        records=tuple(records),
    )


def _module_name(path: Path) -> str:
    relative = path.relative_to(SRC_ROOT).with_suffix("")
    parts = relative.parts
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


def _imported_modules(
    module_name: str,
    path: Path,
    tree: ast.AST,
) -> Iterable[tuple[str, int]]:
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                yield alias.name, node.lineno
        elif isinstance(node, ast.ImportFrom):
            imported = _resolve_import_from(module_name, path, node)
            if imported is None:
                continue
            yield imported, node.lineno
            for alias in node.names:
                yield (
                    f"{imported}.{alias.name}" if imported else alias.name,
                    node.lineno,
                )
        elif isinstance(node, ast.Call):
            imported = _static_dynamic_import_target(node)
            if imported is not None:
                yield imported, node.lineno


def _resolve_import_from(
    module_name: str,
    path: Path,
    node: ast.ImportFrom,
) -> str | None:
    if node.level == 0:
        return node.module
    parts = module_name.split(".")
    package_parts = parts if path.name == "__init__.py" else parts[:-1]
    if node.level > 1:
        package_parts = package_parts[: -(node.level - 1)]
    module_parts = [] if node.module is None else node.module.split(".")
    return ".".join((*package_parts, *module_parts))


def _static_dynamic_import_target(node: ast.Call) -> str | None:
    is_importlib_call = (
        isinstance(node.func, ast.Attribute) and node.func.attr == "import_module"
    )
    is_dunder_import = isinstance(node.func, ast.Name) and node.func.id == "__import__"
    if not (is_importlib_call or is_dunder_import):
        return None
    if not node.args:
        return None
    first_arg = node.args[0]
    if not isinstance(first_arg, ast.Constant) or not isinstance(first_arg.value, str):
        return None
    return first_arg.value


def _resolve_module(imported: str, modules: frozenset[str]) -> str | None:
    parts = imported.split(".")
    for index in range(len(parts), 1, -1):
        candidate = ".".join(parts[:index])
        if candidate in modules:
            return candidate
    return None


def _is_private_science_module(module_name: str) -> bool:
    if module_name in PRIVATE_SCIENCE_MODULES:
        return True
    if any(
        module_name.startswith(f"{private_module}.")
        for private_module in PRIVATE_SCIENCE_MODULES
    ):
        return True
    if not module_name.startswith("phospy.science."):
        return False
    return any(part.startswith("_") for part in module_name.split(".")[2:])


def _is_forbidden_contracts_science_implementation_module(module_name: str) -> bool:
    if not module_name.startswith("phospy.science."):
        return False
    science_parts = module_name.split(".")[2:]
    if any(
        part in CONTRACTS_FORBIDDEN_SCIENCE_STRUCTURE_SEGMENTS for part in science_parts
    ):
        return True
    if any(
        part.endswith("_builder") or part.endswith("_builders")
        for part in science_parts
    ):
        return True
    if any(part.endswith("_executor") for part in science_parts):
        return True
    if any(part in {"internal_view", "internal_views"} for part in science_parts):
        return True
    if science_parts[-1] == "internal_frame_store":
        return True
    if _is_reference_builder_or_validation_module(module_name):
        return True
    if _is_validation_implementation_module(module_name):
        return True
    return _is_workflow_implementation_module(module_name)


def _is_reference_builder_or_validation_module(module_name: str) -> bool:
    if not _module_is_or_under(module_name, "phospy.science.references"):
        return False
    reference_parts = module_name.split(".")[3:]
    return any(
        part in {"builder", "builders", "validation", "validators"}
        or part.endswith("_builder")
        for part in reference_parts
    )


def _is_workflow_implementation_module(module_name: str) -> bool:
    module_leaf = module_name.rsplit(".", maxsplit=1)[-1]
    return (
        module_leaf in CONTRACTS_FORBIDDEN_SCIENCE_WORKFLOW_LEAVES
        or module_leaf.endswith("_assembler")
        or module_leaf.endswith("_interpreter")
        or module_leaf.endswith("_orchestrator")
        or module_leaf.endswith("_runner")
        or module_leaf.endswith("_workflow")
    )


def _is_validation_implementation_module(module_name: str) -> bool:
    if _module_is_or_under(module_name, "phospy.science.configs"):
        return False
    science_parts = module_name.split(".")[2:]
    return any(
        part in {"validation", "validator", "validators", "resolved_validator"}
        or part.endswith("_validation")
        or part.endswith("_validator")
        for part in science_parts
    )


def _is_public_contracts_science_facade_import(
    module_name: str,
    *,
    module_paths: Mapping[str, Path],
) -> bool:
    if _is_private_science_module(module_name):
        return False
    if _is_forbidden_contracts_science_implementation_module(module_name):
        return False
    if any(
        _module_is_or_under(module_name, allowed_prefix)
        for allowed_prefix in CONTRACTS_PUBLIC_SCIENCE_PREFIXES
    ):
        return True
    marker = _contracts_science_facade_role_marker(
        module_name,
        module_paths=module_paths,
    )
    return (
        marker is not None
        and marker.role in CONTRACTS_PUBLIC_SCIENCE_FACADE_ROLES
        and not _contracts_science_facade_role_purity_problems(marker)
    )


def _contracts_science_facade_role_marker(
    module_name: str,
    *,
    module_paths: Mapping[str, Path],
) -> ContractsScienceFacadeRoleMarker | None:
    path = module_paths.get(module_name)
    if path is None:
        return None
    tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
    for node in tree.body:
        if isinstance(node, ast.Assign):
            if not any(
                isinstance(target, ast.Name)
                and target.id == CONTRACTS_SCIENCE_FACADE_ROLE_MARKER
                for target in node.targets
            ):
                continue
            return ContractsScienceFacadeRoleMarker(
                module_name=module_name,
                source_path=path,
                role=_literal_string_value(node.value),
                line=node.lineno,
            )
        if (
            isinstance(node, ast.AnnAssign)
            and isinstance(node.target, ast.Name)
            and node.target.id == CONTRACTS_SCIENCE_FACADE_ROLE_MARKER
        ):
            return ContractsScienceFacadeRoleMarker(
                module_name=module_name,
                source_path=path,
                role=_literal_string_value(node.value),
                line=node.lineno,
            )
    return None


def _literal_string_value(node: ast.AST | None) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    return None


def _contracts_science_facade_role_markers(
    *,
    module_paths: Mapping[str, Path],
) -> tuple[ContractsScienceFacadeRoleMarker, ...]:
    markers: list[ContractsScienceFacadeRoleMarker] = []
    for module_name in module_paths:
        if not module_name.startswith("phospy.science"):
            continue
        marker = _contracts_science_facade_role_marker(
            module_name,
            module_paths=module_paths,
        )
        if marker is not None:
            markers.append(marker)
    return tuple(sorted(markers, key=lambda marker: marker.module_name))


def _contracts_science_facade_role_marker_problems(
    *,
    module_paths: Mapping[str, Path],
    observed_contracts_science_targets: Iterable[str],
) -> list[str]:
    observed = set(observed_contracts_science_targets)
    problems: list[str] = []
    for marker in _contracts_science_facade_role_markers(module_paths=module_paths):
        if marker.role not in CONTRACTS_PUBLIC_SCIENCE_FACADE_ROLES:
            problems.append(
                f"{marker.module_name}:{marker.line}: invalid "
                f"{CONTRACTS_SCIENCE_FACADE_ROLE_MARKER} {marker.role!r}"
            )
        if _is_private_science_module(
            marker.module_name
        ) or _is_forbidden_contracts_science_implementation_module(marker.module_name):
            problems.append(
                f"{marker.module_name}:{marker.line}: forbidden "
                f"{CONTRACTS_SCIENCE_FACADE_ROLE_MARKER}"
            )
        else:
            problems.extend(_contracts_science_facade_role_purity_problems(marker))
        if marker.module_name not in observed:
            problems.append(
                f"{marker.module_name}: unused {CONTRACTS_SCIENCE_FACADE_ROLE_MARKER}"
            )
    return problems


def _contracts_science_facade_role_purity_problems(
    marker: ContractsScienceFacadeRoleMarker,
) -> list[str]:
    if marker.role not in CONTRACTS_PUBLIC_SCIENCE_FACADE_ROLES:
        return []
    tree = ast.parse(
        marker.source_path.read_text(encoding="utf-8-sig"),
        filename=str(marker.source_path),
    )
    problems: list[str] = []
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and _is_forbidden_facade_declaration_name(
            node.name,
        ):
            problems.append(
                _format_facade_role_purity_problem(
                    marker,
                    line=node.lineno,
                    detail=f"offending symbol {node.name!r}",
                )
            )
        elif isinstance(
            node,
            ast.FunctionDef | ast.AsyncFunctionDef,
        ) and _is_forbidden_facade_function_name(node.name):
            problems.append(
                _format_facade_role_purity_problem(
                    marker,
                    line=node.lineno,
                    detail=f"offending symbol {node.name!r}",
                )
            )
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(
                    target,
                    ast.Name,
                ) and _is_forbidden_facade_declaration_name(target.id):
                    problems.append(
                        _format_facade_role_purity_problem(
                            marker,
                            line=node.lineno,
                            detail=f"offending symbol {target.id!r}",
                        )
                    )
        elif (
            isinstance(node, ast.AnnAssign)
            and isinstance(node.target, ast.Name)
            and _is_forbidden_facade_declaration_name(node.target.id)
        ):
            problems.append(
                _format_facade_role_purity_problem(
                    marker,
                    line=node.lineno,
                    detail=f"offending symbol {node.target.id!r}",
                )
            )

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if _is_forbidden_contracts_science_implementation_module(alias.name):
                    problems.append(
                        _format_facade_role_purity_problem(
                            marker,
                            line=node.lineno,
                            detail=f"offending import {alias.name!r}",
                        )
                    )
                imported_name = alias.asname or alias.name.rsplit(".", maxsplit=1)[-1]
                if _is_forbidden_facade_declaration_name(imported_name):
                    problems.append(
                        _format_facade_role_purity_problem(
                            marker,
                            line=node.lineno,
                            detail=f"offending import symbol {imported_name!r}",
                        )
                    )
        elif isinstance(node, ast.ImportFrom):
            for alias in node.names:
                imported = _resolve_import_from(
                    marker.module_name, marker.source_path, node
                )
                if (
                    alias.name == "*"
                    and imported is not None
                    and _is_forbidden_contracts_science_implementation_module(imported)
                ):
                    problems.append(
                        _format_facade_role_purity_problem(
                            marker,
                            line=node.lineno,
                            detail=f"offending import {imported!r}.*",
                        )
                    )
                imported_name = alias.asname or alias.name
                if _is_forbidden_facade_declaration_name(imported_name):
                    problems.append(
                        _format_facade_role_purity_problem(
                            marker,
                            line=node.lineno,
                            detail=f"offending import symbol {imported_name!r}",
                        )
                    )
    return sorted(set(problems))


def _format_facade_role_purity_problem(
    marker: ContractsScienceFacadeRoleMarker,
    *,
    line: int,
    detail: str,
) -> str:
    return (
        f"{marker.module_name}:{line}: "
        f"{CONTRACTS_SCIENCE_FACADE_ROLE_MARKER} role {marker.role!r} "
        f"is not role-pure: {detail}"
    )


def _is_forbidden_facade_function_name(name: str) -> bool:
    return _is_public_facade_name(name) and name.startswith(
        CONTRACTS_FORBIDDEN_FACADE_FUNCTION_PREFIXES
    )


def _is_forbidden_facade_declaration_name(name: str) -> bool:
    return _is_public_facade_name(name) and (
        name.endswith(CONTRACTS_FORBIDDEN_FACADE_DECLARATION_SUFFIXES)
        or any(
            fragment in name
            for fragment in CONTRACTS_FORBIDDEN_FACADE_DECLARATION_FRAGMENTS
        )
    )


def _is_public_facade_name(name: str) -> bool:
    return not name.startswith("_")


def _facade_rule_fixture_module_paths(
    tmp_path: Path,
    *,
    module_name: str,
    source: str,
) -> dict[str, Path]:
    path = tmp_path / f"{module_name.replace('.', '_')}.py"
    path.write_text(source, encoding="utf-8")
    return {module_name: path}


def _module_is_or_under(module_name: str, prefix: str) -> bool:
    return module_name == prefix or module_name.startswith(f"{prefix}.")


def _package_edges(
    module_edges: dict[str, frozenset[str]],
) -> set[tuple[str, str]]:
    edges: set[tuple[str, str]] = set()
    for source, targets in module_edges.items():
        source_package = _package_name(source)
        if source_package is None:
            continue
        for target in targets:
            target_package = _package_name(target)
            if target_package is None or target_package == source_package:
                continue
            edges.add((source_package, target_package))
    return edges


def _package_name(module_name: str) -> str | None:
    parts = module_name.split(".")
    if len(parts) < 2 or parts[0] != "phospy":
        return None
    return ".".join(parts[:2])


def _package_adjacency(
    package_edges: Iterable[tuple[str, str]],
) -> dict[str, frozenset[str]]:
    adjacency: dict[str, set[str]] = {}
    for source, target in package_edges:
        adjacency.setdefault(source, set()).add(target)
        adjacency.setdefault(target, set())
    return {source: frozenset(targets) for source, targets in adjacency.items()}


def _strongly_connected_components(
    graph: dict[str, frozenset[str]],
) -> list[tuple[str, ...]]:
    index = 0
    stack: list[str] = []
    on_stack: set[str] = set()
    indices: dict[str, int] = {}
    lowlinks: dict[str, int] = {}
    components: list[tuple[str, ...]] = []

    def strong_connect(node: str) -> None:
        nonlocal index
        indices[node] = index
        lowlinks[node] = index
        index += 1
        stack.append(node)
        on_stack.add(node)

        for target in graph.get(node, ()):
            if target not in indices:
                strong_connect(target)
                lowlinks[node] = min(lowlinks[node], lowlinks[target])
            elif target in on_stack:
                lowlinks[node] = min(lowlinks[node], indices[target])

        if lowlinks[node] != indices[node]:
            return
        component: list[str] = []
        while True:
            target = stack.pop()
            on_stack.remove(target)
            component.append(target)
            if target == node:
                break
        components.append(tuple(component))

    for node in graph:
        if node not in indices:
            strong_connect(node)
    return components


def _format_cycles(cycles: list[tuple[str, ...]]) -> str:
    if not cycles:
        return ""
    return "\n\n".join("\n".join(cycle) for cycle in cycles)


def _format_package_edges(
    edges: Iterable[tuple[str, str]],
    records: tuple[ImportRecord, ...],
) -> str:
    lines: list[str] = []
    for source, target in edges:
        lines.append(f"{source} -> {target}")
        examples = [
            record
            for record in records
            if _package_name(record.source_module) == source
            and _package_name(record.target) == target
        ][:8]
        for record in examples:
            lines.append(
                f"  {record.source_path.relative_to(PROJECT_ROOT)}:{record.line}"
            )
    return "\n".join(lines)
