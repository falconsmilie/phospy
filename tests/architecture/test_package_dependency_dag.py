from __future__ import annotations

import ast
import os
import subprocess
import sys
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

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

CONTRACTS_PUBLIC_SCIENCE_MODULES = frozenset(
    {
        "phospy.science.activities.models",
        "phospy.science.datasets.models",
        "phospy.science.datasets.preprocessing.batch_correction_models",
        "phospy.science.datasets.preprocessing.protein_aware_preparation",
        "phospy.science.design.contrast_helpers",
        "phospy.science.design.models",
        "phospy.science.differential.models",
        "phospy.science.differential.policy_models",
        "phospy.science.enrichment.models",
        "phospy.science.evidence.dataset_resolution.contracts",
        "phospy.science.prediction.models",
        "phospy.science.references.kinase_library",
        "phospy.science.references.models",
        "phospy.science.result_caveats",
        "phospy.science.signalomes.constants",
        "phospy.science.signalomes.models",
        "phospy.science.tables.kinase",
        "phospy.science.tables.signalome",
        "phospy.science.transformations.models",
    }
)

CONTRACTS_FORBIDDEN_SCIENCE_PREFIXES = frozenset(
    {
        "phospy.science.datasets.builders",
        "phospy.science.datasets.construction",
        "phospy.science.references.builder",
        "phospy.science.references.validation",
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
    module_edges: dict[str, frozenset[str]]
    package_edges: frozenset[tuple[str, str]]
    records: tuple[ImportRecord, ...]


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
        if not _is_public_contracts_science_facade_import(record.target)
    )

    assert offenders == []


def test_contracts_public_science_module_set_has_no_stale_entries() -> None:
    graph = _build_import_graph()
    observed = {
        record.target
        for record in graph.records
        if record.source_module.startswith("phospy.contracts")
        and record.target.startswith("phospy.science")
    }
    stale_modules = sorted(CONTRACTS_PUBLIC_SCIENCE_MODULES - observed)

    assert stale_modules == []


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
    if any(
        _module_is_or_under(module_name, forbidden_prefix)
        for forbidden_prefix in CONTRACTS_FORBIDDEN_SCIENCE_PREFIXES
    ):
        return True
    module_leaf = module_name.rsplit(".", maxsplit=1)[-1]
    if module_leaf in {"executor", "execution", "internal_view"}:
        return True
    if module_leaf.endswith("_executor"):
        return True
    return _is_validation_implementation_module(module_name)


def _is_validation_implementation_module(module_name: str) -> bool:
    if _module_is_or_under(module_name, "phospy.science.configs"):
        return False
    module_leaf = module_name.rsplit(".", maxsplit=1)[-1]
    return (
        module_leaf in {"validation", "validator", "resolved_validator"}
        or module_leaf.endswith("_validation")
        or module_leaf.endswith("_validator")
    )


def _is_public_contracts_science_facade_import(module_name: str) -> bool:
    if _is_private_science_module(module_name):
        return False
    if _is_forbidden_contracts_science_implementation_module(module_name):
        return False
    return module_name in CONTRACTS_PUBLIC_SCIENCE_MODULES or any(
        _module_is_or_under(module_name, allowed_prefix)
        for allowed_prefix in CONTRACTS_PUBLIC_SCIENCE_PREFIXES
    )


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
