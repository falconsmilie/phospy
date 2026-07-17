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
        ("phospy.api", "phospy.contracts"),
        ("phospy.api", "phospy.errors"),
        ("phospy.api", "phospy.io"),
        ("phospy.api", "phospy.science"),
        ("phospy.api", "phospy.tables"),
        ("phospy.api", "phospy.validation"),
        ("phospy.api", "phospy.workflows"),
        ("phospy.contracts", "phospy.errors"),
        ("phospy.contracts", "phospy.frames"),
        ("phospy.contracts", "phospy.policies"),
        ("phospy.contracts", "phospy.provenance"),
        ("phospy.contracts", "phospy.science"),
        ("phospy.contracts", "phospy.tables"),
        ("phospy.frames", "phospy.errors"),
        ("phospy.io", "phospy.contracts"),
        ("phospy.io", "phospy.errors"),
        ("phospy.io", "phospy.provenance"),
        ("phospy.io", "phospy.science"),
        ("phospy.io", "phospy.validation"),
        ("phospy.policies", "phospy.errors"),
        ("phospy.provenance", "phospy.errors"),
        ("phospy.release", "phospy.provenance"),
        ("phospy.science", "phospy.errors"),
        ("phospy.science", "phospy.frames"),
        ("phospy.science", "phospy.policies"),
        ("phospy.science", "phospy.provenance"),
        ("phospy.tables", "phospy.errors"),
        ("phospy.tables", "phospy.frames"),
        ("phospy.tables", "phospy.science"),
        ("phospy.validation", "phospy.contracts"),
        ("phospy.validation", "phospy.errors"),
        ("phospy.validation", "phospy.frames"),
        ("phospy.validation", "phospy.provenance"),
        ("phospy.validation", "phospy.science"),
        ("phospy.workflows", "phospy.contracts"),
        ("phospy.workflows", "phospy.errors"),
        ("phospy.workflows", "phospy.provenance"),
        ("phospy.workflows", "phospy.science"),
        ("phospy.workflows", "phospy.tables"),
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


def test_import_extractor_includes_all_static_import_forms() -> None:
    source = """
import phospy.science.datasets.models
from phospy.validation.datasets import preprocessing
from .science import datasets
import importlib
importlib.import_module("phospy.io.bundles")
__import__("phospy.workflows.kinase")
"""
    path = PACKAGE_ROOT / "__init__.py"
    imports = tuple(_imported_modules("phospy", path, ast.parse(source)))

    assert "phospy.science.datasets.models" in imports
    assert "phospy.validation.datasets" in imports
    assert "phospy.validation.datasets.preprocessing" in imports
    assert "phospy.science" in imports
    assert "phospy.science.datasets" in imports
    assert "phospy.io.bundles" in imports
    assert "phospy.workflows.kinase" in imports


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
        for imported in _imported_modules(module_name, path, tree):
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
                    line=1,
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
) -> Iterable[str]:
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                yield alias.name
        elif isinstance(node, ast.ImportFrom):
            imported = _resolve_import_from(module_name, path, node)
            if imported is None:
                continue
            yield imported
            for alias in node.names:
                yield f"{imported}.{alias.name}" if imported else alias.name
        elif isinstance(node, ast.Call):
            imported = _static_dynamic_import_target(node)
            if imported is not None:
                yield imported


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
            lines.append(f"  {record.source_path.relative_to(PROJECT_ROOT)}")
    return "\n".join(lines)
