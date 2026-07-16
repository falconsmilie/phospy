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
TARGET_DOMAINS = frozenset(
    {
        "api",
        "contracts",
        "errors",
        "io",
        "provenance",
        "science",
        "validation",
        "workflows",
    }
)


@dataclass(frozen=True, slots=True)
class ImportGraph:
    modules: frozenset[str]
    edges: dict[str, frozenset[str]]
    paths: dict[str, Path]


def test_package_dependency_graph_has_no_static_cycles() -> None:
    graph = _build_import_graph()
    cycles = [
        tuple(sorted(component))
        for component in _strongly_connected_components(graph.edges)
        if len(component) > 1
    ]

    assert cycles == [], _format_cycles(cycles)


def test_errors_package_is_a_leaf() -> None:
    graph = _build_import_graph()
    offenders = _forbidden_edges(
        graph,
        source_roots={"phospy.errors"},
        forbidden_targets={"phospy"},
        allowed_targets={"phospy.errors"},
    )

    assert offenders == [], _format_edges(offenders)


def test_contracts_do_not_depend_on_workflow_or_validation_packages() -> None:
    graph = _build_import_graph()
    offenders = _forbidden_edges(
        graph,
        source_roots={"phospy.contracts"},
        forbidden_targets={"phospy.validation", "phospy.workflows"},
    )

    assert offenders == [], _format_edges(offenders)


def test_science_does_not_depend_on_workflow_or_concrete_io_packages() -> None:
    graph = _build_import_graph()
    offenders = _forbidden_edges(
        graph,
        source_roots={"phospy.science"},
        forbidden_targets={"phospy.io", "phospy.workflows"},
    )

    assert offenders == [], _format_edges(offenders)


def test_package_dependency_graph_is_reviewed_in_ci() -> None:
    graph = _build_import_graph()
    package_edges = _package_edges(graph)

    assert package_edges, "package graph is unexpectedly empty"
    assert "phospy.science -> phospy.contracts" in package_edges
    assert "phospy.workflows -> phospy.science" in package_edges
    assert all(
        edge not in package_edges
        for edge in (
            "phospy.errors -> phospy.science",
            "phospy.contracts -> phospy.validation",
            "phospy.contracts -> phospy.workflows",
            "phospy.science -> phospy.io",
            "phospy.science -> phospy.workflows",
        )
    ), "\n".join(package_edges)


def test_public_import_routes_load_in_clean_interpreters() -> None:
    code = "\n".join(
        [
            "from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow",
            "from phospy.api import ReferenceBundleBuilder",
            "from phospy.contracts.requests import DatasetBuildRequest",
            "from phospy.errors import PhosPyError",
            "from phospy.science.references.models import ReferenceBundle",
            "from phospy.validation.references.bundle import ReferenceBundleValidator",
            "assert AnalysisReadyDatasetBuilder",
            "assert ReferenceBundleBuilder",
            "assert DatasetBuildRequest",
            "assert PhosPyError",
            "assert ReferenceBundle",
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
    paths = {
        _module_name(path): path
        for path in PACKAGE_ROOT.rglob("*.py")
        if _domain_name(path) in TARGET_DOMAINS
    }
    modules = frozenset(paths)
    edges: dict[str, frozenset[str]] = {}
    for module_name, path in paths.items():
        resolved_edges = {
            resolved
            for imported in _imported_modules(module_name, path)
            if (resolved := _resolve_module(imported, modules)) is not None
            and _top_level_package(resolved) in TARGET_DOMAINS
        }
        edges[module_name] = frozenset(resolved_edges)
    return ImportGraph(modules=modules, edges=edges, paths=paths)


def _module_name(path: Path) -> str:
    relative = path.relative_to(SRC_ROOT).with_suffix("")
    parts = relative.parts
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


def _domain_name(path: Path) -> str | None:
    relative = path.relative_to(PACKAGE_ROOT)
    if len(relative.parts) < 2:
        return None
    return relative.parts[0]


def _imported_modules(module_name: str, path: Path) -> Iterable[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith("phospy"):
                    yield alias.name
        elif isinstance(node, ast.ImportFrom):
            imported = _resolve_import_from(module_name, path, node)
            if imported is None or not imported.startswith("phospy"):
                continue
            yield imported
            for alias in node.names:
                yield f"{imported}.{alias.name}"


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


def _resolve_module(imported: str, modules: frozenset[str]) -> str | None:
    parts = imported.split(".")
    for index in range(len(parts), 1, -1):
        candidate = ".".join(parts[:index])
        if candidate in modules:
            return candidate
    return None


def _top_level_package(module_name: str) -> str | None:
    parts = module_name.split(".")
    if len(parts) < 2:
        return None
    return parts[1]


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
            if target not in graph:
                continue
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


def _forbidden_edges(
    graph: ImportGraph,
    *,
    source_roots: set[str],
    forbidden_targets: set[str],
    allowed_targets: set[str] | None = None,
) -> list[str]:
    allowed_targets = allowed_targets or set()
    offenders: list[str] = []
    for source, targets in graph.edges.items():
        if not any(_is_under(source, root) for root in source_roots):
            continue
        for target in sorted(targets):
            if any(_is_under(target, allowed) for allowed in allowed_targets):
                continue
            if any(_is_under(target, forbidden) for forbidden in forbidden_targets):
                offenders.append(
                    f"{graph.paths[source].relative_to(PROJECT_ROOT)}: "
                    f"{source} -> {target}"
                )
    return sorted(offenders)


def _is_under(module_name: str, root: str) -> bool:
    return module_name == root or module_name.startswith(f"{root}.")


def _package_edges(graph: ImportGraph) -> list[str]:
    edges: set[str] = set()
    for source, targets in graph.edges.items():
        source_package = ".".join(source.split(".")[:2])
        for target in targets:
            target_package = ".".join(target.split(".")[:2])
            if source_package == target_package:
                continue
            edges.add(f"{source_package} -> {target_package}")
    return sorted(edges)


def _format_cycles(cycles: list[tuple[str, ...]]) -> str:
    if not cycles:
        return ""
    return "\n\n".join("\n".join(cycle) for cycle in cycles)


def _format_edges(edges: list[str]) -> str:
    return "\n".join(edges)
