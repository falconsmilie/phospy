"""Shared static import-graph support for architecture tests.

The analyser intentionally includes imports nested under ``TYPE_CHECKING`` and
function bodies. It resolves absolute and relative imports, plus literal calls
to ``importlib.import_module`` (including directly imported or aliased
bindings) and ``__import__``. Computed module names, aliases introduced through
assignment or re-export, plugin discovery, import hooks, and other
runtime-generated imports are outside this static analyser's scope.
"""

from __future__ import annotations

import ast
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path


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


def build_import_graph(*, package_root: Path, src_root: Path) -> ImportGraph:
    """Build the resolved static import graph for one ``phospy`` source tree."""

    return _build_import_graph(
        package_root=package_root,
        src_root=src_root,
        import_extractor=imported_modules,
    )


def build_eager_import_graph(*, package_root: Path, src_root: Path) -> ImportGraph:
    """Build the complete module-execution import graph for ``phospy``."""

    return _build_import_graph(
        package_root=package_root,
        src_root=src_root,
        import_extractor=eager_imported_modules,
    )


def _build_import_graph(
    *,
    package_root: Path,
    src_root: Path,
    import_extractor: Callable[[str, Path, ast.AST], Iterable[tuple[str, int]]],
) -> ImportGraph:

    paths = {
        module_name(path, src_root=src_root): path
        for path in package_root.rglob("*.py")
    }
    modules = frozenset(paths)
    records: list[ImportRecord] = []
    module_edges: dict[str, frozenset[str]] = {}
    for source_module, path in paths.items():
        resolved_edges: set[str] = set()
        tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
        for imported, line_number in import_extractor(source_module, path, tree):
            if not imported.startswith("phospy"):
                continue
            resolved = resolve_module(imported, modules)
            if resolved is None or resolved == source_module:
                continue
            resolved_edges.add(resolved)
            records.append(
                ImportRecord(
                    source_module=source_module,
                    source_path=path,
                    target=resolved,
                    line=line_number,
                )
            )
        module_edges[source_module] = frozenset(resolved_edges)
    return ImportGraph(
        modules=modules,
        module_paths=paths,
        module_edges=module_edges,
        package_edges=frozenset(package_edges(module_edges)),
        records=tuple(records),
    )


def module_name(path: Path, *, src_root: Path) -> str:
    relative = path.relative_to(src_root).with_suffix("")
    parts = relative.parts
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


def imported_modules(
    source_module: str,
    path: Path,
    tree: ast.AST,
) -> Iterable[tuple[str, int]]:
    """Yield every statically visible import target and source line."""

    importlib_aliases, import_module_aliases = _importlib_bindings(tree)
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                yield alias.name, node.lineno
        elif isinstance(node, ast.ImportFrom):
            imported = resolve_import_from(source_module, path, node)
            if imported is None:
                continue
            yield imported, node.lineno
            for alias in node.names:
                yield (
                    f"{imported}.{alias.name}" if imported else alias.name,
                    node.lineno,
                )
        elif isinstance(node, ast.Call):
            imported = static_dynamic_import_target(
                node,
                importlib_aliases=importlib_aliases,
                import_module_aliases=import_module_aliases,
            )
            if imported is not None:
                yield imported, node.lineno


def eager_imported_modules(
    source_module: str,
    path: Path,
    tree: ast.AST,
) -> Iterable[tuple[str, int]]:
    """Yield imports reachable while a module and its class bodies execute."""

    visitor = _EagerNodeVisitor()
    visitor.visit(tree)
    importlib_aliases, import_module_aliases = _importlib_bindings_from_nodes(
        visitor.nodes
    )
    for node in visitor.nodes:
        if isinstance(node, ast.Import):
            for alias in node.names:
                yield alias.name, node.lineno
        elif isinstance(node, ast.ImportFrom):
            imported = resolve_import_from(source_module, path, node)
            if imported is None:
                continue
            yield imported, node.lineno
            for alias in node.names:
                yield (
                    f"{imported}.{alias.name}" if imported else alias.name,
                    node.lineno,
                )
        elif isinstance(node, ast.Call):
            imported = static_dynamic_import_target(
                node,
                importlib_aliases=importlib_aliases,
                import_module_aliases=import_module_aliases,
            )
            if imported is not None:
                yield imported, node.lineno


class _EagerNodeVisitor(ast.NodeVisitor):
    """Collect nodes evaluated during import, conservatively across branches."""

    def __init__(self) -> None:
        self.nodes: list[ast.AST] = []

    def visit(self, node: ast.AST) -> None:
        self.nodes.append(node)
        super().visit(node)

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        self._visit_function_declaration(node)

    def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:
        self._visit_function_declaration(node)

    def visit_Lambda(self, node: ast.Lambda) -> None:
        self.visit(node.args)

    def visit_If(self, node: ast.If) -> None:
        self.visit(node.test)
        if _is_type_checking_guard(node.test):
            for child in node.orelse:
                self.visit(child)
            return
        for child in (*node.body, *node.orelse):
            self.visit(child)

    def _visit_function_declaration(
        self,
        node: ast.FunctionDef | ast.AsyncFunctionDef,
    ) -> None:
        for decorator in node.decorator_list:
            self.visit(decorator)
        self.visit(node.args)
        if node.returns is not None:
            self.visit(node.returns)


def _is_type_checking_guard(node: ast.AST) -> bool:
    return (isinstance(node, ast.Name) and node.id == "TYPE_CHECKING") or (
        isinstance(node, ast.Attribute)
        and isinstance(node.value, ast.Name)
        and node.value.id == "typing"
        and node.attr == "TYPE_CHECKING"
    )


def _importlib_bindings(tree: ast.AST) -> tuple[frozenset[str], frozenset[str]]:
    """Return module and callable names bound by explicit importlib imports."""

    return _importlib_bindings_from_nodes(ast.walk(tree))


def _importlib_bindings_from_nodes(
    nodes: Iterable[ast.AST],
) -> tuple[frozenset[str], frozenset[str]]:
    """Return importlib bindings present in an already selected node set."""

    module_aliases: set[str] = set()
    import_module_aliases: set[str] = set()
    for node in nodes:
        if isinstance(node, ast.Import):
            module_aliases.update(
                alias.asname or alias.name
                for alias in node.names
                if alias.name == "importlib"
            )
        elif (
            isinstance(node, ast.ImportFrom)
            and node.level == 0
            and node.module == "importlib"
        ):
            import_module_aliases.update(
                alias.asname or alias.name
                for alias in node.names
                if alias.name == "import_module"
            )
    return frozenset(module_aliases), frozenset(import_module_aliases)


def resolve_import_from(
    source_module: str,
    path: Path,
    node: ast.ImportFrom,
) -> str | None:
    if node.level == 0:
        return node.module
    parts = source_module.split(".")
    package_parts = parts if path.name == "__init__.py" else parts[:-1]
    if node.level > 1:
        package_parts = package_parts[: -(node.level - 1)]
    module_parts = [] if node.module is None else node.module.split(".")
    return ".".join((*package_parts, *module_parts))


def resolve_module(imported: str, modules: frozenset[str]) -> str | None:
    parts = imported.split(".")
    for index in range(len(parts), 0, -1):
        candidate = ".".join(parts[:index])
        if candidate in modules:
            return candidate
    return None


def static_dynamic_import_target(
    node: ast.Call,
    *,
    importlib_aliases: frozenset[str] = frozenset({"importlib"}),
    import_module_aliases: frozenset[str] = frozenset(),
) -> str | None:
    is_importlib_call = (
        isinstance(node.func, ast.Attribute)
        and node.func.attr == "import_module"
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id in importlib_aliases
    )
    is_importlib_binding_call = (
        isinstance(node.func, ast.Name) and node.func.id in import_module_aliases
    )
    is_dunder_import = isinstance(node.func, ast.Name) and node.func.id == "__import__"
    if (
        not (is_importlib_call or is_importlib_binding_call or is_dunder_import)
        or not node.args
    ):
        return None
    first_arg = node.args[0]
    if not isinstance(first_arg, ast.Constant) or not isinstance(first_arg.value, str):
        return None
    return first_arg.value


def package_edges(
    module_edges: Mapping[str, frozenset[str]],
) -> set[tuple[str, str]]:
    edges: set[tuple[str, str]] = set()
    for source, targets in module_edges.items():
        source_package = package_name(source)
        if source_package is None:
            continue
        for target in targets:
            target_package = package_name(target)
            if target_package is None or target_package == source_package:
                continue
            edges.add((source_package, target_package))
    return edges


def package_name(name: str) -> str | None:
    parts = name.split(".")
    if len(parts) < 2 or parts[0] != "phospy":
        return None
    return ".".join(parts[:2])
