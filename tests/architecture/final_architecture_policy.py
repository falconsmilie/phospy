"""Exact permanent architecture policy."""

from __future__ import annotations

import ast
from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import PurePosixPath


@dataclass(frozen=True, slots=True)
class OwnerEdgeRule:
    """One exact observed owner edge and its governing decision."""

    kind: str
    authority: str


def is_structural_forwarding_candidate(
    source: str,
    *,
    filename: str = "<unknown>",
) -> bool:
    """Identify a deliberately narrow class of forwarding-shaped modules.

    The detector ignores a module docstring and ``__future__`` imports, then
    accepts a runtime body made only from imports, ``__all__`` maintenance, and
    direct name aliases of imported bindings. At least one absolute ``phospy``
    or relative import must be present. It identifies structure only; policy
    separately decides whether a candidate is a supported facade, a deliberate
    internal facade, or an obsolete compatibility route.

    This is not an exhaustive Python re-export analyser. It intentionally does
    not interpret computed exports, forwarding hidden inside functions or
    classes, module ``__getattr__``, import hooks, or arbitrary runtime binding.
    """

    tree = ast.parse(source, filename=filename)
    body = list(tree.body)
    if body and _is_docstring_statement(body[0]):
        body = body[1:]

    imported_bindings: set[str] = set()
    alias_bindings: set[str] = set()
    saw_local_forwarding_import = False
    saw_runtime_statement = False

    for node in body:
        if isinstance(node, ast.ImportFrom) and node.module == "__future__":
            continue
        saw_runtime_statement = True
        if isinstance(node, ast.ImportFrom):
            module = node.module or ""
            is_local = (
                node.level > 0 or module == "phospy" or module.startswith("phospy.")
            )
            saw_local_forwarding_import |= is_local
            imported_bindings.update(
                alias.asname or alias.name for alias in node.names if alias.name != "*"
            )
            continue
        if isinstance(node, ast.Import):
            saw_local_forwarding_import |= any(
                alias.name == "phospy" or alias.name.startswith("phospy.")
                for alias in node.names
            )
            imported_bindings.update(
                alias.asname or alias.name.split(".", maxsplit=1)[0]
                for alias in node.names
            )
            continue
        if _is_dunder_all_assignment(node):
            continue
        if _is_direct_imported_name_alias(
            node,
            available_bindings=imported_bindings | alias_bindings,
            alias_bindings=alias_bindings,
        ):
            continue
        return False

    return saw_runtime_statement and saw_local_forwarding_import


def _is_docstring_statement(node: ast.stmt) -> bool:
    return (
        isinstance(node, ast.Expr)
        and isinstance(node.value, ast.Constant)
        and isinstance(node.value.value, str)
    )


def _is_dunder_all_assignment(node: ast.stmt) -> bool:
    if isinstance(node, ast.Assign):
        return all(
            isinstance(target, ast.Name) and target.id == "__all__"
            for target in node.targets
        )
    if isinstance(node, ast.AnnAssign):
        return isinstance(node.target, ast.Name) and node.target.id == "__all__"
    if isinstance(node, ast.AugAssign):
        return isinstance(node.target, ast.Name) and node.target.id == "__all__"
    return False


def _is_direct_imported_name_alias(
    node: ast.stmt,
    *,
    available_bindings: set[str],
    alias_bindings: set[str],
) -> bool:
    if isinstance(node, ast.Assign):
        if (
            not isinstance(node.value, ast.Name)
            or node.value.id not in available_bindings
        ):
            return False
        targets = node.targets
    elif isinstance(node, ast.AnnAssign):
        if (
            node.value is None
            or not isinstance(node.value, ast.Name)
            or node.value.id not in available_bindings
        ):
            return False
        targets = (node.target,)
    else:
        return False
    if not all(isinstance(target, ast.Name) for target in targets):
        return False
    alias_bindings.update(
        target.id for target in targets if isinstance(target, ast.Name)
    )
    return True


TOP_LEVEL_ROOT_CLASSIFICATIONS: Mapping[str, str] = {
    ".": "public_facade_machinery",
    "activities": "passive_shared_scientific_foundation",
    "advanced": "public_facade",
    "api": "public_facade",
    "data": "packaged_resource_payload",
    "datasets": "capability",
    "differential": "capability",
    "enrichment": "capability",
    "errors": "shared_technical_foundation",
    "evidence": "shared_scientific_foundation",
    "frames": "shared_technical_foundation",
    "io": "adapter",
    "kinase": "capability",
    "phosphosites": "shared_scientific_foundation",
    "policies": "shared_technical_foundation",
    "preparation_records": "passive_shared_record_foundation",
    "preprocessing": "capability",
    "provenance": "shared_scientific_foundation",
    "quantitative": "shared_scientific_foundation",
    "references": "shared_scientific_foundation",
    "run_lineage": "composition_foundation",
    "signalome": "capability",
    "statistics": "shared_scientific_foundation",
    "taxonomy": "shared_scientific_foundation",
    "validation": "shared_technical_foundation",
}

RETIRED_TOP_LEVEL_ROOTS = frozenset({"contracts", "science", "tables", "workflows"})

PERMANENT_SPECIAL_MODULES: Mapping[str, str] = {
    "src/phospy/_api_inventory.py": (
        "authoritative aggregate inventory data; explicitly free of compatibility routes"
    ),
    "src/phospy/differential/aggregation/experimental.py": (
        "permanent fail-closed scientific safety seam governed by ADR-0041; it is "
        "not a supported compatibility import route"
    ),
}

PERMANENT_LAZY_LOADERS: Mapping[str, str] = {
    "src/phospy/__init__.py": "three-name supported root facade",
    "src/phospy/io/__init__.py": "IO adapter composition",
    "src/phospy/io/bundles/__init__.py": "IO bundle adapter composition",
    "src/phospy/io/readers/__init__.py": "IO reader adapter composition",
    "src/phospy/preprocessing/contracts/plan.py": (
        "cycle-free owner-internal plan/interpreter composition"
    ),
}

_OBSERVED_TARGETS: Mapping[str, frozenset[str]] = {
    "activities": frozenset({"errors", "frames", "quantitative"}),
    "datasets": frozenset(
        {
            "errors",
            "evidence",
            "frames",
            "io",
            "phosphosites",
            "preparation_records",
            "preprocessing",
            "provenance",
            "quantitative",
            "references",
            "run_lineage",
            "taxonomy",
            "validation",
        }
    ),
    "differential": frozenset(
        {
            "datasets",
            "errors",
            "evidence",
            "frames",
            "phosphosites",
            "policies",
            "preparation_records",
            "provenance",
            "quantitative",
            "references",
            "run_lineage",
            "statistics",
        }
    ),
    "enrichment": frozenset(
        {
            "errors",
            "frames",
            "provenance",
            "quantitative",
            "run_lineage",
            "statistics",
        }
    ),
    "evidence": frozenset(
        {"errors", "frames", "phosphosites", "policies", "quantitative"}
    ),
    "frames": frozenset({"errors"}),
    "io": frozenset(
        {
            "activities",
            "datasets",
            "enrichment",
            "errors",
            "evidence",
            "frames",
            "kinase",
            "phosphosites",
            "preparation_records",
            "preprocessing",
            "provenance",
            "quantitative",
            "references",
            "run_lineage",
            "signalome",
            "taxonomy",
            "validation",
        }
    ),
    "kinase": frozenset(
        {
            "activities",
            "datasets",
            "errors",
            "evidence",
            "frames",
            "phosphosites",
            "policies",
            "provenance",
            "quantitative",
            "references",
            "run_lineage",
            "statistics",
            "taxonomy",
        }
    ),
    "phosphosites": frozenset({"errors", "frames", "taxonomy"}),
    "policies": frozenset({"errors"}),
    "preparation_records": frozenset({"errors", "frames", "policies", "provenance"}),
    "preprocessing": frozenset(
        {
            "errors",
            "evidence",
            "frames",
            "phosphosites",
            "policies",
            "preparation_records",
            "provenance",
            "quantitative",
            "references",
            "taxonomy",
        }
    ),
    "provenance": frozenset({"errors"}),
    "public_facades": frozenset(
        {
            "datasets",
            "differential",
            "enrichment",
            "errors",
            "evidence",
            "io",
            "kinase",
            "preprocessing",
            "provenance",
            "quantitative",
            "references",
            "signalome",
            "statistics",
            "taxonomy",
        }
    ),
    "quantitative": frozenset({"errors", "provenance"}),
    "references": frozenset(
        {"errors", "frames", "phosphosites", "policies", "provenance", "taxonomy"}
    ),
    "run_lineage": frozenset(
        {"errors", "evidence", "preprocessing", "provenance", "references"}
    ),
    "signalome": frozenset(
        {
            "datasets",
            "errors",
            "evidence",
            "frames",
            "kinase",
            "phosphosites",
            "policies",
            "provenance",
            "quantitative",
            "references",
            "run_lineage",
            "validation",
        }
    ),
    "validation": frozenset({"errors"}),
}

_CAPABILITY_SOURCES = frozenset(
    {"datasets", "differential", "enrichment", "kinase", "preprocessing", "signalome"}
)
_SHARED_SOURCES = frozenset(
    {
        "activities",
        "evidence",
        "frames",
        "phosphosites",
        "policies",
        "preparation_records",
        "provenance",
        "quantitative",
        "references",
        "validation",
    }
)


def _edge_rule(source: str, target: str) -> OwnerEdgeRule:
    if (source, target) == ("datasets", "io"):
        return OwnerEdgeRule(
            "permanent_classification",
            "ADR-0053 D05B-1 dataset builder path-reader authority",
        )
    if (source, target) == ("datasets", "validation"):
        return OwnerEdgeRule(
            "permanent_classification",
            "ADR-0053 D05B-2 dataset input-source validation authority",
        )
    if source == "public_facades":
        return OwnerEdgeRule("public_facade", "ADR-0054 aggregate facade")
    if source == "io":
        return OwnerEdgeRule("adapter", "ADR-0053 inward IO adapter")
    if source == "run_lineage":
        return OwnerEdgeRule("composition", "ADR-0053 run-lineage composition row")
    if source in _CAPABILITY_SOURCES:
        return OwnerEdgeRule("capability", f"ADR-0053 {source} direction")
    if source in _SHARED_SOURCES:
        authority = (
            "ACT-10 exact activities row"
            if source == "activities"
            else f"ADR-0053 {source} row"
        )
        return OwnerEdgeRule("shared", authority)
    raise AssertionError(f"unclassified owner-edge source: {source!r}")


EXACT_OWNER_EDGE_RULES: Mapping[tuple[str, str], OwnerEdgeRule] = {
    (source, target): _edge_rule(source, target)
    for source, targets in _OBSERVED_TARGETS.items()
    for target in targets
}


def owner_edge_policy_problems(
    *,
    observed_edges: Iterable[tuple[str, str]],
    policy_filtered_edges: Iterable[tuple[str, str]],
    rules: Mapping[tuple[str, str], OwnerEdgeRule] = EXACT_OWNER_EDGE_RULES,
) -> list[str]:
    """Compare complete owner edges and check the policy-filtered DAG."""

    observed = frozenset(observed_edges)
    permitted = frozenset(rules)
    problems = [
        f"observed owner edge has no exact rule: {source} -> {target}"
        for source, target in sorted(observed - permitted)
    ]
    problems.extend(
        f"unused or missing exact owner-edge permission: {source} -> {target}"
        for source, target in sorted(permitted - observed)
    )
    for (source, target), rule in sorted(rules.items()):
        expected = _edge_rule(source, target)
        if rule != expected:
            problems.append(
                f"incorrect owner-edge classification: {source} -> {target}; "
                f"expected {expected!r}, got {rule!r}"
            )
    activities_targets = {
        target for source, target in observed if source == "activities"
    }
    expected_activities = {"errors", "frames", "quantitative"}
    if activities_targets != expected_activities:
        problems.append(
            "activities outbound edges differ from ACT-10: "
            f"expected {sorted(expected_activities)!r}, got {sorted(activities_targets)!r}"
        )
    problems.extend(_cycle_problems(frozenset(policy_filtered_edges)))
    return problems


def tracked_source_inventory_problems(
    *,
    tracked_paths: Iterable[str],
    owner_paths: Iterable[str],
    owner_modules: Iterable[str],
) -> list[str]:
    """Validate complete tracked production-file and Python-owner accounting."""

    paths = tuple(tracked_paths)
    owners = tuple(owner_paths)
    modules = tuple(owner_modules)
    problems: list[str] = []
    path_counts = Counter(paths)
    owner_counts = Counter(owners)
    module_counts = Counter(modules)
    problems.extend(
        f"duplicate tracked production path: {path}"
        for path, count in sorted(path_counts.items())
        if count != 1
    )
    problems.extend(
        f"duplicate owner-map source path: {path}"
        for path, count in sorted(owner_counts.items())
        if count != 1
    )
    problems.extend(
        f"duplicate owner-map module definition: {module}"
        for module, count in sorted(module_counts.items())
        if count != 1
    )
    if len(owners) != len(modules):
        problems.append(
            "owner-map source-path and module counts differ: "
            f"{len(owners)} paths != {len(modules)} modules"
        )
    python_paths = {path for path in paths if path.endswith(".py")}
    owner_path_set = set(owners)
    problems.extend(
        f"tracked Python module lacks exactly one owner: {path}"
        for path in sorted(python_paths - owner_path_set)
    )
    problems.extend(
        f"owner-map path is not a tracked Python module: {path}"
        for path in sorted(owner_path_set - python_paths)
    )
    observed_roots = {_top_level_root(path) for path in paths}
    problems.extend(
        f"unclassified top-level production root: {root}"
        for root in sorted(observed_roots - set(TOP_LEVEL_ROOT_CLASSIFICATIONS))
    )
    problems.extend(
        f"stale top-level root classification: {root}"
        for root in sorted(set(TOP_LEVEL_ROOT_CLASSIFICATIONS) - observed_roots)
    )
    for retired in sorted(RETIRED_TOP_LEVEL_ROOTS & observed_roots):
        problems.append(f"retired top-level production root is tracked: {retired}")
    non_python = {path for path in paths if not path.endswith(".py")}
    problems.extend(
        f"unclassified non-Python production resource: {path}"
        for path in sorted(non_python)
        if not path.startswith("src/phospy/data/")
    )
    return problems


def old_path_literal_inventory(
    source_by_path: Mapping[str, str],
) -> dict[str, tuple[str, ...]]:
    """Return exact retired-route string literals without treating them as imports."""

    result: dict[str, tuple[str, ...]] = {}
    for path, source in sorted(source_by_path.items()):
        tree = ast.parse(source, filename=path)
        values = sorted(
            {
                node.value
                for node in ast.walk(tree)
                if isinstance(node, ast.Constant)
                and isinstance(node.value, str)
                and _contains_retired_route(node.value)
            }
        )
        if values:
            result[path] = tuple(values)
    return result


def _contains_retired_route(value: str) -> bool:
    return any(f"phospy.{root}" in value for root in RETIRED_TOP_LEVEL_ROOTS)


def _top_level_root(path: str) -> str:
    normalized = PurePosixPath(path)
    parts = normalized.parts
    if parts[:2] != ("src", "phospy"):
        raise ValueError(f"production path is outside src/phospy: {path!r}")
    return "." if len(parts) == 3 else parts[2]


def _cycle_problems(edges: Iterable[tuple[str, str]]) -> list[str]:
    adjacency: dict[str, set[str]] = {}
    for source, target in edges:
        adjacency.setdefault(source, set()).add(target)
        adjacency.setdefault(target, set())
    visiting: set[str] = set()
    visited: set[str] = set()
    stack: list[str] = []
    cycles: set[tuple[str, ...]] = set()

    def visit(owner: str) -> None:
        if owner in visited:
            return
        if owner in visiting:
            start = stack.index(owner)
            cycles.add(tuple(stack[start:] + [owner]))
            return
        visiting.add(owner)
        stack.append(owner)
        for target in sorted(adjacency[owner]):
            visit(target)
        stack.pop()
        visiting.remove(owner)
        visited.add(owner)

    for owner in sorted(adjacency):
        visit(owner)
    return [f"owner dependency cycle: {' -> '.join(cycle)}" for cycle in sorted(cycles)]
