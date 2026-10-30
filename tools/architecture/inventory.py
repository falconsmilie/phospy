"""Render the maintained architecture inventory."""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

from tools.architecture.final_policy import (
    EXACT_OWNER_EDGE_RULES,
    PERMANENT_LAZY_LOADERS,
    PERMANENT_SPECIAL_MODULES,
    TOP_LEVEL_ROOT_CLASSIFICATIONS,
)
from tools.architecture.import_graph import build_import_graph
from tools.architecture.logical_ownership import (
    load_permanent_condensation_classifications,
    owner_graph_views,
)


def render_inventory(repo_root: Path) -> str:
    """Render the complete tracked production and ownership inventory."""

    root = repo_root.resolve()
    tracked = _git_paths(root, "src/phospy")
    python_paths = tuple(path for path in tracked if path.endswith(".py"))
    owners = _owner_rows(root)
    owners_by_module = {row["module"]: row["owner"] for row in owners}
    observed_graph = build_import_graph(
        package_root=root / "src/phospy",
        src_root=root / "src",
    )
    permanent_classifications = load_permanent_condensation_classifications(
        root / "architecture/permanent-condensation-classifications.json"
    )
    owner_views = owner_graph_views(
        observed_graph.module_edges,
        owners_by_module=owners_by_module,
        classifications=permanent_classifications,
    )
    classifications = json.loads(
        (root / "architecture/final-state-classifications.json").read_text(
            encoding="utf-8"
        )
    )
    forwarding_classifications = classifications["structural_forwarding_modules"]

    lines = [
        "# Architecture inventory",
        "",
        "This generated inventory records tracked production files, module ownership,",
        "and the complete and policy-filtered dependency views. Executable architecture",
        "tests compare it with tracked source, the one-row-per-module owner map, the",
        "exact owner-edge policy, and the permanent classification registry.",
        "",
        "## Counts and architecture assertions",
        "",
        f"- Tracked production files: **{len(tracked)}**.",
        f"- Tracked Python modules: **{len(python_paths)}**.",
        f"- Owner-map entries: **{len(owners)}**.",
        "- Temporary architecture exceptions: **0**.",
        "- Retired tracked package roots (`contracts`, `errors`, `frames`, `policies`, `science`, `tables`, `validation`, `workflows`): **0**.",
        "- Unclassified structural forwarding candidates, old-path literals, compatibility behavior, and lazy loading: **0**.",
        "- The complete eager module graph and policy-filtered owner graph are required to be acyclic by maintained tests.",
        "",
        "## Top-level roots",
        "",
        "| Root | Current classification |",
        "| --- | --- |",
    ]
    lines.extend(
        f"| `{root_name}` | `{classification}` |"
        for root_name, classification in sorted(TOP_LEVEL_ROOT_CLASSIFICATIONS.items())
    )

    lines.extend(
        [
            "",
            "## Complete statically observed owner graph (before classifications)",
            "",
            f"This view contains **{len(owner_views.complete.edges)}** owner edges before",
            "accepted condensation classifications are applied. Every edge remains",
            "subject to the maintained exact-rule comparison; a classified edge is an",
            "auditable accepted dependency, not an unauthorised import.",
            "",
            "### Strongly connected components",
            "",
        ]
    )
    if owner_views.complete.cyclic_components:
        lines.extend(
            f"- `{', '.join(component)}`"
            for component in owner_views.complete.cyclic_components
        )
    else:
        lines.append("- No multi-owner strongly connected components.")
    singleton_count = sum(
        len(component) == 1
        for component in owner_views.complete.strongly_connected_components
    )
    lines.extend(
        [
            f"- Single-owner components: **{singleton_count}**.",
            "",
            "### Complete observed edges",
            "",
            "| Source | Target | Classification | Authority |",
            "| --- | --- | --- | --- |",
        ]
    )
    lines.extend(
        f"| `{source}` | `{target}` | `{rule.kind}` | {rule.authority} |"
        for source, target in sorted(owner_views.complete.edges)
        for rule in (EXACT_OWNER_EDGE_RULES[source, target],)
    )

    lines.extend(
        [
            "",
            "## Excluded imports from accepted classifications",
            "",
            "These exact module imports remain visible in the complete view and are",
            "removed only before policy condensation and the DAG gate.",
            "",
            "| Classification ID | Source path | Source owner | Target owner | Exact import | Rationale |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
    )
    lines.extend(
        f"| `{item.classification.classification_id}` | "
        f"`{item.classification.source_path}` | `{item.source_owner}` | "
        f"`{item.target_owner}` | `{item.classification.source_module} -> "
        f"{item.classification.target_module}.{item.classification.symbol}` | "
        f"{_escape(item.classification.rationale)} |"
        for item in owner_views.excluded_imports
    )

    lines.extend(
        [
            "",
            "## Policy-filtered owner graph used by the DAG gate",
            "",
            f"This view contains **{len(owner_views.policy_filtered.edges)}** owner edges",
            "after the exact accepted imports above are removed. It is the graph checked",
            "for owner-level acyclicity.",
            "",
            "### Strongly connected components",
            "",
        ]
    )
    if owner_views.policy_filtered.cyclic_components:
        lines.extend(
            f"- `{', '.join(component)}`"
            for component in owner_views.policy_filtered.cyclic_components
        )
    else:
        lines.append(
            "- No multi-owner strongly connected components; the view is acyclic."
        )
    filtered_singleton_count = sum(
        len(component) == 1
        for component in owner_views.policy_filtered.strongly_connected_components
    )
    lines.extend(
        [
            f"- Single-owner components: **{filtered_singleton_count}**.",
            "",
            "### Policy-filtered edges",
            "",
            "| Source | Target | Classification | Authority |",
            "| --- | --- | --- | --- |",
        ]
    )
    lines.extend(
        f"| `{source}` | `{target}` | `{rule.kind}` | {rule.authority} |"
        for source, target in sorted(owner_views.policy_filtered.edges)
        for rule in (EXACT_OWNER_EDGE_RULES[source, target],)
    )

    lines.extend(
        [
            "",
            "## Static-analysis scope",
            "",
            "The analyser includes absolute and relative imports, imports nested under",
            "`TYPE_CHECKING` or function bodies, and literal `importlib.import_module` or",
            "`__import__` calls. It does not resolve computed module names, bindings created",
            "through assignment or re-export, plugin discovery, import hooks, or other",
            "runtime-generated imports. The complete owner view therefore reports static",
            "dependency potential; the separately enforced eager-module graph describes",
            "imports reachable while modules and class bodies execute.",
        ]
    )

    lines.extend(
        [
            "",
            "## Public facades, adapters, and foundations",
            "",
            "The root table is authoritative for every top-level owner. The public",
            "supported facades are `phospy.api` and `phospy.api.advanced`; the",
            "`phospy` root provides metadata only, and `phospy.io`",
            "is the inward adapter. Shared scientific, record, composition, and",
            "technical foundations are the roots carrying those classifications above.",
            "",
            "Exact permanent facade/safety modules:",
            "",
        ]
    )
    lines.extend(
        f"- `{path}`: {explanation}"
        for path, explanation in sorted(PERMANENT_SPECIAL_MODULES.items())
    )
    lines.extend(["", "Exact permanent lazy loaders:", ""])
    lines.extend(
        f"- `{path}`: {explanation}"
        for path, explanation in sorted(PERMANENT_LAZY_LOADERS.items())
    )

    lines.extend(
        [
            "",
            "## Structurally detected forwarding candidates",
            "",
            "The AST detector identifies modules whose executable top-level body is",
            "limited to imports, `__all__` maintenance, and direct aliases of imported",
            "names. A docstring and `__future__` imports do not affect detection. The",
            "detector establishes forwarding-shaped structure only; each candidate",
            "still requires the explicit policy disposition below.",
            "",
            "The declared scope does not interpret computed exports, forwarding inside",
            "functions or classes, module `__getattr__`, import hooks, or arbitrary",
            "runtime binding, and it does not attempt to model arbitrary Python behavior.",
            "",
            "| Source path | Policy disposition | Rationale |",
            "| --- | --- | --- |",
        ]
    )
    lines.extend(
        f"| `{path}` | `{item['disposition']}` | {_escape(item['rationale'])} |"
        for path, item in sorted(forwarding_classifications.items())
    )

    lines.extend(
        [
            "",
            "## Classified versioned and historical behavior",
            "",
            "These records distinguish legitimate versioned schema or legacy-data",
            "identifiers from live compatibility routes. Their exact sets are enforced.",
            "",
            "### Retired-route string literals",
            "",
        ]
    )
    for path, values in sorted(classifications["retired_route_literals"].items()):
        lines.append(f"- `{path}`")
        lines.extend(f"  - `{_escape(value)}`" for value in values)

    lines.extend(["", "### Legacy behavior markers", ""])
    for behavior, paths in sorted(classifications["legacy_behavior_modules"].items()):
        lines.append(f"- `{behavior}`")
        lines.extend(f"  - `{path}`" for path in paths)

    lines.extend(
        [
            "",
            "## Complete tracked production file inventory",
            "",
        ]
    )
    lines.extend(f"- `{path}`" for path in tracked)

    lines.extend(
        [
            "",
            "## Complete Python module ownership inventory",
            "",
            "| Source path | Module | Owner | Rationale |",
            "| --- | --- | --- | --- |",
        ]
    )
    lines.extend(
        "| `{source_path}` | `{module}` | `{owner}` | {rationale} |".format(
            **{key: _escape(value) for key, value in row.items()}
        )
        for row in owners
    )
    return "\n".join(lines) + "\n"


def _escape(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ")


def _git_paths(repo_root: Path, pathspec: str) -> tuple[str, ...]:
    result = subprocess.run(
        ["git", "ls-files", "-z", "--", pathspec],
        cwd=repo_root,
        check=True,
        capture_output=True,
    )
    paths = sorted(item.decode("utf-8") for item in result.stdout.split(b"\0") if item)
    return tuple(path for path in paths if (repo_root / path).is_file())


def _owner_rows(repo_root: Path) -> tuple[dict[str, str], ...]:
    path = repo_root / "architecture/module-owners.csv"
    with path.open(encoding="utf-8", newline="") as handle:
        return tuple(sorted(csv.DictReader(handle), key=lambda row: row["source_path"]))
