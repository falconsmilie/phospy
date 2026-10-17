"""Generate the maintained architecture inventory."""

from __future__ import annotations

import argparse
import csv
import importlib
import json
import subprocess
import sys
from pathlib import Path

OUTPUT = Path("docs/architecture_inventory.md")


def render_inventory(repo_root: Path) -> str:
    """Render the complete tracked production and ownership inventory."""

    root = repo_root.resolve()
    architecture_tests = root / "tests" / "architecture"
    sys.path.insert(0, str(architecture_tests))
    try:
        policy = importlib.import_module("final_architecture_policy")
    finally:
        sys.path.remove(str(architecture_tests))

    exact_owner_edge_rules = policy.EXACT_OWNER_EDGE_RULES
    permanent_lazy_loaders = policy.PERMANENT_LAZY_LOADERS
    permanent_special_modules = policy.PERMANENT_SPECIAL_MODULES
    top_level_root_classifications = policy.TOP_LEVEL_ROOT_CLASSIFICATIONS

    tracked = _git_paths(root, "src/phospy")
    python_paths = tuple(path for path in tracked if path.endswith(".py"))
    owners = _owner_rows(root)
    classifications = json.loads(
        (root / "architecture/final-state-classifications.json").read_text(
            encoding="utf-8"
        )
    )

    lines = [
        "# Architecture inventory",
        "",
        "This generated inventory records the current production tree and its",
        "The executable architecture tests regenerate and compare these facts",
        "architecture authority. Maintained tests compare it against tracked source,",
        "the one-row-per-module owner map, the complete",
        "observed import graph, and the permanent classification registry.",
        "",
        "## Counts and architecture assertions",
        "",
        f"- Tracked production files: **{len(tracked)}**.",
        f"- Tracked Python modules: **{len(python_paths)}**.",
        f"- Owner-map entries: **{len(owners)}**.",
        "- Temporary architecture exceptions: **0**.",
        "- Retired tracked package roots (`contracts`, `science`, `tables`, `workflows`): **0**.",
        "- Unclassified facade forwarding, old-path literals, compatibility behavior, and lazy loading: **0**.",
        "- The complete eager module graph and the owner condensation graph are required to be acyclic by maintained tests.",
        "",
        "## Top-level roots",
        "",
        "| Root | Final classification |",
        "| --- | --- |",
    ]
    lines.extend(
        f"| `{root_name}` | `{classification}` |"
        for root_name, classification in sorted(top_level_root_classifications.items())
    )

    lines.extend(
        [
            "",
            "## Observed owner edges and exact authorities",
            "",
            "Every row is both observed and explicitly allowed. A maintained exact-set",
            "comparison rejects added edges and unused broad permissions.",
            "",
            "| Source | Target | Classification | Authority |",
            "| --- | --- | --- | --- |",
        ]
    )
    lines.extend(
        f"| `{source}` | `{target}` | `{rule.kind}` | {rule.authority} |"
        for (source, target), rule in sorted(exact_owner_edge_rules.items())
    )

    lines.extend(
        [
            "",
            "## Public facades, adapters, and foundations",
            "",
            "The root table is authoritative for every top-level owner. The public",
            "facades are `phospy`, `phospy.api`, and `phospy.advanced`; `phospy.io`",
            "is the inward adapter. Shared scientific, record, composition, and",
            "technical foundations are the roots carrying those classifications above.",
            "",
            "Exact permanent facade/safety modules:",
            "",
        ]
    )
    lines.extend(
        f"- `{path}`: {explanation}"
        for path, explanation in sorted(permanent_special_modules.items())
    )
    lines.extend(["", "Exact permanent lazy loaders:", ""])
    lines.extend(
        f"- `{path}`: {explanation}"
        for path, explanation in sorted(permanent_lazy_loaders.items())
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
    return tuple(
        sorted(item.decode("utf-8") for item in result.stdout.split(b"\0") if item)
    )


def _owner_rows(repo_root: Path) -> tuple[dict[str, str], ...]:
    path = repo_root / "architecture/module-owners.csv"
    with path.open(encoding="utf-8", newline="") as handle:
        return tuple(sorted(csv.DictReader(handle), key=lambda row: row["source_path"]))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root", type=Path, default=Path(__file__).resolve().parents[2]
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    output = args.repo_root.resolve() / OUTPUT
    rendered = render_inventory(args.repo_root)
    if args.check:
        if not output.is_file() or output.read_text(encoding="utf-8") != rendered:
            print(f"architecture inventory is stale: {output}", file=sys.stderr)
            return 1
        print(f"architecture inventory is current: {output}")
        return 0
    output.write_text(rendered, encoding="utf-8", newline="\n")
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
