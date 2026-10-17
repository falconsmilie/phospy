from __future__ import annotations

import ast
import csv
import json
import re
import subprocess
import tomllib
from pathlib import Path

from final_architecture_policy import (
    EXACT_OWNER_EDGE_RULES,
    PERMANENT_LAZY_LOADERS,
    PERMANENT_SPECIAL_MODULES,
    OwnerEdgeRule,
    old_path_literal_inventory,
    owner_edge_policy_problems,
    tracked_source_inventory_problems,
)
from import_graph_support import build_import_graph
from logical_ownership_policy import load_module_owners, logical_owner_edges

from tools.architecture.generate_architecture_inventory import render_inventory

ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = ROOT / "src"
PACKAGE_ROOT = SRC_ROOT / "phospy"
OWNER_MAP = ROOT / "architecture/module-owners.csv"
CLASSIFICATIONS = ROOT / "architecture/final-state-classifications.json"
FINAL_INVENTORY = ROOT / "docs/architecture_inventory.md"
PERMANENT_CLASSIFICATIONS = (
    ROOT / "architecture/permanent-condensation-classifications.json"
)


def _git_paths(pathspec: str) -> tuple[str, ...]:
    result = subprocess.run(
        ["git", "ls-files", "-z", "--", pathspec],
        cwd=ROOT,
        check=True,
        capture_output=True,
    )
    return tuple(item.decode("utf-8") for item in result.stdout.split(b"\0") if item)


def _live_owner_edges() -> frozenset[tuple[str, str]]:
    entries = load_module_owners(OWNER_MAP)
    graph = build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    return frozenset(
        logical_owner_edges(
            graph.module_edges,
            owners_by_module={entry.module: entry.owner for entry in entries},
        )
    )


def test_complete_observed_owner_graph_matches_exact_adr_0053_rules() -> None:
    assert owner_edge_policy_problems(observed_edges=_live_owner_edges()) == []


def test_owner_edge_verifier_rejects_added_missing_and_broadened_edges() -> None:
    observed = _live_owner_edges()
    added = set(observed) | {("activities", "kinase")}
    missing = set(observed) - {("activities", "frames")}
    broadened = dict(EXACT_OWNER_EDGE_RULES)
    broadened[("activities", "statistics")] = OwnerEdgeRule(
        "shared", "invented broad permission"
    )

    assert any(
        "has no exact rule: activities -> kinase" in problem
        for problem in owner_edge_policy_problems(observed_edges=added)
    )
    assert any(
        "unused or missing exact owner-edge permission: activities -> frames" in problem
        for problem in owner_edge_policy_problems(observed_edges=missing)
    )
    assert any(
        "unused or missing exact owner-edge permission: activities -> statistics"
        in problem
        for problem in owner_edge_policy_problems(
            observed_edges=observed,
            rules=broadened,
        )
    )


def test_owner_edge_verifier_rejects_reversed_cyclic_and_misclassified_edges() -> None:
    observed = _live_owner_edges()
    reversed_edges = (set(observed) - {("frames", "errors")}) | {("errors", "frames")}
    cyclic_edges = set(observed) | {("errors", "frames")}
    misclassified = dict(EXACT_OWNER_EDGE_RULES)
    misclassified[("activities", "errors")] = OwnerEdgeRule(
        "capability", "incorrect fixture"
    )

    assert any(
        "has no exact rule: errors -> frames" in problem
        for problem in owner_edge_policy_problems(observed_edges=reversed_edges)
    )
    assert any(
        "owner dependency cycle" in problem
        for problem in owner_edge_policy_problems(observed_edges=cyclic_edges)
    )
    assert any(
        "incorrect owner-edge classification: activities -> errors" in problem
        for problem in owner_edge_policy_problems(
            observed_edges=observed,
            rules=misclassified,
        )
    )


def test_final_tracked_source_inventory_is_complete_and_classified() -> None:
    tracked = _git_paths("src/phospy")
    owners = load_module_owners(OWNER_MAP)

    assert len(tracked) == 647
    assert sum(path.endswith(".py") for path in tracked) == 641
    assert len(owners) == 641
    assert (
        tracked_source_inventory_problems(
            tracked_paths=tracked,
            owner_paths=(entry.source_path for entry in owners),
            owner_modules=(entry.module for entry in owners),
        )
        == []
    )


def test_retired_package_roots_are_absent_from_git_and_pyright() -> None:
    tracked = _git_paths("src/phospy")
    retired = {"contracts", "science", "tables", "workflows"}
    assert (
        not {path.split("/")[2] for path in tracked if len(path.split("/")) > 3}
        & retired
    )

    configuration = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    includes = configuration["tool"]["pyright"]["include"]
    assert "src/phospy/tables" not in includes
    assert all((ROOT / include).exists() for include in includes)


def test_final_state_has_no_temporary_exception_mechanism() -> None:
    policy = (ROOT / "tests/architecture/logical_ownership_policy.py").read_text(
        encoding="utf-8"
    )
    assert "ArchitectureException" not in policy
    assert "load_architecture_exceptions" not in policy
    assert not (ROOT / "architecture/temporary-exceptions.json").exists()
    permanent = json.loads(PERMANENT_CLASSIFICATIONS.read_text(encoding="utf-8"))
    assert {item["classification_id"] for item in permanent["classifications"]} == {
        "dataset-builder-path-reader",
        "dataset-input-source-validation",
        "dataset-request-source-validation",
    }


def test_forwarding_safety_and_lazy_composition_are_exactly_classified() -> None:
    payload = json.loads(CLASSIFICATIONS.read_text(encoding="utf-8"))
    assert set(payload["forwarding_and_safety_modules"]) == set(
        PERMANENT_SPECIAL_MODULES
    )
    assert set(payload["lazy_loader_modules"]) == set(PERMANENT_LAZY_LOADERS)

    lazy_paths: set[str] = set()
    special_paths: set[str] = set()
    for path in PACKAGE_ROOT.rglob("*.py"):
        relative = path.relative_to(ROOT).as_posix()
        tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
        if any(
            isinstance(node, ast.FunctionDef) and node.name == "__getattr__"
            for node in tree.body
        ):
            lazy_paths.add(relative)
        docstring = ast.get_docstring(tree) or ""
        if (
            "compatibility route" in docstring
            or "re-exported from their direct owners" in docstring
        ):
            special_paths.add(relative)

    assert lazy_paths == set(PERMANENT_LAZY_LOADERS)
    assert special_paths == set(PERMANENT_SPECIAL_MODULES)


def test_retired_route_literals_are_exact_versioned_identifiers() -> None:
    payload = json.loads(CLASSIFICATIONS.read_text(encoding="utf-8"))
    expected = {
        path: tuple(values)
        for path, values in payload["retired_route_literals"].items()
    }
    source_by_path = {
        path.relative_to(ROOT).as_posix(): path.read_text(encoding="utf-8-sig")
        for path in PACKAGE_ROOT.rglob("*.py")
    }
    assert old_path_literal_inventory(source_by_path) == expected


def test_every_legacy_marker_is_bound_to_a_permanent_behavior_classification() -> None:
    payload = json.loads(CLASSIFICATIONS.read_text(encoding="utf-8"))
    classified = {
        path for paths in payload["legacy_behavior_modules"].values() for path in paths
    }
    observed = {
        path.relative_to(ROOT).as_posix()
        for path in PACKAGE_ROOT.rglob("*.py")
        if re.search(r"legacy", path.read_text(encoding="utf-8-sig"), re.IGNORECASE)
    }
    assert observed == classified


def test_owner_map_has_permanent_schema_and_rationales() -> None:
    with OWNER_MAP.open(encoding="utf-8", newline="") as handle:
        rows = tuple(csv.DictReader(handle))
    assert rows
    assert set(rows[0]) == {"source_path", "module", "owner", "rationale"}
    assert all(row["rationale"].strip() for row in rows)


def test_final_architecture_inventory_is_complete_and_current() -> None:
    text = FINAL_INVENTORY.read_text(encoding="utf-8")
    assert text == render_inventory(ROOT)
    assert "Tracked production files: **647**" in text
    assert "Tracked Python modules: **641**" in text
    assert "Owner-map entries: **641**" in text
