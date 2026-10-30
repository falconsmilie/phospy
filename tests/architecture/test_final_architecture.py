from __future__ import annotations

import ast
import csv
import json
import os
import re
import subprocess
import sys
import tomllib
from pathlib import Path

from tools.architecture.final_policy import (
    EXACT_OWNER_EDGE_RULES,
    PERMANENT_LAZY_LOADERS,
    PERMANENT_SPECIAL_MODULES,
    OwnerEdgeRule,
    is_structural_forwarding_candidate,
    old_path_literal_inventory,
    owner_edge_policy_problems,
    tracked_source_inventory_problems,
)
from tools.architecture.import_graph import build_import_graph
from tools.architecture.inventory import render_inventory
from tools.architecture.logical_ownership import (
    OwnerGraphViews,
    load_module_owners,
    load_permanent_condensation_classifications,
    owner_graph_views,
)

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
    paths = tuple(item.decode("utf-8") for item in result.stdout.split(b"\0") if item)
    return tuple(path for path in paths if (ROOT / path).is_file())


def _live_owner_graph_views() -> OwnerGraphViews:
    entries = load_module_owners(OWNER_MAP)
    graph = build_import_graph(package_root=PACKAGE_ROOT, src_root=SRC_ROOT)
    return owner_graph_views(
        graph.module_edges,
        owners_by_module={entry.module: entry.owner for entry in entries},
        classifications=load_permanent_condensation_classifications(
            PERMANENT_CLASSIFICATIONS
        ),
    )


def test_complete_observed_owner_graph_matches_exact_adr_0053_rules() -> None:
    views = _live_owner_graph_views()

    assert (
        owner_edge_policy_problems(
            observed_edges=views.complete.edges,
            policy_filtered_edges=views.policy_filtered.edges,
        )
        == []
    )


def test_owner_edge_verifier_rejects_added_missing_and_broadened_edges() -> None:
    views = _live_owner_graph_views()
    observed = views.complete.edges
    policy_filtered = views.policy_filtered.edges
    added = set(observed) | {("activities", "kinase")}
    missing = set(observed) - {("activities", "frames")}
    broadened = dict(EXACT_OWNER_EDGE_RULES)
    broadened[("activities", "statistics")] = OwnerEdgeRule(
        "shared", "invented broad permission"
    )

    assert any(
        "has no exact rule: activities -> kinase" in problem
        for problem in owner_edge_policy_problems(
            observed_edges=added,
            policy_filtered_edges=set(policy_filtered) | {("activities", "kinase")},
        )
    )
    assert any(
        "unused or missing exact owner-edge permission: activities -> frames" in problem
        for problem in owner_edge_policy_problems(
            observed_edges=missing,
            policy_filtered_edges=set(policy_filtered) - {("activities", "frames")},
        )
    )
    assert any(
        "unused or missing exact owner-edge permission: activities -> statistics"
        in problem
        for problem in owner_edge_policy_problems(
            observed_edges=observed,
            policy_filtered_edges=policy_filtered,
            rules=broadened,
        )
    )


def test_owner_edge_verifier_rejects_reversed_cyclic_and_misclassified_edges() -> None:
    views = _live_owner_graph_views()
    observed = views.complete.edges
    policy_filtered = views.policy_filtered.edges
    reversed_edges = (set(observed) - {("frames", "errors")}) | {("errors", "frames")}
    cyclic_edges = set(observed) | {("errors", "frames")}
    misclassified = dict(EXACT_OWNER_EDGE_RULES)
    misclassified[("activities", "errors")] = OwnerEdgeRule(
        "capability", "incorrect fixture"
    )

    assert any(
        "has no exact rule: errors -> frames" in problem
        for problem in owner_edge_policy_problems(
            observed_edges=reversed_edges,
            policy_filtered_edges=(set(policy_filtered) - {("frames", "errors")})
            | {("errors", "frames")},
        )
    )
    assert any(
        "owner dependency cycle" in problem
        for problem in owner_edge_policy_problems(
            observed_edges=cyclic_edges,
            policy_filtered_edges=set(policy_filtered) | {("errors", "frames")},
        )
    )
    assert any(
        "incorrect owner-edge classification: activities -> errors" in problem
        for problem in owner_edge_policy_problems(
            observed_edges=observed,
            policy_filtered_edges=policy_filtered,
            rules=misclassified,
        )
    )


def test_final_tracked_source_inventory_is_complete_and_classified() -> None:
    tracked = _git_paths("src/phospy")
    owners = load_module_owners(OWNER_MAP)

    assert len(tracked) == 646
    assert sum(path.endswith(".py") for path in tracked) == 640
    assert len(owners) == 640
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
    production_roots = {
        path.split("/")[2] for path in tracked if len(path.split("/")) > 3
    }
    retired = {
        "contracts",
        "errors",
        "frames",
        "policies",
        "science",
        "tables",
        "validation",
        "workflows",
    }
    expected_intermediate_roots = {
        "activities",
        "api",
        "core",
        "data",
        "datasets",
        "differential",
        "enrichment",
        "evidence",
        "io",
        "kinase",
        "phosphosites",
        "preparation_records",
        "preprocessing",
        "provenance",
        "quantitative",
        "references",
        "run_lineage",
        "signalome",
        "statistics",
        "taxonomy",
    }
    assert not production_roots & retired
    assert production_roots == expected_intermediate_roots

    configuration = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    includes = configuration["tool"]["pyright"]["include"]
    assert {
        "src/phospy/errors",
        "src/phospy/frames",
        "src/phospy/policies",
        "src/phospy/tables",
        "src/phospy/validation",
    }.isdisjoint(includes)
    assert "src/phospy/core/__init__.py" in includes
    assert all((ROOT / include).exists() for include in includes)


def test_core_namespace_initializer_is_inert_and_narrowly_owned() -> None:
    initializer = PACKAGE_ROOT / "core/__init__.py"
    tree = ast.parse(initializer.read_text(encoding="utf-8"), filename=str(initializer))
    assert tree.body == []

    entries = load_module_owners(OWNER_MAP)
    core_entries = [entry for entry in entries if entry.module == "phospy.core"]
    assert len(core_entries) == 1
    assert core_entries[0].source_path == "src/phospy/core/__init__.py"
    assert core_entries[0].owner == "namespace_containers"


def test_final_state_has_no_temporary_exception_mechanism() -> None:
    policy = (ROOT / "tools/architecture/logical_ownership.py").read_text(
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


def test_reusable_architecture_support_has_one_maintainer_owner() -> None:
    support_paths = (
        ROOT / "tools/architecture/final_policy.py",
        ROOT / "tools/architecture/import_graph.py",
        ROOT / "tools/architecture/inventory.py",
        ROOT / "tools/architecture/logical_ownership.py",
    )
    retired_test_paths = (
        ROOT / "tests/architecture/final_architecture_policy.py",
        ROOT / "tests/architecture/import_graph_support.py",
        ROOT / "tests/architecture/logical_ownership_policy.py",
    )

    assert all(path.is_file() for path in support_paths)
    assert not any(path.exists() for path in retired_test_paths)
    for path in support_paths:
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source, filename=str(path))
        imported_modules = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        }
        imported_modules.update(
            node.module or ""
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
        )
        assert not any(
            module == "tests" or module.startswith("tests.")
            for module in imported_modules
        )
        assert "sys.path" not in source


def test_safety_and_lazy_composition_are_exactly_classified() -> None:
    payload = json.loads(CLASSIFICATIONS.read_text(encoding="utf-8"))
    assert set(payload["forwarding_and_safety_modules"]) == set(
        PERMANENT_SPECIAL_MODULES
    )
    assert set(payload["lazy_loader_modules"]) == set(PERMANENT_LAZY_LOADERS)

    lazy_paths: set[str] = set()
    for path in PACKAGE_ROOT.rglob("*.py"):
        relative = path.relative_to(ROOT).as_posix()
        tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
        if any(
            isinstance(node, ast.FunctionDef) and node.name == "__getattr__"
            for node in tree.body
        ):
            lazy_paths.add(relative)

    assert lazy_paths == set(PERMANENT_LAZY_LOADERS)


def test_structural_forwarding_candidates_have_exact_policy_dispositions() -> None:
    payload = json.loads(CLASSIFICATIONS.read_text(encoding="utf-8"))
    classified = payload["structural_forwarding_modules"]
    observed = {
        path.relative_to(ROOT).as_posix()
        for path in PACKAGE_ROOT.rglob("*.py")
        if is_structural_forwarding_candidate(
            path.read_text(encoding="utf-8-sig"),
            filename=str(path),
        )
    }

    assert observed == set(classified)
    assert {item["disposition"] for item in classified.values()} == {
        "supported_aggregate_facade",
        "deliberate_internal_facade",
    }
    assert {
        path
        for path, item in classified.items()
        if item["disposition"] == "supported_aggregate_facade"
    } == {
        "src/phospy/api/advanced/__init__.py",
        "src/phospy/api/__init__.py",
    }
    assert all(item["rationale"].strip() for item in classified.values())


def test_structural_forwarding_detection_ignores_docstring_wording_and_case() -> None:
    bodies = (
        '"""Compatibility route."""\n',
        '"""PERMANENT PUBLIC FACADE."""\n',
        "",
    )

    for prefix in bodies:
        assert is_structural_forwarding_candidate(
            prefix
            + "from phospy.core.frames.validation import require_dataframe\n"
            + "__all__ = ['require_dataframe']\n"
        )


def test_structural_forwarding_detection_covers_direct_imported_name_aliases() -> None:
    source = """
from phospy.core.frames.validation import require_dataframe as _require_dataframe

require_dataframe = _require_dataframe
__all__ = ["require_dataframe"]
"""

    assert is_structural_forwarding_candidate(source)


def test_structural_forwarding_detection_rejects_behavioral_modules() -> None:
    source = """
from phospy.core.frames.validation import require_dataframe

def validate(value):
    return require_dataframe(value)
"""

    assert not is_structural_forwarding_candidate(source)


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
    assert "Tracked production files: **646**" in text
    assert "Tracked Python modules: **640**" in text
    assert "Owner-map entries: **640**" in text
    assert "Complete statically observed owner graph (before classifications)" in text
    assert "`datasets, io, kinase, signalome`" in text
    assert "Policy-filtered owner graph used by the DAG gate" in text
    assert "No multi-owner strongly connected components; the view is acyclic." in text
    assert "`dataset-builder-path-reader`" in text
    assert "`dataset-input-source-validation`" in text
    assert "`dataset-request-source-validation`" in text


def test_architecture_inventory_cli_supports_module_and_direct_invocation(
    tmp_path: Path,
) -> None:
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    arguments = ("--repo-root", str(ROOT), "--check")

    module_result = subprocess.run(
        [
            sys.executable,
            "-m",
            "tools.architecture.generate_architecture_inventory",
            *arguments,
        ],
        cwd=ROOT,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )
    assert module_result.returncode == 0, module_result.stderr

    direct_result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "tools/architecture/generate_architecture_inventory.py"),
            *arguments,
        ],
        cwd=tmp_path,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )
    assert direct_result.returncode == 0, direct_result.stderr
