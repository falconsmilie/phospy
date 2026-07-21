from __future__ import annotations

import json
import os
import shlex
import subprocess
import sys
import tempfile
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

COLLECTION_OUTPUT_ENV = "PHOSPY_RELEASE_SELECTOR_COLLECTION_OUTPUT"

RELEASE_BLOCKING_STANDALONE_MARKERS = frozenset(
    {"release_gate", "golden", "reproducibility"}
)
PERFORMANCE_RELEASE_SELECTOR_MARKERS = frozenset({"performance", "release_gate"})


@dataclass(frozen=True, order=True)
class CollectedNode:
    nodeid: str
    markers: frozenset[str]

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> CollectedNode:
        nodeid = payload.get("nodeid")
        markers = payload.get("markers")
        if not isinstance(nodeid, str):
            raise ValueError(f"collection payload has invalid nodeid: {payload!r}")
        if not isinstance(markers, list) or not all(
            isinstance(marker, str) for marker in markers
        ):
            raise ValueError(f"collection payload has invalid markers: {payload!r}")
        return cls(nodeid=_normalize_nodeid(nodeid), markers=frozenset(markers))


@dataclass(frozen=True)
class CollectionTarget:
    name: str
    paths: tuple[str, ...]
    marker_expression: str | None


@dataclass(frozen=True)
class ReleaseSelectorCoverage:
    release_blocking_nodes: Mapping[str, CollectedNode]
    selected_nodes_by_target: Mapping[str, frozenset[str]]
    target_definitions: tuple[CollectionTarget, ...]

    @property
    def missing_nodes(self) -> dict[str, CollectedNode]:
        return compare_release_blocking_coverage(
            self.release_blocking_nodes.values(),
            self.selected_nodes_by_target,
        )

    def format_failure(self) -> str:
        missing = self.missing_nodes
        if not missing:
            return "All release-blocking nodes are selected by release targets."

        lines = [
            "Release-blocking tests were not selected by any authoritative "
            "release test target.",
            "",
            "Missing release-blocking nodes:",
        ]
        for nodeid, node in sorted(missing.items()):
            marker_list = ", ".join(sorted(node.markers)) or "<none>"
            lines.append(f"- {nodeid} markers=[{marker_list}]")

        lines.extend(
            [
                "",
                "Authoritative target inventory used for comparison:",
            ]
        )
        for target in self.target_definitions:
            selected = self.selected_nodes_by_target.get(target.name, frozenset())
            paths = " ".join(target.paths) if target.paths else "<configured testpaths>"
            marker_expression = target.marker_expression or "<none>"
            lines.append(
                "- "
                f"{target.name}: paths={paths}; "
                f"selector={marker_expression!r}; selected_nodes={len(selected)}"
            )

        lines.extend(
            [
                "",
                "Release-blocking inventory rules:",
                "- any collected node marked release_gate, golden, or reproducibility",
                "- any tests/parity node marked parity and not parity_diagnostic",
                "- any tests/performance node selected by performance or release_gate",
                f"- total_release_blocking_nodes={len(self.release_blocking_nodes)}",
            ]
        )
        return "\n".join(lines)


class PytestCollectionError(RuntimeError):
    pass


INVENTORY_TARGET = CollectionTarget(
    name="complete-release-blocking-inventory",
    paths=("tests",),
    marker_expression=None,
)
AUTHORITATIVE_RELEASE_TARGETS = (
    CollectionTarget(
        name="default non-parity target",
        paths=(),
        marker_expression="not parity",
    ),
    CollectionTarget(
        name="test-parity",
        paths=("tests/parity",),
        marker_expression="parity and not parity_diagnostic",
    ),
    CollectionTarget(
        name="test-performance",
        paths=("tests/performance",),
        marker_expression="performance or release_gate",
    ),
    CollectionTarget(
        name="test-release-gates",
        paths=("tests/release", "tests/golden"),
        marker_expression="release_gate or golden or reproducibility",
    ),
)


def pytest_collection_finish(session: Any) -> None:
    output_path = os.environ.get(COLLECTION_OUTPUT_ENV)
    if output_path is None:
        return

    nodes = [
        {
            "nodeid": _normalize_nodeid(item.nodeid),
            "markers": sorted({marker.name for marker in item.iter_markers()}),
        }
        for item in session.items
    ]
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps({"nodes": nodes}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def collect_pytest_nodes(
    repo_root: Path,
    target: CollectionTarget,
) -> tuple[CollectedNode, ...]:
    with tempfile.TemporaryDirectory() as tmp_dir:
        output_path = Path(tmp_dir) / "pytest-collection.json"
        env = os.environ.copy()
        env[COLLECTION_OUTPUT_ENV] = str(output_path)
        command = [
            sys.executable,
            "-m",
            "pytest",
            "--collect-only",
            "-q",
            "-o",
            "addopts=",
            "-p",
            "tools.testing.release_selector_coverage",
            "--color=no",
        ]
        if target.marker_expression is not None:
            command.extend(["-m", target.marker_expression])
        command.extend(target.paths)

        result = subprocess.run(
            command,
            cwd=repo_root,
            env=env,
            text=True,
            capture_output=True,
            check=False,
        )
        if result.returncode != 0:
            raise PytestCollectionError(
                "pytest collection failed for "
                f"{target.name!r}: {_format_command(command)}\n\n"
                f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"
            )
        if not output_path.is_file():
            raise PytestCollectionError(
                "pytest collection did not write node inventory for "
                f"{target.name!r}: {_format_command(command)}\n\n"
                f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"
            )
        payload = json.loads(output_path.read_text(encoding="utf-8"))

    nodes = payload.get("nodes")
    if not isinstance(nodes, list):
        raise PytestCollectionError(
            f"pytest collection inventory is missing a nodes list: {payload!r}"
        )
    return tuple(sorted(CollectedNode.from_payload(node) for node in nodes))


def audit_authoritative_release_coverage(repo_root: Path) -> ReleaseSelectorCoverage:
    inventory_nodes = collect_pytest_nodes(repo_root, INVENTORY_TARGET)
    release_blocking_nodes = release_blocking_inventory(inventory_nodes)
    selected_nodes_by_target = {
        target.name: frozenset(
            node.nodeid for node in collect_pytest_nodes(repo_root, target)
        )
        for target in AUTHORITATIVE_RELEASE_TARGETS
    }
    return ReleaseSelectorCoverage(
        release_blocking_nodes=release_blocking_nodes,
        selected_nodes_by_target=selected_nodes_by_target,
        target_definitions=AUTHORITATIVE_RELEASE_TARGETS,
    )


def release_blocking_inventory(
    nodes: Iterable[CollectedNode],
) -> dict[str, CollectedNode]:
    return {node.nodeid: node for node in nodes if is_release_blocking_node(node)}


def is_release_blocking_node(node: CollectedNode) -> bool:
    if node.markers & RELEASE_BLOCKING_STANDALONE_MARKERS:
        return True

    node_path = _node_path(node.nodeid)
    if (
        node_path.startswith("tests/parity/")
        and "parity" in node.markers
        and "parity_diagnostic" not in node.markers
    ):
        return True

    return node_path.startswith("tests/performance/") and bool(
        node.markers & PERFORMANCE_RELEASE_SELECTOR_MARKERS
    )


def compare_release_blocking_coverage(
    release_blocking_nodes: Iterable[CollectedNode],
    selected_nodes_by_target: Mapping[str, Iterable[str]],
) -> dict[str, CollectedNode]:
    inventory = {node.nodeid: node for node in release_blocking_nodes}
    selected_nodes: set[str] = set()
    for nodeids in selected_nodes_by_target.values():
        selected_nodes.update(_normalize_nodeid(nodeid) for nodeid in nodeids)
    return {
        nodeid: node
        for nodeid, node in sorted(inventory.items())
        if nodeid not in selected_nodes
    }


def _node_path(nodeid: str) -> str:
    return _normalize_nodeid(nodeid).split("::", maxsplit=1)[0]


def _normalize_nodeid(nodeid: str) -> str:
    return nodeid.replace("\\", "/")


def _format_command(command: list[str]) -> str:
    if os.name == "nt":
        return subprocess.list2cmdline(command)
    return shlex.join(command)
