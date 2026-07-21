from __future__ import annotations

from tools.testing.release_selector_coverage import (
    CollectedNode,
    compare_release_blocking_coverage,
    is_release_blocking_node,
    release_blocking_inventory,
)


def _node(nodeid: str, *markers: str) -> CollectedNode:
    return CollectedNode(nodeid=nodeid, markers=frozenset(markers))


def test_compare_release_blocking_coverage_reports_missing_node_with_markers() -> None:
    covered = _node("tests/release/test_policy.py::test_covered", "release_gate")
    missing = _node("tests/golden/test_fixture.py::test_missing", "golden")

    missing_nodes = compare_release_blocking_coverage(
        [covered, missing],
        {
            "test-release-gates": [covered.nodeid],
            "test-performance": [],
        },
    )

    assert missing_nodes == {missing.nodeid: missing}
    assert missing_nodes[missing.nodeid].markers == frozenset({"golden"})


def test_parity_diagnostic_node_is_not_blocking_only_because_it_is_in_parity_suite() -> (
    None
):
    diagnostic = _node(
        "tests/parity/test_l6_prediction_parity.py::test_diagnostic",
        "parity",
        "parity_diagnostic",
    )

    assert is_release_blocking_node(diagnostic) is False


def test_non_diagnostic_parity_node_under_parity_suite_is_blocking() -> None:
    parity = _node(
        "tests/parity/test_prediction_science_parity.py::test_threshold",
        "parity",
    )

    assert is_release_blocking_node(parity) is True


def test_golden_only_node_under_golden_suite_is_release_blocking() -> None:
    golden = _node("tests/golden/test_fixture.py::test_schema", "golden")

    assert release_blocking_inventory([golden]) == {golden.nodeid: golden}


def test_reproducibility_only_node_is_release_blocking() -> None:
    reproducibility = _node(
        "tests/golden/test_fixture.py::test_replay",
        "reproducibility",
    )

    assert release_blocking_inventory([reproducibility]) == {
        reproducibility.nodeid: reproducibility
    }
