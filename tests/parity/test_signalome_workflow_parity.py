from __future__ import annotations

from functools import lru_cache

import pandas.testing as pdt
import pytest

from phospy import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    ReferencePreset,
    SignalomeConfig,
    SignalomeWorkflow,
    SignalomeWorkflowRequest,
)
from tests.support.rewrite_fixture_data import (
    build_rat_l6_dataset,
    load_signalome_rewrite_l6_contract,
    load_signalome_rewrite_l6_module_assignments_selected,
    load_signalome_rewrite_l6_modules,
    load_signalome_rewrite_l6_network_edges_selected,
    load_signalome_rewrite_l6_network_nodes,
)

# Fixture provenance:
# - Input dataset fixture:
#   tests/fixtures/rewrite_parity/r_reference_l6/l6_phospho_matrix.csv
# - Expected signalome regression tables:
#   tests/fixtures/public_workflow_reference/signalome_rewrite_l6_*.{csv,json}
# These expectations are scoped to the supported rewrite lane only:
# dataset -> KinaseWorkflow -> SignalomeWorkflow with the config in
# `_run_signalome_l6_supported_slice`.
pytestmark = pytest.mark.parity


@lru_cache(maxsize=1)
def _run_signalome_l6_supported_slice():
    dataset = build_rat_l6_dataset(n_sites=260)
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(top_k=6, ensemble_size=12),
            activity_config=None,
        )
    )
    return SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(substrate_support_cutoff=0.5),
        )
    )


def test_signalome_module_assignments_match_selected_l6_regression_points() -> None:
    # Contract: selected point checks + structural distribution checks.
    assignments = _run_signalome_l6_supported_slice().module_assignments.table
    contract = load_signalome_rewrite_l6_contract()

    assert assignments.index.name == "site_id"
    assert int(assignments.shape[0]) == int(contract["n_assignments"])
    assert {
        "protein_id",
        "module_id",
        "top_kinase",
        "top_score",
        "top_kinase_candidates",
        "top_kinase_weights",
        "top_kinase_tie_count",
        "top_kinase_is_ambiguous",
        "top_kinase_selection_policy",
        "module_top_kinase",
        "module_top_kinase_candidates",
        "module_top_kinase_tie_count",
        "module_top_kinase_is_ambiguous",
        "module_top_kinase_selection_policy",
    }.issubset(set(assignments.columns))

    expected_points = load_signalome_rewrite_l6_module_assignments_selected()
    observed_points = assignments.loc[expected_points.index, expected_points.columns]
    pdt.assert_frame_equal(
        observed_points,
        expected_points,
        check_dtype=False,
    )

    expected_module_id_counts = {
        int(module_id): int(count)
        for module_id, count in dict(contract["module_id_counts"]).items()
    }
    observed_module_id_counts = (
        assignments.loc[:, "module_id"]
        .value_counts()
        .sort_index()
        .astype("int64")
        .to_dict()
    )
    assert observed_module_id_counts == expected_module_id_counts

    expected_tie_distribution = {
        int(tie_count): int(count)
        for tie_count, count in dict(contract["tie_count_distribution"]).items()
    }
    observed_tie_distribution = (
        assignments.loc[:, "top_kinase_tie_count"]
        .value_counts()
        .sort_index()
        .astype("int64")
        .to_dict()
    )
    assert observed_tie_distribution == expected_tie_distribution
    assert int(assignments.loc[:, "top_kinase_is_ambiguous"].sum()) == int(
        contract["ambiguous_assignment_count"]
    )


def test_signalome_modules_match_l6_fixture_table_exactly() -> None:
    # Contract: exact equality for this small deterministic module table.
    modules = _run_signalome_l6_supported_slice().signalome_modules.table
    expected_modules = load_signalome_rewrite_l6_modules()
    contract = load_signalome_rewrite_l6_contract()

    assert int(modules.shape[0]) == int(contract["n_modules"])
    assert int(modules.shape[1]) == int(contract["n_module_kinases"])
    pdt.assert_frame_equal(modules, expected_modules, check_dtype=False)
    assert (modules.sum(axis=1) - 100.0).abs().le(0.01).all()


def test_signalome_network_nodes_match_l6_fixture_counts_and_selected_rows() -> None:
    # Contract: structural comparison + selected row equality (not full graph parity).
    nodes = _run_signalome_l6_supported_slice().kinase_network.nodes
    assert nodes is not None
    expected_nodes = load_signalome_rewrite_l6_network_nodes()
    contract = load_signalome_rewrite_l6_contract()

    assert int(nodes.shape[0]) == int(contract["n_nodes"])
    assert {"degree", "n_substrates"} == set(nodes.columns)

    selected_kinases = ["AKT1", "PRKAA1", "RPS6KB1", "Yang.Ampk", "Yang.Akt"]
    pdt.assert_frame_equal(
        nodes.loc[selected_kinases, :],
        expected_nodes.loc[selected_kinases, :],
        check_dtype=False,
    )


def test_signalome_network_edges_match_l6_fixture_pairs_and_sign_counts() -> None:
    # Contract: structural counts + selected pair correlations with tolerance.
    edges = _run_signalome_l6_supported_slice().kinase_network.edges
    expected_edges = load_signalome_rewrite_l6_network_edges_selected()
    contract = load_signalome_rewrite_l6_contract()

    assert int(edges.shape[0]) == int(contract["n_edges"])
    assert {"source_kinase", "target_kinase", "correlation"} == set(edges.columns)
    assert int((edges.loc[:, "correlation"] < 0.0).sum()) == int(
        contract["negative_edge_count"]
    )
    assert int((edges.loc[:, "correlation"] > 0.0).sum()) == int(
        contract["positive_edge_count"]
    )

    expected_indexed = expected_edges.set_index(["source_kinase", "target_kinase"])
    observed_indexed = edges.set_index(["source_kinase", "target_kinase"]).loc[
        expected_indexed.index
    ]
    for pair, expected_row in expected_indexed.iterrows():
        assert observed_indexed.at[pair, "correlation"] == pytest.approx(
            expected_row["correlation"],
            rel=1e-9,
            abs=1e-12,
        )
