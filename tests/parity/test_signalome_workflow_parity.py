from __future__ import annotations

from functools import lru_cache

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import (
    KinaseWorkflow,
    SignalomeWorkflow,
)
from phospy.api import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    ReferencePreset,
    SignalomeConfig,
    SignalomeWorkflowRequest,
)
from phospy.signalomes.science import build_kinase_network
from tests.support.parity_reporting import format_shape, record_parity_metrics
from tests.support.rewrite_fixture_data import (
    build_rat_l6_dataset,
    load_signalome_rewrite_l6_contract,
    load_signalome_rewrite_l6_expanded_signalome,
    load_signalome_rewrite_l6_module_assignments,
    load_signalome_rewrite_l6_modules,
    load_signalome_rewrite_l6_network_edges,
    load_signalome_rewrite_l6_network_nodes,
    normalize_signalome_expanded_signalome_for_parity,
    normalize_signalome_module_assignments_for_parity,
    normalize_signalome_modules_for_parity,
    normalize_signalome_network_edges_for_parity,
    normalize_signalome_network_nodes_for_parity,
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

NUMERIC_RTOL = 1e-9
NUMERIC_ATOL = 1e-12


@lru_cache(maxsize=1)
def _run_signalome_l6_supported_slice():
    dataset = build_rat_l6_dataset(n_sites=260)
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                deterministic_max_selected_kinases=12,
                adaptive_ensemble_runs=12,
            ),
            activity_config=None,
        )
    )
    return SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(substrate_support_cutoff=0.5),
        )
    )


def _assert_supported_signalome_contract(contract: dict[str, object]) -> None:
    assert contract["supported_outputs"] == [
        "signalome_modules",
        "module_assignments.table",
        "kinase_network.nodes",
        "kinase_network.edges",
        "expanded_signalome",
    ]
    assert contract["comparison_rules"] == {
        "signalome_modules": {
            "row_order": "fixture_order_semantic",
            "column_order": "fixture_order_semantic",
            "comparison": "exact_equality",
        },
        "module_assignments.table": {
            "row_order": "site_id_sorted",
            "column_order": "fixture_order_semantic",
            "comparison": "numeric_tolerance",
            "numeric_rtol": NUMERIC_RTOL,
            "numeric_atol": NUMERIC_ATOL,
            "normalization": "canonicalize_collection_columns",
        },
        "kinase_network.nodes": {
            "row_order": "kinase_sorted",
            "column_order": "fixture_order_semantic",
            "comparison": "exact_equality",
        },
        "kinase_network.edges": {
            "row_order": "source_target_sorted",
            "column_order": "fixture_order_semantic",
            "comparison": "numeric_tolerance",
            "numeric_rtol": NUMERIC_RTOL,
            "numeric_atol": NUMERIC_ATOL,
        },
        "expanded_signalome": {
            "row_order": "kinase_row_kind_site_sorted",
            "column_order": "fixture_order_semantic",
            "comparison": "numeric_tolerance",
            "numeric_rtol": NUMERIC_RTOL,
            "numeric_atol": NUMERIC_ATOL,
        },
    }


def test_signalome_contract_declares_full_supported_output_parity_lane() -> None:
    contract = load_signalome_rewrite_l6_contract()
    _assert_supported_signalome_contract(contract)


def test_signalome_module_assignments_match_l6_full_fixture_table(
    request: pytest.FixtureRequest,
) -> None:
    contract = load_signalome_rewrite_l6_contract()
    _assert_supported_signalome_contract(contract)
    observed = normalize_signalome_module_assignments_for_parity(
        _run_signalome_l6_supported_slice().module_assignments.table
    )
    expected = normalize_signalome_module_assignments_for_parity(
        load_signalome_rewrite_l6_module_assignments()
    )

    assert observed.columns.tolist() == expected.columns.tolist()
    assert int(observed.shape[0]) == int(contract["n_assignments"])
    pdt.assert_frame_equal(
        observed,
        expected,
        check_dtype=False,
        check_exact=False,
        rtol=NUMERIC_RTOL,
        atol=NUMERIC_ATOL,
    )

    record_parity_metrics(
        request.config,
        family="signalome_workflow",
        metrics=[
            ("module assignments full-table parity", "pass"),
            ("module assignment table shape", format_shape(*observed.shape)),
            ("module assignment rows compared", int(observed.shape[0])),
        ],
    )


def test_signalome_modules_match_l6_full_fixture_table_exactly(
    request: pytest.FixtureRequest,
) -> None:
    contract = load_signalome_rewrite_l6_contract()
    _assert_supported_signalome_contract(contract)
    observed = normalize_signalome_modules_for_parity(
        _run_signalome_l6_supported_slice().signalome_modules.table
    )
    expected = normalize_signalome_modules_for_parity(
        load_signalome_rewrite_l6_modules()
    )

    assert observed.index.tolist() == expected.index.tolist()
    assert observed.columns.tolist() == expected.columns.tolist()
    assert int(observed.shape[0]) == int(contract["n_modules"])
    assert int(observed.shape[1]) == int(contract["n_module_kinases"])
    pdt.assert_frame_equal(
        observed,
        expected,
        check_dtype=False,
    )

    record_parity_metrics(
        request.config,
        family="signalome_workflow",
        metrics=[
            ("signalome modules full-table parity", "pass"),
            ("module table shape", format_shape(*observed.shape)),
            ("module rows compared", int(observed.shape[0])),
        ],
    )


def test_signalome_network_nodes_match_l6_full_fixture_table(
    request: pytest.FixtureRequest,
) -> None:
    contract = load_signalome_rewrite_l6_contract()
    _assert_supported_signalome_contract(contract)
    nodes = _run_signalome_l6_supported_slice().kinase_network.nodes
    assert nodes is not None
    observed = normalize_signalome_network_nodes_for_parity(nodes)
    expected = normalize_signalome_network_nodes_for_parity(
        load_signalome_rewrite_l6_network_nodes()
    )

    assert observed.columns.tolist() == expected.columns.tolist()
    assert int(observed.shape[0]) == int(contract["n_nodes"])
    pdt.assert_frame_equal(
        observed,
        expected,
        check_dtype=False,
    )

    record_parity_metrics(
        request.config,
        family="signalome_workflow",
        metrics=[
            ("network nodes full-table parity", "pass"),
            ("network node table shape", format_shape(*observed.shape)),
            ("network node rows compared", int(observed.shape[0])),
        ],
    )


def test_signalome_network_edges_match_l6_full_fixture_table_with_tolerance(
    request: pytest.FixtureRequest,
) -> None:
    contract = load_signalome_rewrite_l6_contract()
    _assert_supported_signalome_contract(contract)
    observed = normalize_signalome_network_edges_for_parity(
        _run_signalome_l6_supported_slice().kinase_network.edges
    )
    expected = normalize_signalome_network_edges_for_parity(
        load_signalome_rewrite_l6_network_edges()
    )

    assert observed.columns.tolist() == expected.columns.tolist()
    assert int(observed.shape[0]) == int(contract["n_edges"])
    pdt.assert_frame_equal(
        observed,
        expected,
        check_dtype=False,
        check_exact=False,
        rtol=NUMERIC_RTOL,
        atol=NUMERIC_ATOL,
    )
    assert int((observed.loc[:, "correlation"] > 0.0).sum()) == int(
        contract["positive_edge_count"]
    )
    assert int((observed.loc[:, "correlation"] < 0.0).sum()) == int(
        contract["negative_edge_count"]
    )

    record_parity_metrics(
        request.config,
        family="signalome_workflow",
        metrics=[
            ("network edges full-table parity", "pass"),
            ("network edge table shape", format_shape(*observed.shape)),
            ("network edge rows compared", int(observed.shape[0])),
        ],
    )


def test_signalome_expanded_signalome_matches_l6_full_fixture_table_with_tolerance(
    request: pytest.FixtureRequest,
) -> None:
    contract = load_signalome_rewrite_l6_contract()
    _assert_supported_signalome_contract(contract)
    expanded = _run_signalome_l6_supported_slice().expanded_signalome
    assert expanded is not None
    observed = normalize_signalome_expanded_signalome_for_parity(expanded)
    expected = normalize_signalome_expanded_signalome_for_parity(
        load_signalome_rewrite_l6_expanded_signalome()
    )

    assert observed.columns.tolist() == expected.columns.tolist()
    assert int(observed.shape[0]) == int(contract["n_expanded_rows"])
    pdt.assert_frame_equal(
        observed,
        expected,
        check_dtype=False,
        check_exact=False,
        rtol=NUMERIC_RTOL,
        atol=NUMERIC_ATOL,
    )
    assert observed.loc[:, "row_kind"].value_counts().sort_index().astype(
        "int64"
    ).to_dict() == dict(contract["expanded_row_kind_counts"])

    record_parity_metrics(
        request.config,
        family="signalome_workflow",
        metrics=[
            ("expanded signalome full-table parity", "pass"),
            ("expanded signalome table shape", format_shape(*observed.shape)),
            ("expanded signalome rows compared", int(observed.shape[0])),
        ],
    )


def test_signalome_network_policy_variants_match_fixed_matrix_expectations() -> None:
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0, 2.0, 3.0, 4.0],
            "K2": [4.0, 3.0, 2.0, 1.0],
            "K3": [1.0, 2.0, 2.0, 3.0],
        },
        index=pd.Index(["S1", "S2", "S3", "S4"], name="site_id"),
    )
    common_kwargs = {
        "downstream_score_matrix": downstream_scores,
        "kinase_order": ["K1", "K2", "K3"],
        "kinase_substrates": {"K1": (), "K2": (), "K3": ()},
        "threshold": 0.9,
    }

    positive_only_edges, _ = build_kinase_network(
        **common_kwargs,
        network_policy="positive_only",
    )
    absolute_threshold_edges, _ = build_kinase_network(
        **common_kwargs,
        network_policy="absolute_threshold",
    )
    signed_edges, _ = build_kinase_network(
        **common_kwargs,
        network_policy="signed",
    )

    assert positive_only_edges.to_dict("records") == [
        {
            "source_kinase": "K1",
            "target_kinase": "K3",
            "correlation": pytest.approx(0.9486832980505138),
        }
    ]
    assert absolute_threshold_edges.to_dict("records") == [
        {
            "source_kinase": "K1",
            "target_kinase": "K2",
            "correlation": pytest.approx(1.0),
        },
        {
            "source_kinase": "K1",
            "target_kinase": "K3",
            "correlation": pytest.approx(0.9486832980505138),
        },
        {
            "source_kinase": "K2",
            "target_kinase": "K3",
            "correlation": pytest.approx(0.9486832980505138),
        },
    ]
    assert signed_edges.to_dict("records") == [
        {
            "source_kinase": "K1",
            "target_kinase": "K2",
            "correlation": pytest.approx(-1.0),
        },
        {
            "source_kinase": "K1",
            "target_kinase": "K3",
            "correlation": pytest.approx(0.9486832980505138),
        },
        {
            "source_kinase": "K2",
            "target_kinase": "K3",
            "correlation": pytest.approx(-0.9486832980505138),
        },
    ]
