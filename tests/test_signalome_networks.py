from __future__ import annotations

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas as pd

from phospy.signalomes import SignalomeNetworkData, build_signalome_result
from phospy.signalomes.networks import (
    SignalomeNetworkEdge,
    SignalomeNetworkNode,
    build_signalome_network_data,
)


def _build_signalome_result_with_network_edge():
    site_ids = [
        "PROTEIN_1;S1;",
        "PROTEIN_1;S2;",
        "PROTEIN_2;S3;",
        "PROTEIN_2;S4;",
    ]
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 1.0, 4.0, 4.0],
            "KINASE_B": [1.1, 1.1, 4.1, 4.1],
            "KINASE_C": [4.0, 4.0, 1.0, 1.0],
        },
        index=site_ids,
    )
    pred_mat = pd.DataFrame(
        {
            "KINASE_A": [0.95, 0.93, 0.20, 0.25],
            "KINASE_B": [0.90, 0.89, 0.94, 0.92],
            "KINASE_C": [0.10, 0.15, 0.91, 0.90],
        },
        index=site_ids,
    )
    expression_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 1.1, 3.0, 3.1],
            "sample_2": [1.2, 1.0, 2.9, 3.0],
            "sample_3": [0.9, 1.2, 3.2, 2.8],
        },
        index=site_ids,
    )
    return build_signalome_result(
        scoring_matrix=scoring_matrix,
        pred_mat=pred_mat,
        expression_matrix=expression_matrix,
        kinases_of_interest=["KINASE_A", "KINASE_B"],
        kinase_network_threshold=0.9,
        signalome_cutoff=0.75,
        module_count=2,
    )


def _load_example_module(path: Path):
    spec = spec_from_file_location(path.stem, path)
    assert spec is not None
    assert spec.loader is not None
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_signalome_result_to_network_data_returns_canonical_graph_tables() -> None:
    result = _build_signalome_result_with_network_edge()

    network_data = result.to_network_data()

    assert isinstance(network_data, SignalomeNetworkData)
    assert list(network_data.nodes().columns) == [
        "degree",
        "n_substrates",
        "module_count",
        "total_share_percent",
        "is_kinase_of_interest",
    ]
    assert list(network_data.edges().columns) == [
        "source_kinase",
        "target_kinase",
        "correlation",
        "shared_module_count",
        "shared_modules",
        "source_is_kinase_of_interest",
        "target_is_kinase_of_interest",
    ]
    assert network_data.adjacency().equals(result.kinase_adjacency_matrix)
    assert not network_data.adjacency().equals(result.kinase_correlation_matrix)

    assert network_data.nodes().loc["KINASE_A", "degree"] == 1
    assert network_data.nodes().loc["KINASE_A", "module_count"] == 1
    assert network_data.nodes().loc["KINASE_B", "module_count"] == 2
    assert bool(network_data.nodes().loc["KINASE_A", "is_kinase_of_interest"])
    assert not bool(network_data.nodes().loc["KINASE_C", "is_kinase_of_interest"])

    assert network_data.edges().to_dict("records") == [
        {
            "source_kinase": "KINASE_A",
            "target_kinase": "KINASE_B",
            "correlation": 1.0,
            "shared_module_count": 1,
            "shared_modules": "[1]",
            "source_is_kinase_of_interest": True,
            "target_is_kinase_of_interest": True,
        }
    ]

    assert isinstance(network_data.node_models[0], SignalomeNetworkNode)
    assert isinstance(network_data.edge_models[0], SignalomeNetworkEdge)
    assert network_data.edge_models[0].shared_modules == (1,)


def test_build_signalome_network_data_matches_result_method() -> None:
    result = _build_signalome_result_with_network_edge()

    via_method = result.to_network_data()
    via_function = build_signalome_network_data(result)

    pd.testing.assert_frame_equal(via_method.adjacency(), via_function.adjacency())
    pd.testing.assert_frame_equal(via_method.nodes(), via_function.nodes())
    pd.testing.assert_frame_equal(via_method.edges(), via_function.edges())
    assert via_method.node_models == via_function.node_models
    assert via_method.edge_models == via_function.edge_models


def test_signalome_network_data_to_frames_and_csv_exports_stable_tables(
    tmp_path: Path,
) -> None:
    result = _build_signalome_result_with_network_edge()

    network_data = result.to_network_data()
    frames = network_data.to_frames()

    assert list(frames) == [
        "signalome_network_nodes",
        "signalome_network_edges",
        "signalome_network_adjacency",
    ]

    written = network_data.to_csv(tmp_path)
    assert sorted(written) == [
        "signalome_network_adjacency",
        "signalome_network_edges",
        "signalome_network_nodes",
    ]

    reloaded_nodes = pd.read_csv(written["signalome_network_nodes"], index_col=0)
    reloaded_edges = pd.read_csv(written["signalome_network_edges"])
    reloaded_adjacency = pd.read_csv(
        written["signalome_network_adjacency"],
        index_col=0,
    )

    reloaded_nodes.index.name = network_data.nodes().index.name
    reloaded_adjacency.index.name = network_data.adjacency().index.name
    reloaded_adjacency.columns.name = network_data.adjacency().columns.name

    reloaded_nodes = reloaded_nodes.astype(network_data.nodes().dtypes.to_dict())
    reloaded_edges = reloaded_edges.astype(network_data.edges().dtypes.to_dict())
    reloaded_adjacency = reloaded_adjacency.astype(
        network_data.adjacency().dtypes.to_dict()
    )

    pd.testing.assert_frame_equal(reloaded_nodes, network_data.nodes())
    pd.testing.assert_frame_equal(reloaded_edges, network_data.edges())
    pd.testing.assert_frame_equal(reloaded_adjacency, network_data.adjacency())


def test_signalome_network_data_defaults_to_zero_copy_with_explicit_safe_copy() -> None:
    network_data = _build_signalome_result_with_network_edge().to_network_data()

    shared_frames = network_data.to_frames(copy=False)
    detached_frames = network_data.to_frames(copy=True)

    assert SignalomeNetworkData.__dataclass_params__.frozen is False
    assert shared_frames["signalome_network_nodes"] is network_data.node_table
    assert shared_frames["signalome_network_edges"] is network_data.edge_table
    assert shared_frames["signalome_network_adjacency"] is network_data.adjacency_matrix
    assert detached_frames["signalome_network_nodes"] is not network_data.node_table
    assert detached_frames["signalome_network_edges"] is not network_data.edge_table
    assert (
        detached_frames["signalome_network_adjacency"]
        is not network_data.adjacency_matrix
    )


def test_kinase_network_demo_runs_end_to_end(tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parents[1]
    module = _load_example_module(repo_root / "examples" / "kinase_network_demo.py")

    signalome_result, network_data, written = module.run_demo(tmp_path)

    assert set(written) == {
        "signalome_network_nodes",
        "signalome_network_edges",
        "signalome_network_adjacency",
    }
    assert signalome_result.signalome_modules.shape == (2, 3)
    assert network_data.nodes().shape[0] == 3
    assert written["signalome_network_edges"].exists()

    reloaded_edges = pd.read_csv(written["signalome_network_edges"])
    reloaded_edges = reloaded_edges.astype(network_data.edges().dtypes.to_dict())
    pd.testing.assert_frame_equal(reloaded_edges, network_data.edges())
