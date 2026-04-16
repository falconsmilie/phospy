from __future__ import annotations

import json
from collections import OrderedDict
from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api import (
    PredictionRunConfig,
    SignalomeRunConfig,
    SignalomeWorkflow,
)
from phospy.internal.kinase_workflows import KinaseWorkflow

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_DATA = ROOT / "examples" / "data"
WORKFLOW_FIXTURES = ROOT / "tests" / "fixtures" / "public_workflow_reference"

PREDMAT_BENCHMARKS = {
    "default": "predmat_default.csv",
    "r_parity": "predmat_r_parity.csv",
}

SIGNALOME_BENCHMARKS = {
    "modules": "signalome_modules.csv",
    "map_modules": "signalome_map_modules.csv",
    "network_nodes": "signalome_network_nodes.csv",
    "network_edges": "signalome_network_edges.csv",
}


def _load_demo_inputs() -> tuple[
    pd.DataFrame,
    dict[str, list[str]],
    dict[str, str],
    dict[str, list[str]],
]:
    phospho_matrix = pd.read_csv(
        EXAMPLE_DATA / "predmat_phospho_matrix.csv",
        index_col=0,
    )
    phospho_matrix.index = phospho_matrix.index.map(str)
    substrate_map = json.loads(
        (EXAMPLE_DATA / "predmat_substrate_map.json").read_text(encoding="utf-8")
    )
    site_sequences = json.loads(
        (EXAMPLE_DATA / "predmat_site_sequences.json").read_text(encoding="utf-8")
    )
    motif_sequences = json.loads(
        (EXAMPLE_DATA / "predmat_motif_sequences.json").read_text(encoding="utf-8")
    )
    return phospho_matrix, substrate_map, site_sequences, motif_sequences


def _run_public_predmat_workflow(*, svm_mode: str) -> pd.DataFrame:
    phospho_matrix, substrate_map, site_sequences, motif_sequences = _load_demo_inputs()
    prediction_config = PredictionRunConfig(
        min_substrates=2,
        min_motif_size=2,
        ensemble_size=3,
        top=4,
        score_threshold=0.75,
        inclusion=3,
        n_iterations=2,
        random_state=17,
    )
    result = KinaseWorkflow(flank_size=2, svm_mode=svm_mode).run(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        site_sequences=site_sequences,
        motif_sequences=motif_sequences,
        prediction_config=prediction_config,
    )
    return result.prediction_result.pred_mat_result.to_frame(copy=False)


def _run_public_signalome_workflow(
    *, svm_mode: str
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    phospho_matrix, substrate_map, site_sequences, motif_sequences = _load_demo_inputs()
    prediction_config = PredictionRunConfig(
        min_substrates=2,
        min_motif_size=2,
        ensemble_size=3,
        top=4,
        score_threshold=0.75,
        inclusion=3,
        n_iterations=2,
        random_state=17,
    )
    pred_mat_result = KinaseWorkflow(flank_size=2, svm_mode=svm_mode).run(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        site_sequences=site_sequences,
        motif_sequences=motif_sequences,
        prediction_config=prediction_config,
    )
    signalome_result = SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A", "KINASE_B"],
        site_to_protein={
            str(site_id): str(site_id) for site_id in phospho_matrix.index
        },
        config=SignalomeRunConfig(signalome_cutoff=0.5),
    )
    map_data = signalome_result.to_map_data()
    network_data = signalome_result.to_network_data()
    return (
        signalome_result.modules.to_frame(),
        map_data.modules(),
        network_data.nodes(),
        network_data.edges(),
    )


def _read_indexed_fixture(name: str) -> pd.DataFrame:
    return pd.read_csv(WORKFLOW_FIXTURES / name, index_col=0)


def _read_unindexed_fixture(name: str) -> pd.DataFrame:
    return pd.read_csv(WORKFLOW_FIXTURES / name)


@pytest.mark.parity
@pytest.mark.parametrize("svm_mode", ["default", "r_parity"])
def test_public_predmat_workflow_matches_committed_benchmark(
    svm_mode: str,
) -> None:
    actual = _run_public_predmat_workflow(svm_mode=svm_mode)
    expected = _read_indexed_fixture(PREDMAT_BENCHMARKS[svm_mode])

    pdt.assert_frame_equal(actual, expected)
    assert list(actual.columns) == ["KINASE_A", "KINASE_B"]
    assert list(actual.index) == [f"SITE_{index}" for index in range(1, 9)]

    dominant = actual.idxmax(axis=1).to_dict()
    assert dominant == {
        "SITE_1": "KINASE_A",
        "SITE_2": "KINASE_A",
        "SITE_3": "KINASE_A",
        "SITE_4": "KINASE_A",
        "SITE_5": "KINASE_B",
        "SITE_6": "KINASE_B",
        "SITE_7": "KINASE_B",
        "SITE_8": "KINASE_B",
    }


@pytest.mark.parity
def test_public_predmat_workflow_default_mode_is_order_invariant_end_to_end() -> None:
    phospho_matrix, substrate_map, site_sequences, motif_sequences = _load_demo_inputs()
    prediction_config = PredictionRunConfig(
        min_substrates=2,
        min_motif_size=2,
        ensemble_size=3,
        top=4,
        score_threshold=0.75,
        inclusion=3,
        n_iterations=2,
        random_state=17,
    )
    reference = (
        KinaseWorkflow(flank_size=2, svm_mode="default")
        .run(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            motif_sequences=motif_sequences,
            prediction_config=prediction_config,
        )
        .prediction_result.pred_mat_result.to_frame(copy=False)
    )

    reversed_items = list(substrate_map.items())[::-1]
    reordered_substrate_map = OrderedDict(reversed_items)
    reordered_motif_sequences = OrderedDict(
        (kinase, motif_sequences[kinase]) for kinase, _ in reversed_items
    )

    reordered = (
        KinaseWorkflow(flank_size=2, svm_mode="default")
        .run(
            phospho_matrix=phospho_matrix,
            substrate_map=reordered_substrate_map,
            site_sequences=site_sequences,
            motif_sequences=reordered_motif_sequences,
            prediction_config=prediction_config,
        )
        .prediction_result.pred_mat_result.to_frame(copy=False)
    )

    pdt.assert_frame_equal(
        reference.sort_index().sort_index(axis=1),
        reordered.sort_index().sort_index(axis=1),
    )


@pytest.mark.parity
@pytest.mark.parametrize("svm_mode", ["default", "r_parity"])
def test_public_signalome_workflow_matches_committed_benchmark(
    svm_mode: str,
) -> None:
    actual_modules, actual_map_modules, actual_nodes, actual_edges = (
        _run_public_signalome_workflow(svm_mode=svm_mode)
    )

    expected_modules = _read_indexed_fixture(SIGNALOME_BENCHMARKS["modules"])
    expected_map_modules = _read_indexed_fixture(SIGNALOME_BENCHMARKS["map_modules"])
    expected_nodes = _read_indexed_fixture(SIGNALOME_BENCHMARKS["network_nodes"])
    expected_edges = _read_unindexed_fixture(SIGNALOME_BENCHMARKS["network_edges"])

    expected_modules.index.name = actual_modules.index.name
    expected_modules.columns.name = actual_modules.columns.name
    expected_map_modules.index.name = actual_map_modules.index.name
    expected_nodes.index.name = actual_nodes.index.name

    pdt.assert_frame_equal(actual_modules, expected_modules)
    pdt.assert_frame_equal(actual_map_modules, expected_map_modules)
    pdt.assert_frame_equal(actual_nodes, expected_nodes)
    pdt.assert_frame_equal(actual_edges, expected_edges, check_dtype=False)

    assert actual_map_modules["dominant_kinase"].to_dict() == {
        1: "KINASE_A",
        2: "KINASE_B",
    }
    assert actual_nodes["n_substrates"].to_dict() == {
        "KINASE_A": 4,
        "KINASE_B": 4,
    }
    assert actual_edges.empty
