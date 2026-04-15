from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd

from phospy.api import (
    PredictionRunConfig,
    SignalomeRunConfig,
    SignalomeWorkflow,
    SimpleKinaseWorkflow,
)
from phospy.datasets import AnalysisReadyPhosphoDataset, PhosphoDataset
from phospy.preprocessing import CorePreprocessingConfig

ROOT = Path(__file__).resolve().parents[1]
SIMPLE_WORKFLOW_FIXTURE_DIR = ROOT / "examples" / "data" / "simple_workflow"
EXAMPLE_COMPARISONS = [("group1", "group4"), ("group2", "group5"), ("group3", "group6")]


def make_total_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "genes": ["Prkaca", "Prkaca", "Btk", "Lyn"],
            "group1": [1.0, 5.0, 2.0, 3.0],
            "group2": [1.0, 5.0, 2.0, 3.0],
            "group3": [1.0, 5.0, 2.0, 3.0],
            "group4": [1.0, 5.0, 2.0, 3.0],
            "group5": [1.0, 5.0, 2.0, 3.0],
            "group6": [1.0, 5.0, 2.0, 3.0],
        }
    )


def make_phospho_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "uid": ["u1", "u2", "u3", "u4"],
            "gene_names": ["PRKACA", "BTK", "LYN", "PRKACA"],
            "gene_p_site": ["PRKACA_S339", "BTK_Y551", "LYN_Y397", "PRKACA_S339"],
            "localization_prob": [0.95, 0.95, 0.95, 0.95],
            "centralized_sequence": ["AAAAAA", "BBBBBB", "CCCCCC", "DDDDDD"],
            "p_group1": [8.0, 6.0, 7.0, 9.0],
            "p_group2": [7.0, 5.0, 6.0, 8.0],
            "p_group3": [6.0, 4.0, 5.0, 7.0],
            "p_group4": [5.0, 3.0, 4.0, 6.0],
            "p_group5": [4.0, 2.0, 3.0, 5.0],
            "p_group6": [3.0, 1.0, 2.0, 4.0],
        }
    )


def test_public_root_exports() -> None:
    import phospy

    assert phospy.__all__ == []
    assert not hasattr(phospy, "PhosphoDataset")
    assert not hasattr(phospy, "SimpleKinaseWorkflow")


def test_supported_public_surface_is_small_and_explicit() -> None:
    import phospy.api as api

    assert set(api.__all__) == {
        "DatasetLoadOptions",
        "KinaseActivityConfig",
        "PredictionRunConfig",
        "SignalomeRunConfig",
        "SignalomeWorkflow",
        "SimpleKinaseWorkflow",
    }
    assert not hasattr(api, "PredMatWorkflow")
    assert not hasattr(api, "KinaseWorkflow")


def test_removed_workflow_entrypoints_are_not_importable() -> None:
    assert importlib.util.find_spec("phospy.pipeline") is None
    assert importlib.util.find_spec("phospy.api.kinase_workflows") is None
    assert importlib.util.find_spec("phospy.matrices") is None
    assert importlib.util.find_spec("phospy.motifs") is None
    assert importlib.util.find_spec("phospy.profiles") is None
    assert importlib.util.find_spec("phospy.orchestration") is None


def test_phospho_dataset_preprocessing_run() -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
        comparisons=EXAMPLE_COMPARISONS,
    )
    result = dataset.preprocessing.run(config=CorePreprocessingConfig())

    assert sorted(result.total_unique["genes"].tolist()) == ["BTK", "LYN", "PRKACA"]
    assert "p_group1_group4" in result.phospho_corrected.columns
    assert "PRKACA;S339;" in result.site_matrix.matrix.index
    assert result.site_matrix.row_drop_stats["retained_rows"] == len(
        result.site_matrix.matrix
    )


def test_simple_kinase_workflow_runs_from_public_supported_path() -> None:
    total_df = pd.read_csv(SIMPLE_WORKFLOW_FIXTURE_DIR / "total.tsv", sep="\t")
    phospho_df = pd.read_csv(SIMPLE_WORKFLOW_FIXTURE_DIR / "phospho.tsv", sep="\t")

    with SimpleKinaseWorkflow(flank_size=7).run(
        total=total_df,
        phospho=phospho_df,
        species="rat",
        prediction_config=PredictionRunConfig(
            min_substrates=1,
            min_motif_size=1,
            ensemble_size=2,
            top=3,
            inclusion=2,
            n_iterations=2,
            random_state=7,
        ),
    ) as result:
        assert isinstance(result.analysis_ready_dataset, AnalysisReadyPhosphoDataset)
        assert result.reference_bundle.species == "rat"
        assert result.reference_bundle.source_metadata.reference == "l6_native"
        assert (
            result.analysis_ready_dataset.provenance.source == "simple kinase workflow"
        )
        assert result.pred_mat_result.to_frame(copy=False).shape == (5, 8)
        assert result.scoring_result is not None
        assert result.prediction_result is not None


def test_signalome_workflow_runs_from_simple_workflow_outputs() -> None:
    total_df = pd.read_csv(SIMPLE_WORKFLOW_FIXTURE_DIR / "total.tsv", sep="\t")
    phospho_df = pd.read_csv(SIMPLE_WORKFLOW_FIXTURE_DIR / "phospho.tsv", sep="\t")
    simple_result = SimpleKinaseWorkflow(flank_size=7).run(
        total=total_df,
        phospho=phospho_df,
        species="rat",
        prediction_config=PredictionRunConfig(
            min_substrates=1,
            min_motif_size=1,
            ensemble_size=2,
            top=3,
            inclusion=2,
            n_iterations=2,
            random_state=7,
        ),
    )
    site_to_protein = {
        str(site_id): str(site_id).split(";", 1)[0]
        for site_id in simple_result.analysis_ready_dataset.phospho_matrix.index
    }

    signalome_result = SignalomeWorkflow().run_from_analysis_ready(
        dataset=simple_result.analysis_ready_dataset,
        scoring_result=simple_result.scoring_result,
        prediction_result=simple_result.prediction_result,
        kinases_of_interest=list(
            simple_result.pred_mat_result.to_frame(copy=False).columns[:2]
        ),
        site_to_protein=site_to_protein,
        config=SignalomeRunConfig(signalome_cutoff=0.5),
    )

    assert signalome_result.modules.to_frame(copy=False).shape[0] >= 1
    assert signalome_result.network.nodes().shape[0] >= 1
    simple_result.close()
