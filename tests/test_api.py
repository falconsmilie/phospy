from __future__ import annotations

import pandas as pd
import pytest

from phosrpy import (
    KinaseActivityAnalyzer,
    KinasePredictor,
    KinaseProfileBuilder,
    KinaseScorer,
    PhosphoDataset,
    PhosRPipeline,
    kinase_substrate_score,
)

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


def make_pred_mat() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "PRKACA": [0.9, 0.8, 0.7],
            "BTK": [0.2, 0.85, 0.75],
        },
        index=["PRKACA;S339;", "BTK;Y551;", "LYN;Y397;"],
    )


def test_phospho_dataset_process_core() -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
        comparisons=EXAMPLE_COMPARISONS,
    )
    result = dataset.process_core()

    assert sorted(result.total_unique["genes"].tolist()) == ["BTK", "LYN", "PRKACA"]
    assert "p_group1_group4" in result.phospho_corrected.columns
    assert "PRKACA;S339;" in result.site_matrix.matrix.index
    assert result.site_matrix.sequences.loc["PRKACA;S339;"] == "DDDDDD"


def test_kinase_activity_analyzer() -> None:
    analyzer = KinaseActivityAnalyzer(pred_mat=make_pred_mat())
    phospho_matrix = pd.DataFrame(
        {
            "phospho_corrected_1": [4.0, 4.0, 4.0],
            "phospho_corrected_2": [5.0, 5.0, 5.0],
        },
        index=["PRKACA;S339;", "BTK;Y551;", "LYN;Y397;"],
    )
    result = analyzer.analyze(
        phospho_matrix=phospho_matrix, threshold=0.6, min_substrates=2
    )

    assert set(result.weighted_activity.index) == {"PRKACA", "BTK"}
    assert int(result.target_counts.loc["BTK"]) == 2
    assert set(result.ksea_scores.columns) == {
        "phospho_corrected_1",
        "phospho_corrected_2",
    }


def test_pipeline_runs_with_class_api(tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    pred_path = tmp_path / "predMat.csv"
    outdir = tmp_path / "out"

    make_total_df().to_csv(total_path, sep="\t", index=False)
    make_phospho_df().to_csv(phospho_path, sep="\t", index=False)
    make_pred_mat().to_csv(pred_path)

    pipeline = PhosRPipeline.from_files(
        total_path=total_path, phospho_path=phospho_path, pred_mat_path=pred_path
    )
    outputs = pipeline.run(outdir=outdir)

    assert outputs.kinase_activity is not None
    assert (outdir / "df_phospho_corrected.csv").exists()
    assert (outdir / "kinase_activity_matrix.csv").exists()


def test_kinase_scorer_profile_api() -> None:
    scorer = KinaseScorer.from_profile_dict(
        {
            "PRKACA": pd.Series([1.0, 2.0, 3.0], index=["s1", "s2", "s3"]),
            "BTK": pd.Series([3.0, 2.0, 1.0], index=["s1", "s2", "s3"]),
        }
    )
    phospho_matrix = pd.DataFrame(
        {"s1": [1.0], "s2": [2.0], "s3": [3.0]},
        index=["PRKACA;S339;"],
    )

    result = scorer.score(phospho_matrix)

    assert float(result.profile_scores.loc["PRKACA;S339;", "PRKACA"]) == pytest.approx(
        1.0
    )
    assert result.combined_scores is None


def test_kinase_predictor_api() -> None:
    predictor = KinasePredictor()
    combined_scores = pd.DataFrame(
        {
            "PRKACA": [0.95, 0.91, 0.87, 0.82, 0.20, 0.10],
            "BTK": [0.10, 0.20, 0.25, 0.30, 0.90, 0.88],
        },
        index=[
            "SITE_1",
            "SITE_2",
            "SITE_3",
            "SITE_4",
            "SITE_5",
            "SITE_6",
        ],
    )

    result = predictor.predict(
        combined_scores=combined_scores,
        ensemble_size=2,
        top=4,
        score_threshold=0.8,
        inclusion=2,
        n_iterations=2,
        random_state=7,
    )

    assert set(result.pred_matrix.columns) == {"PRKACA", "BTK"}
    assert 0.0 <= float(result.pred_matrix.loc["SITE_1", "PRKACA"]) <= 1.0
    assert float(result.pred_matrix.loc["SITE_1", "PRKACA"]) > float(
        result.pred_matrix.loc["SITE_6", "PRKACA"]
    )


def test_kinase_profile_builder_and_scorer_api() -> None:
    phospho_matrix = pd.DataFrame(
        {"s1": [1.0, 3.0], "s2": [2.0, 4.0], "s3": [3.0, 5.0]},
        index=["PRKACA;S339;", "BTK;Y551;"],
    )
    builder = KinaseProfileBuilder()
    profile_result = builder.build(
        substrate_map={"PRKACA": ["PRKACA;S339;", "BTK;Y551;"]},
        phospho_matrix=phospho_matrix,
    )
    scorer = KinaseScorer.from_profile_result(profile_result)

    result = scorer.score(phospho_matrix)

    assert list(profile_result.profile_matrix.index) == ["PRKACA"]
    assert float(profile_result.profile_matrix.loc["PRKACA", "s1"]) == pytest.approx(
        2.0
    )
    assert list(result.profile_scores.columns) == ["PRKACA"]


def test_kinase_scorer_from_substrate_map_api() -> None:
    phospho_matrix = pd.DataFrame(
        {"s1": [1.0, 3.0], "s2": [2.0, 4.0], "s3": [3.0, 5.0]},
        index=["PRKACA;S339;", "BTK;Y551;"],
    )

    scorer = KinaseScorer.from_substrate_map(
        substrate_map={"PRKACA": ["PRKACA;S339;", "BTK;Y551;"]},
        phospho_matrix=phospho_matrix,
    )
    result = scorer.score(phospho_matrix)

    assert list(result.profile_scores.columns) == ["PRKACA"]


def test_kinase_substrate_score_api() -> None:
    phospho_matrix = pd.DataFrame(
        {"s1": [1.0, 3.0], "s2": [2.0, 2.0], "s3": [3.0, 1.0]},
        index=["SITE_A", "SITE_B"],
    )
    motif_scores = pd.DataFrame(
        {"KINASE_A": [0.8, 0.3], "KINASE_B": [0.4, 0.9]},
        index=phospho_matrix.index.copy(),
    )
    motif_sizes = pd.Series({"KINASE_A": 4, "KINASE_B": 2})

    result = kinase_substrate_score(
        substrate_map={"KINASE_A": ["SITE_A"], "KINASE_B": ["SITE_B"]},
        phospho_matrix=phospho_matrix,
        motif_scores=motif_scores,
        motif_sizes=motif_sizes,
    )

    assert list(result.combined_scores.columns) == ["KINASE_A", "KINASE_B"]
    assert list(result.ks_activity_matrix.index) == ["KINASE_A", "KINASE_B"]
