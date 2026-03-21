from __future__ import annotations

import pandas as pd
import pytest

from phosrpy import (
    KinaseActivityAnalyzer,
    KinaseMotifScorer,
    KinaseScorer,
    PhosphoDataset,
    PhosRPipeline,
    create_frequency_matrix,
    score_phosphosite_motifs,
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


def test_kinase_motif_scorer_api() -> None:
    scorer = KinaseMotifScorer.from_substrate_sequences(
        {
            "KINASE_A": ["AAAAA"],
            "KINASE_B": ["TTTTT"],
        },
        flank_size=2,
    )

    result = scorer.score_sequences(
        seqs={"SITE_A": "QQAAAAAYY", "SITE_B": "QQTTTTTYY"},
        min_motif_size=1,
    )

    assert float(result.motif_scores.loc["SITE_A", "KINASE_A"]) == pytest.approx(1.0)
    assert result.sequence_windows.loc["SITE_B"] == "TTTTT"


def test_score_phosphosite_motifs_api() -> None:
    result = score_phosphosite_motifs(
        seqs={"SITE_A": "QQAAAAAYY"},
        motif_frequency_matrices={
            "KINASE_A": create_frequency_matrix(["AAAAA"], flank_size=2),
        },
        motif_sizes=pd.Series({"KINASE_A": 4}, dtype=float),
        flank_size=2,
    )

    assert list(result.motif_scores.columns) == ["KINASE_A"]
    assert float(result.motif_scores.loc["SITE_A", "KINASE_A"]) != float("inf")


def test_phospho_dataset_process_core_does_not_multiply_rows_for_mixed_case_total_duplicates() -> (
    None
):
    total_df = pd.DataFrame(
        {
            "genes": ["Prkaca", "PRKACA", "Btk", "Lyn"],
            "group1": [1.0, 5.0, 2.0, 3.0],
            "group2": [1.0, 5.0, 2.0, 3.0],
            "group3": [1.0, 5.0, 2.0, 3.0],
            "group4": [1.0, 5.0, 2.0, 3.0],
            "group5": [1.0, 5.0, 2.0, 3.0],
            "group6": [1.0, 5.0, 2.0, 3.0],
        }
    )

    dataset = PhosphoDataset(
        total_df=total_df,
        phospho_df=make_phospho_df(),
        comparisons=EXAMPLE_COMPARISONS,
    )
    result = dataset.process_core()

    assert result.total_unique["genes"].is_unique
    assert sorted(result.total_unique["genes"].tolist()) == ["BTK", "LYN", "PRKACA"]
    assert len(result.phospho_corrected) == 4
    assert (result.phospho_corrected["gene_names"] == "PRKACA").sum() == 2
