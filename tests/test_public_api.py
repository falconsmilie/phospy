from __future__ import annotations

import pandas as pd

from phospy import KinaseActivityAnalyzer, KinaseWorkflow, PhosphoDataset, PhosRPipeline

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


def test_public_root_exports() -> None:
    import phospy

    expected = {
        "CoreOutputs",
        "CoreProcessingResult",
        "KinaseActivityAnalyzer",
        "KinaseActivityResult",
        "KinasePredictionResult",
        "KinaseWorkflow",
        "KinaseWorkflowResult",
        "PhosphoDataset",
        "PhosRPipeline",
        "SiteMatrixResult",
    }
    assert set(phospy.__all__) == expected


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


def test_pipeline_runs_with_class_api(tmp_path) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    pred_path = tmp_path / "predMat.csv"
    outdir = tmp_path / "out"

    make_total_df().to_csv(total_path, sep="\t", index=False)
    make_phospho_df().to_csv(phospho_path, sep="\t", index=False)
    make_pred_mat().to_csv(pred_path)

    pipeline = PhosRPipeline.from_files(
        total_path=total_path,
        phospho_path=phospho_path,
        pred_mat_path=pred_path,
        kinase_activity_threshold=0.6,
        kinase_activity_min_substrates=2,
        kinase_activity_top_n_substrates=2,
    )
    outputs = pipeline.run(outdir=outdir)

    assert outputs.kinase_activity is not None
    assert pipeline.kinase_activity_threshold == 0.6
    assert pipeline.kinase_activity_min_substrates == 2
    assert pipeline.kinase_activity_top_n_substrates == 2


def test_kinase_workflow_runs_with_class_api() -> None:
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 2.0, 10.0, 11.0],
            "sample_2": [1.5, 2.5, 10.5, 11.5],
        },
        index=["SITE_1", "SITE_2", "SITE_3", "SITE_4"],
    )

    workflow = KinaseWorkflow()
    result = workflow.run(
        phospho_matrix=phospho_matrix,
        substrate_map={
            "KINASE_A": ["SITE_1", "SITE_2"],
            "KINASE_B": ["SITE_3", "SITE_4"],
        },
        site_sequences={
            "SITE_1": "QQAAAAAYY",
            "SITE_2": "QQAAAAAYY",
            "SITE_3": "QQTTTTTYY",
            "SITE_4": "QQTTTTTYY",
        },
        motif_sequences={
            "KINASE_A": ["QQAAAAAYY", "QQAAAAAYY"],
            "KINASE_B": ["QQTTTTTYY", "QQTTTTTYY"],
        },
        min_substrates=2,
        min_motif_size=1,
        top=2,
        score_threshold=0.5,
        inclusion=1,
        ensemble_size=2,
        n_iterations=1,
        random_state=7,
    )

    assert result.prediction_result.pred_matrix.shape[1] == 2


def test_pipeline_propagates_max_unmatched_fraction(tmp_path) -> None:
    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )
    phospho_df = pd.DataFrame(
        {
            "uid": ["u1", "u2"],
            "gene_names": ["PRKACA", "BTK"],
            "gene_p_site": ["PRKACA_S339", "BTK_Y551"],
            "localization_prob": [0.95, 0.95],
            "centralized_sequence": ["AAAAAA", "BBBBBB"],
            "p_group1": [8.0, 6.0],
            "p_group2": [7.0, 5.0],
            "p_group3": [6.0, 4.0],
            "p_group4": [5.0, 3.0],
            "p_group5": [4.0, 2.0],
            "p_group6": [3.0, 1.0],
        }
    )

    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    total_df.to_csv(total_path, sep="	", index=False)
    phospho_df.to_csv(phospho_path, sep="	", index=False)

    pipeline = PhosRPipeline.from_files(
        total_path=total_path,
        phospho_path=phospho_path,
        max_unmatched_fraction=0.5,
    )
    outputs = pipeline.run()

    assert outputs.core.phospho_corrected.shape[0] == 1
    assert pipeline.max_unmatched_fraction == 0.5
