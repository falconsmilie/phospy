from __future__ import annotations

import pandas as pd
import pytest

from phosrpy import KinaseWorkflow, run_kinase_workflow


def make_workflow_inputs() -> tuple[
    pd.DataFrame,
    dict[str, list[str]],
    dict[str, str],
    dict[str, list[str]],
]:
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 1.1, 0.9, 1.2, 3.0, 2.9, 3.1, 2.8],
            "sample_2": [2.0, 2.1, 1.9, 2.2, 2.0, 2.1, 1.9, 2.2],
            "sample_3": [3.0, 3.1, 2.9, 3.2, 1.0, 1.1, 0.9, 1.2],
        },
        index=[f"SITE_{i}" for i in range(1, 9)],
    )
    substrate_map = {
        "KINASE_A": ["SITE_1", "SITE_2", "SITE_3", "SITE_4"],
        "KINASE_B": ["SITE_5", "SITE_6", "SITE_7", "SITE_8"],
    }
    site_sequences = {
        "SITE_1": "QQAAAAAYY",
        "SITE_2": "QQAAAAAYY",
        "SITE_3": "QQAAAAAYY",
        "SITE_4": "QQAAAAAYY",
        "SITE_5": "QQTTTTTYY",
        "SITE_6": "QQTTTTTYY",
        "SITE_7": "QQTTTTTYY",
        "SITE_8": "QQTTTTTYY",
    }
    motif_sequences = {
        "KINASE_A": ["QQAAAAAYY", "QQAAAAAYY", "QQAAAAAYY"],
        "KINASE_B": ["QQTTTTTYY", "QQTTTTTYY", "QQTTTTTYY"],
    }
    return phospho_matrix, substrate_map, site_sequences, motif_sequences


def test_run_kinase_workflow_runs_native_end_to_end_path() -> None:
    phospho_matrix, substrate_map, site_sequences, motif_sequences = (
        make_workflow_inputs()
    )

    result = run_kinase_workflow(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        site_sequences=site_sequences,
        motif_sequences=motif_sequences,
        min_substrates=2,
        min_motif_size=2,
        ensemble_size=3,
        top=4,
        score_threshold=0.75,
        inclusion=3,
        n_iterations=2,
        random_state=17,
        flank_size=2,
    )

    assert list(result.profile_result.profile_matrix.index) == ["KINASE_A", "KINASE_B"]
    assert list(result.motif_result.motif_scores.columns) == ["KINASE_A", "KINASE_B"]
    assert list(result.scoring_result.combined_scores.columns) == [
        "KINASE_A",
        "KINASE_B",
    ]
    assert list(result.prediction_result.pred_matrix.columns) == [
        "KINASE_A",
        "KINASE_B",
    ]
    assert set(result.prediction_result.substrate_list) == {"KINASE_A", "KINASE_B"}
    assert (
        result.prediction_result.pred_matrix.loc[
            ["SITE_1", "SITE_2", "SITE_3", "SITE_4"], "KINASE_A"
        ].mean()
        > result.prediction_result.pred_matrix.loc[
            ["SITE_5", "SITE_6", "SITE_7", "SITE_8"], "KINASE_A"
        ].mean()
    )


def test_kinase_workflow_can_fall_back_to_profile_only_prediction() -> None:
    phospho_matrix, substrate_map, _, _ = make_workflow_inputs()
    workflow = KinaseWorkflow()

    result = workflow.run(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        motif_sequences=None,
        allow_profile_only_fallback=True,
        ensemble_size=3,
        top=4,
        score_threshold=0.75,
        inclusion=3,
        n_iterations=2,
        random_state=23,
    )

    assert result.motif_result is None
    assert result.scoring_result.combined_scores is None
    assert list(result.prediction_result.pred_matrix.columns) == [
        "KINASE_A",
        "KINASE_B",
    ]


def test_kinase_workflow_requires_motif_sequences_without_profile_fallback() -> None:
    phospho_matrix, substrate_map, _, _ = make_workflow_inputs()
    workflow = KinaseWorkflow()

    with pytest.raises(ValueError, match="motif_sequences are required"):
        workflow.run(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
        )


def test_kinase_workflow_rejects_empty_substrate_map() -> None:
    phospho_matrix, _, _, _ = make_workflow_inputs()
    workflow = KinaseWorkflow()

    with pytest.raises(ValueError, match="substrate_map must not be empty"):
        workflow.run(
            phospho_matrix=phospho_matrix,
            substrate_map={},
            motif_sequences=None,
            allow_profile_only_fallback=True,
        )


def test_kinase_workflow_rejects_empty_motif_sequences_mapping() -> None:
    phospho_matrix, substrate_map, site_sequences, _ = make_workflow_inputs()
    workflow = KinaseWorkflow()

    with pytest.raises(ValueError, match="motif_sequences must not be empty"):
        workflow.run(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            motif_sequences={},
        )


def test_kinase_workflow_requires_site_sequences_when_motifs_are_provided() -> None:
    phospho_matrix, substrate_map, _, motif_sequences = make_workflow_inputs()
    workflow = KinaseWorkflow()

    with pytest.raises(ValueError, match="site_sequences are required"):
        workflow.run(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=None,
            motif_sequences=motif_sequences,
        )
