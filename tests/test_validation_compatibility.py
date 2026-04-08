from __future__ import annotations

import pandas as pd
import pytest

from phospy import KinaseActivityAnalyzer, PhosphoDataset
from phospy.validation import (
    validate_protein_correction_inputs,
    validate_workflow_request,
)
from phospy.validation.errors import InputCompatibilityError, TableSchemaError
from phospy.workflow import KinaseWorkflow


def test_kinase_activity_analyzer_rejects_zero_overlap() -> None:
    pred_mat = pd.DataFrame(
        {
            "PRKACA": [0.9],
        },
        index=["SITE_A"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0],
        },
        index=["SITE_B"],
    )

    with pytest.raises(InputCompatibilityError, match="no overlapping phosphosite IDs"):
        KinaseActivityAnalyzer().run(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
        )


def test_kinase_workflow_accepts_partial_site_sequence_coverage() -> None:
    workflow = KinaseWorkflow()
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 2.0],
            "sample_2": [1.1, 2.1],
        },
        index=["SITE_1", "SITE_2"],
    )

    result = workflow.run(
        phospho_matrix=phospho_matrix,
        substrate_map={"KINASE_A": ["SITE_1", "SITE_2"]},
        site_sequences={"SITE_1": "QQAAAAAYY"},
        motif_sequences={"KINASE_A": ["QQAAAAAYY", "QQAAAAAYY"]},
        min_substrates=1,
        min_motif_size=1,
        ensemble_size=2,
        top=1,
        score_threshold=0.0,
        inclusion=1,
        n_iterations=1,
        random_state=7,
    )

    assert list(result.motif_result.motif_scores.index) == ["SITE_1"]
    assert list(result.scoring_result.profile_scores.index) == ["SITE_1"]
    assert list(result.scoring_result.combined_scores.index) == ["SITE_1"]
    assert list(result.prediction_result.pred_matrix.index) == ["SITE_1"]


def test_validate_workflow_request_limits_sequence_validation_to_scoring_subset() -> (
    None
):
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 2.0],
            "sample_2": [1.1, 2.1],
        },
        index=["SITE_1", "SITE_2"],
    )

    validated = validate_workflow_request(
        phospho_matrix=phospho_matrix,
        substrate_map={"KINASE_A": ["SITE_1", "SITE_2"]},
        site_sequences={"SITE_1": "QQAAAAAYY"},
        motif_sequences={"KINASE_A": ["QQAAAAAYY", "QQAAAAAYY"]},
        min_substrates=1,
        min_motif_size=1,
        ensemble_size=2,
        top=1,
        score_threshold=0.0,
        inclusion=1,
        n_iterations=1,
        random_state=7,
    )

    assert validated.scoring_site_index == ("SITE_1",)


def test_kinase_workflow_rejects_zero_substrate_overlap() -> None:
    workflow = KinaseWorkflow()
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0],
        },
        index=["SITE_1"],
    )

    with pytest.raises(
        InputCompatibilityError,
        match="no overlap between substrate_map and phospho_matrix",
    ):
        workflow.run(
            phospho_matrix=phospho_matrix,
            substrate_map={"KINASE_A": ["SITE_9"]},
            motif_sequences=None,
            allow_profile_only_fallback=True,
            ensemble_size=2,
            top=2,
            score_threshold=0.5,
            inclusion=1,
            n_iterations=1,
        )


def _make_total_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "genes": ["PRKACA", "BTK"],
            "group1": [1.0, 2.0],
            "group2": [1.0, 2.0],
            "group3": [1.0, 2.0],
            "group4": [1.0, 2.0],
            "group5": [1.0, 2.0],
            "group6": [1.0, 2.0],
        }
    )


def _make_phospho_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "uid": ["u1", "u2"],
            "gene_names": ["PRKACA", "MISSING1"],
            "gene_p_site": ["PRKACA_S339", "MISSING1_S1"],
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


def test_core_processing_rejects_zero_gene_overlap_before_correction() -> None:
    dataset = PhosphoDataset(
        total_df=_make_total_df().loc[[1]].reset_index(drop=True),
        phospho_df=_make_phospho_df().loc[[0]].assign(gene_names="PRKACA_MISSING"),
    )

    with pytest.raises(
        InputCompatibilityError,
        match="no overlapping gene identifiers",
    ):
        dataset.preprocessing.run()


def test_core_processing_rejects_row_loss_before_correction() -> None:
    dataset = PhosphoDataset(
        total_df=_make_total_df(),
        phospho_df=_make_phospho_df(),
    )

    with pytest.raises(
        InputCompatibilityError,
        match=r"would drop 1 of 2 phosphosite rows \(50.0%\)",
    ):
        dataset.preprocessing.run(min_observed=1)


def test_core_processing_allows_row_loss_with_explicit_tolerance() -> None:
    dataset = PhosphoDataset(
        total_df=_make_total_df(),
        phospho_df=_make_phospho_df(),
    )

    result = dataset.preprocessing.run(min_observed=1, max_unmatched_fraction=0.5)

    assert result.phospho_corrected.shape[0] == 1
    assert result.phospho_corrected["gene_names"].tolist() == ["PRKACA"]


def test_kinase_activity_analyzer_rejects_insufficient_overlap_fraction() -> None:
    pred_mat = pd.DataFrame(
        {
            "PRKACA": [0.9],
        },
        index=["SITE_1"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0] * 20,
        },
        index=[f"SITE_{idx}" for idx in range(1, 21)],
    )

    with pytest.raises(
        InputCompatibilityError, match="insufficient overlapping phosphosite IDs"
    ):
        KinaseActivityAnalyzer().run(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
        )


def test_validate_protein_correction_inputs_reports_missing_value_columns() -> None:
    phospho_df = pd.DataFrame({"gene_names": ["PRKACA"]})
    total_df = pd.DataFrame({"genes": ["PRKACA"], "group1": [1.0]})

    with pytest.raises(
        TableSchemaError,
        match=(
            r"Protein correction inputs phospho input is missing required columns: p_group1"
        ),
    ):
        validate_protein_correction_inputs(
            phospho_df,
            total_df,
            phospho_gene_col="gene_names",
            total_gene_col="genes",
            phospho_cols=["p_group1"],
            protein_cols=["group1"],
        )


def test_validate_protein_correction_inputs_uses_same_normalization_as_merge() -> None:
    phospho_df = pd.DataFrame(
        {
            "gene_names": [" prkaca ", "BTK", "missing "],
            "p_group1": [8.0, 6.0, 2.0],
        }
    )
    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA", " btk "],
            "group1": [1.0, 2.0],
        }
    )

    with pytest.raises(
        InputCompatibilityError,
        match=r"would drop 1 of 3 phosphosite rows \(33.3%\).*MISSING",
    ):
        validate_protein_correction_inputs(
            phospho_df,
            total_df,
            phospho_gene_col="gene_names",
            total_gene_col="genes",
            phospho_cols=["p_group1"],
            protein_cols=["group1"],
            max_unmatched_fraction=0.0,
        )


def test_kinase_workflow_rejects_zero_site_sequence_overlap_when_motifs_are_enabled() -> (
    None
):
    workflow = KinaseWorkflow()
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0],
        },
        index=["SITE_1"],
    )

    with pytest.raises(
        InputCompatibilityError,
        match="sequence coverage required for scoring and prediction",
    ):
        workflow.run(
            phospho_matrix=phospho_matrix,
            substrate_map={"KINASE_A": ["SITE_1"]},
            site_sequences={"SITE_9": "QQAAAAAYY"},
            motif_sequences={"KINASE_A": ["QQAAAAAYY", "QQAAAAAYY"]},
            min_substrates=1,
            min_motif_size=1,
            ensemble_size=2,
            top=1,
            score_threshold=0.0,
            inclusion=1,
            n_iterations=1,
        )
