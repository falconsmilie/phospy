from __future__ import annotations

import pandas as pd
import pytest

from phospy.activities.models import KinaseActivityInputs, PredMatOverlapSummary
from phospy.activities.scoring import compute_activity_from_inputs
from phospy.errors.workflows import WorkflowBoundaryError


def _inputs(
    *,
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    threshold: float,
    min_substrates: int,
    top_n_substrates: int,
) -> KinaseActivityInputs:
    overlap_count = int(pred_mat.index.intersection(phospho_matrix.index).size)
    return KinaseActivityInputs(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=threshold,
        min_substrates=min_substrates,
        top_n_substrates=top_n_substrates,
        overlap_summary=PredMatOverlapSummary(
            overlap_count=overlap_count,
            pred_mat_rows=int(pred_mat.index.size),
            phospho_rows=int(phospho_matrix.index.size),
        ),
    )


def test_weighted_activity_ignores_missing_values_per_sample() -> None:
    pred_mat = pd.DataFrame(
        {"PRKACA": [0.9, 0.8, 0.7]},
        index=["A;S1;", "B;S2;", "C;S3;"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "phospho_corrected_1": [10.0, float("nan"), 1.0],
            "phospho_corrected_2": [20.0, 6.0, float("nan")],
        },
        index=["A;S1;", "B;S2;", "C;S3;"],
    )

    result = compute_activity_from_inputs(
        _inputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=0.6,
            min_substrates=3,
            top_n_substrates=3,
        )
    )

    assert result.weighted_activity.at[
        "PRKACA", "phospho_corrected_1"
    ] == pytest.approx(6.0625)
    assert result.weighted_activity.at[
        "PRKACA", "phospho_corrected_2"
    ] == pytest.approx((20 * 0.9 + 6 * 0.8) / (0.9 + 0.8))


def test_ksea_scoring_respects_threshold_and_min_substrates() -> None:
    pred_mat = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.8, 0.2],
            "AKT1": [0.95, 0.7, 0.61],
        },
        index=["A;S1;", "B;S2;", "C;S3;"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [2.0, 4.0, 6.0],
        },
        index=pred_mat.index.copy(),
    )

    result = compute_activity_from_inputs(
        _inputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=0.6,
            min_substrates=2,
            top_n_substrates=3,
        )
    )

    assert result.ksea_counts.to_dict() == {"AKT1": 3, "MAP2K6": 2}
    assert result.ksea_scores.at["MAP2K6", "sample_a"] == pytest.approx(1.5)
    assert result.ksea_scores.at["AKT1", "sample_b"] == pytest.approx(4.0)


def test_ksea_ignores_missing_values_per_sample() -> None:
    pred_mat = pd.DataFrame(
        {"PRKACA": [0.9, 0.8, 0.7]},
        index=["A;S1;", "B;S2;", "C;S3;"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "phospho_corrected_1": [10.0, float("nan"), 1.0],
            "phospho_corrected_2": [20.0, 6.0, float("nan")],
        },
        index=pred_mat.index.copy(),
    )

    result = compute_activity_from_inputs(
        _inputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=0.6,
            min_substrates=3,
            top_n_substrates=3,
        )
    )

    assert result.ksea_scores.at["PRKACA", "phospho_corrected_1"] == pytest.approx(5.5)
    assert result.ksea_scores.at["PRKACA", "phospho_corrected_2"] == pytest.approx(13.0)
    assert result.ksea_counts.to_dict() == {"PRKACA": 3}


def test_top_n_substrate_selection_is_deterministic_for_ties() -> None:
    pred_mat = pd.DataFrame(
        {"MAP2K6": [0.9, 0.9, 0.2]},
        index=["SITE_1", "SITE_2", "SITE_3"],
    )
    phospho_matrix = pd.DataFrame(
        {"sample_a": [10.0, 1.0, 100.0]},
        index=pred_mat.index.copy(),
    )

    result = compute_activity_from_inputs(
        _inputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=0.6,
            min_substrates=2,
            top_n_substrates=2,
        )
    )

    assert result.weighted_activity.at["MAP2K6", "sample_a"] == pytest.approx(5.5)


def test_target_count_and_target_table_outputs_are_consistent() -> None:
    pred_mat = pd.DataFrame(
        {
            "MAP2K6": [0.9, 0.5, 0.0],
            "AKT1": [0.4, 0.2, 0.1],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [2.0, 4.0, 6.0],
        },
        index=pred_mat.index.copy(),
    )

    result = compute_activity_from_inputs(
        _inputs(
            pred_mat=pred_mat,
            phospho_matrix=phospho_matrix,
            threshold=0.3,
            min_substrates=2,
            top_n_substrates=2,
        )
    )

    assert result.target_counts.to_dict() == {"MAP2K6": 2, "AKT1": 1}
    assert set(result.target_table.columns) == {"site_id", "kinase", "score"}
    assert int(result.target_table.shape[0]) == 3


def test_activity_stage_raises_when_all_candidates_are_filtered() -> None:
    pred_mat = pd.DataFrame(
        {"MAP2K6": [0.8, 0.7]},
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )
    phospho_matrix = pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [2.0, 4.0]},
        index=pred_mat.index.copy(),
    )

    with pytest.raises(
        WorkflowBoundaryError, match="seam=kinase.activity.valid_candidates"
    ):
        compute_activity_from_inputs(
            _inputs(
                pred_mat=pred_mat,
                phospho_matrix=phospho_matrix,
                threshold=0.95,
                min_substrates=3,
                top_n_substrates=2,
            )
        )
