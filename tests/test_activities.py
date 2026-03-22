from __future__ import annotations

import pandas as pd

from phospy import (
    build_kinase_target_table,
    compute_ksea_scores,
    compute_weighted_kinase_activity,
    count_predicted_targets,
)


def make_pred_mat() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "PRKACA": [0.9, 0.8, 0.7, 0.1],
            "BTK": [0.2, 0.85, 0.75, 0.65],
        },
        index=["A;S1;", "B;S2;", "C;S3;", "D;S4;"],
    )


def make_phospho_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "phospho_corrected_1": [10.0, 4.0, 1.0, 8.0],
            "phospho_corrected_2": [20.0, 6.0, 2.0, 10.0],
        },
        index=["A;S1;", "B;S2;", "C;S3;", "D;S4;"],
    )


def test_compute_weighted_kinase_activity() -> None:
    pred_mat = make_pred_mat()
    phospho_matrix = make_phospho_matrix()
    result = compute_weighted_kinase_activity(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        top_n_substrates=3,
        min_substrates=3,
    )
    assert set(result.index) == {"PRKACA", "BTK"}
    assert round(float(result.loc["PRKACA", "phospho_corrected_1"]), 6) == round(
        (10 * 0.9 + 4 * 0.8 + 1 * 0.7) / (0.9 + 0.8 + 0.7), 6
    )


def test_compute_ksea_scores_and_target_counts() -> None:
    pred_mat = make_pred_mat()
    phospho_matrix = make_phospho_matrix()
    scores, counts = compute_ksea_scores(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=0.6,
        min_substrates=3,
    )
    assert set(scores.index) == {"PRKACA", "BTK"}
    assert float(scores.loc["PRKACA", "phospho_corrected_1"]) == 5.0
    target_counts = count_predicted_targets(pred_mat, threshold=0.6)
    assert int(target_counts.loc["BTK"]) == 3
    assert int(counts.loc["PRKACA"]) == 3


def test_build_kinase_target_table() -> None:
    table = build_kinase_target_table(make_pred_mat(), threshold=0.6)
    assert {"site_id", "kinase", "score"} <= set(table.columns)
    assert table.shape[0] == 6
