from __future__ import annotations

import numpy as np
import pandas as pd
import pandas.testing as pdt

from phospy.activities import (
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


def _reference_weighted_kinase_activity(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    top_n_substrates: int,
    min_substrates: int,
) -> pd.DataFrame:
    kinases = pred_mat.columns.tolist()
    samples = phospho_matrix.columns.tolist()
    kinase_mat = pd.DataFrame(index=kinases, columns=samples, dtype=float)

    available_sites = set(phospho_matrix.index)

    for kinase in kinases:
        top_substrates = pred_mat[kinase].nlargest(top_n_substrates)
        substrates = [site for site in top_substrates.index if site in available_sites]
        if len(substrates) < min_substrates:
            continue

        weights = top_substrates.loc[substrates].to_numpy(dtype=float)
        weight_sum = float(weights.sum())
        if weight_sum <= 0.0:
            continue

        values = phospho_matrix.loc[substrates, samples].to_numpy(dtype=float)
        valid_mask = ~np.isnan(values)
        broadcast_weights = weights[:, np.newaxis]
        weighted_values = np.where(valid_mask, values * broadcast_weights, 0.0)
        weight_totals = np.where(valid_mask, broadcast_weights, 0.0).sum(axis=0)
        result = np.full(values.shape[1], np.nan, dtype=float)
        np.divide(
            weighted_values.sum(axis=0),
            weight_totals,
            out=result,
            where=weight_totals > 0.0,
        )
        if np.isnan(result).all():
            continue
        kinase_mat.loc[kinase, :] = result

    return kinase_mat.dropna(how="all")


def _reference_ksea_scores(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    threshold: float,
    min_substrates: int,
) -> tuple[pd.DataFrame, pd.Series]:
    common_sites = pred_mat.index.intersection(phospho_matrix.index)
    aligned_pred_mat = pred_mat.loc[common_sites]
    aligned_matrix = phospho_matrix.loc[common_sites]
    substrate_mask = aligned_pred_mat.gt(threshold)
    substrate_counts = substrate_mask.sum(axis=0)
    candidate_kinases = substrate_counts.index[substrate_counts >= min_substrates]

    if len(candidate_kinases) == 0:
        empty_scores = pd.DataFrame(columns=list(phospho_matrix.columns), dtype=float)
        empty_counts = pd.Series(dtype=int, name="n_substrates")
        empty_counts.index.name = "kinase"
        return empty_scores, empty_counts

    score_dict: dict[str, pd.Series] = {}
    counts: dict[str, int] = {}

    for kinase in candidate_kinases:
        selected_sites = aligned_pred_mat.index[substrate_mask[kinase].to_numpy()]
        kinase_scores = aligned_matrix.loc[selected_sites].mean(axis=0, skipna=True)
        if kinase_scores.isna().all():
            continue

        counts[kinase] = len(selected_sites)
        score_dict[kinase] = kinase_scores

    score_frame = pd.DataFrame.from_dict(score_dict, orient="index")
    if score_frame.empty:
        score_frame = pd.DataFrame(columns=list(phospho_matrix.columns), dtype=float)
    count_series = pd.Series(counts, name="n_substrates").sort_values(ascending=False)
    count_series.index.name = "kinase"
    return score_frame, count_series


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


def test_compute_weighted_kinase_activity_skips_zero_weight_kinases() -> None:
    pred_mat = pd.DataFrame(
        {
            "PRKACA": [0.0, 0.0, 0.0],
            "BTK": [0.7, 0.6, 0.5],
        },
        index=["A;S1;", "B;S2;", "C;S3;"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "phospho_corrected_1": [10.0, 4.0, 1.0],
            "phospho_corrected_2": [20.0, 6.0, 2.0],
        },
        index=["A;S1;", "B;S2;", "C;S3;"],
    )

    result = compute_weighted_kinase_activity(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        top_n_substrates=3,
        min_substrates=3,
    )

    assert list(result.index) == ["BTK"]


def test_weighted_activity_ignores_missing_values_per_sample() -> None:
    pred_mat = pd.DataFrame(
        {
            "PRKACA": [0.9, 0.8, 0.7],
        },
        index=["A;S1;", "B;S2;", "C;S3;"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "phospho_corrected_1": [10.0, float("nan"), 1.0],
            "phospho_corrected_2": [20.0, 6.0, float("nan")],
        },
        index=["A;S1;", "B;S2;", "C;S3;"],
    )

    result = compute_weighted_kinase_activity(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        top_n_substrates=3,
        min_substrates=3,
    )

    assert round(float(result.loc["PRKACA", "phospho_corrected_1"]), 6) == 6.0625
    assert round(float(result.loc["PRKACA", "phospho_corrected_2"]), 6) == round(
        (20 * 0.9 + 6 * 0.8) / (0.9 + 0.8),
        6,
    )


def test_ksea_scores_ignore_missing_values_per_sample() -> None:
    pred_mat = pd.DataFrame(
        {
            "PRKACA": [0.9, 0.8, 0.7],
        },
        index=["A;S1;", "B;S2;", "C;S3;"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "phospho_corrected_1": [10.0, float("nan"), 1.0],
            "phospho_corrected_2": [20.0, 6.0, float("nan")],
        },
        index=["A;S1;", "B;S2;", "C;S3;"],
    )

    scores, counts = compute_ksea_scores(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=0.6,
        min_substrates=3,
    )

    assert float(scores.loc["PRKACA", "phospho_corrected_1"]) == 5.5
    assert float(scores.loc["PRKACA", "phospho_corrected_2"]) == 13.0
    assert int(counts.loc["PRKACA"]) == 3


def test_activity_methods_drop_kinases_with_all_missing_substrate_values() -> None:
    pred_mat = pd.DataFrame(
        {
            "PRKACA": [0.9, 0.8, 0.7],
            "BTK": [0.2, 0.1, 0.05],
        },
        index=["A;S1;", "B;S2;", "C;S3;"],
    )
    phospho_matrix = pd.DataFrame(
        {
            "phospho_corrected_1": [float("nan"), float("nan"), float("nan")],
            "phospho_corrected_2": [float("nan"), float("nan"), float("nan")],
        },
        index=["A;S1;", "B;S2;", "C;S3;"],
    )

    weighted = compute_weighted_kinase_activity(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        top_n_substrates=3,
        min_substrates=3,
    )
    scores, counts = compute_ksea_scores(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=0.6,
        min_substrates=3,
    )

    assert weighted.empty
    assert scores.empty
    assert counts.empty


def test_weighted_activity_matches_reference_on_large_sparse_input() -> None:
    rng = np.random.default_rng(42)
    site_ids = [f"SITE_{i}" for i in range(300)]
    pred_mat = pd.DataFrame(
        rng.random((300, 12)),
        index=site_ids,
        columns=[f"K{i}" for i in range(12)],
    )
    phospho_index = site_ids[20:] + [f"EXTRA_{i}" for i in range(10)]
    phospho_matrix = pd.DataFrame(
        rng.normal(size=(290, 8)),
        index=phospho_index,
        columns=[f"sample_{i}" for i in range(8)],
    )
    phospho_matrix[rng.random(phospho_matrix.shape) < 0.15] = np.nan

    result = compute_weighted_kinase_activity(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        top_n_substrates=25,
        min_substrates=5,
    )
    expected = _reference_weighted_kinase_activity(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        top_n_substrates=25,
        min_substrates=5,
    )

    pdt.assert_frame_equal(result.sort_index(), expected.sort_index())


def test_ksea_scores_match_reference_on_large_sparse_input() -> None:
    rng = np.random.default_rng(7)
    common_sites = [f"SITE_{i}" for i in range(240)]
    pred_mat = pd.DataFrame(
        rng.random((260, 10)),
        index=common_sites + [f"ONLY_PRED_{i}" for i in range(20)],
        columns=[f"K{i}" for i in range(10)],
    )
    phospho_matrix = pd.DataFrame(
        rng.normal(size=(255, 6)),
        index=common_sites[10:] + [f"ONLY_PHOS_{i}" for i in range(25)],
        columns=[f"sample_{i}" for i in range(6)],
    )
    phospho_matrix[rng.random(phospho_matrix.shape) < 0.2] = np.nan

    scores, counts = compute_ksea_scores(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=0.65,
        min_substrates=4,
    )
    expected_scores, expected_counts = _reference_ksea_scores(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=0.65,
        min_substrates=4,
    )

    pdt.assert_frame_equal(scores.sort_index(), expected_scores.sort_index())
    pdt.assert_series_equal(counts.sort_index(), expected_counts.sort_index())
