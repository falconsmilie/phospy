from __future__ import annotations

import numpy as np
import pandas as pd


def compute_weighted_kinase_activity(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    top_n_substrates: int = 20,
    min_substrates: int = 3,
) -> pd.DataFrame:
    """Compute weighted downstream kinase activity scores.

    Missing phosphosite values are ignored on a per-sample basis. For each sample,
    weights are re-normalized across the observed substrates only. Kinases whose
    selected substrates contain no observed phosphosite values in any sample are
    omitted from the returned activity matrix.
    """

    kinases = pred_mat.columns.tolist()
    samples = phospho_matrix.columns.tolist()
    kinase_rows: list[np.ndarray] = []
    kinase_names: list[str] = []

    phospho_values = phospho_matrix.to_numpy(dtype=float, copy=False)
    phospho_site_positions = _build_site_position_lookup(phospho_matrix.index)

    for kinase in kinases:
        top_substrates = pred_mat[kinase].nlargest(top_n_substrates)
        substrate_positions: list[int] = []
        weights: list[float] = []

        for site_id, weight in top_substrates.items():
            site_position = phospho_site_positions.get(site_id)
            if site_position is None:
                continue
            substrate_positions.append(site_position)
            weights.append(float(weight))

        if len(substrate_positions) < min_substrates:
            continue

        weights_array = np.asarray(weights, dtype=float)
        weight_sum = float(weights_array.sum())
        if weight_sum <= 0.0:
            continue

        substrate_values = phospho_values[np.asarray(substrate_positions, dtype=int), :]
        weighted_values = _nan_aware_weighted_average(substrate_values, weights_array)
        if np.isnan(weighted_values).all():
            continue

        kinase_names.append(kinase)
        kinase_rows.append(weighted_values)

    if not kinase_rows:
        return pd.DataFrame(columns=samples, dtype=float)

    return pd.DataFrame(kinase_rows, index=kinase_names, columns=samples, dtype=float)


def build_kinase_target_table(
    pred_mat: pd.DataFrame,
    threshold: float = 0.6,
) -> pd.DataFrame:
    """Materialize a kinase-target edge table for reporting and export."""

    filtered = pred_mat.where(pred_mat > threshold)
    try:
        edges = filtered.stack(future_stack=True).rename("score").reset_index()
        edges = edges.loc[edges["score"].notna()]
    except TypeError:
        edges = filtered.stack(dropna=True).rename("score").reset_index()
    edges.columns = ["site_id", "kinase", "score"]
    return edges.sort_values(["kinase", "score"], ascending=[True, False])


def count_predicted_targets(
    pred_mat: pd.DataFrame,
    threshold: float = 0.6,
) -> pd.Series:
    """Count predicted kinase targets using matrix-native thresholding."""

    counts = _prediction_mask(pred_mat, threshold=threshold).sum(axis=0)
    counts = counts.reindex(pred_mat.columns, fill_value=0).astype(int)
    counts.index.name = "kinase"
    return counts.rename("n_targets").sort_values(ascending=False)


def compute_ksea_scores(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    threshold: float = 0.6,
    min_substrates: int = 3,
) -> tuple[pd.DataFrame, pd.Series]:
    """Compute KSEA-style downstream kinase scores.

    Missing phosphosite values are ignored on a per-sample basis. Each sample score
    is the arithmetic mean across the observed substrates only. Kinases whose
    selected substrates contain no observed phosphosite values in any sample are
    omitted from both returned outputs.
    """

    aligned_pred_mat, aligned_matrix = _align_activity_inputs(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
    )
    substrate_mask = _prediction_mask_array(aligned_pred_mat, threshold=threshold)
    substrate_counts = substrate_mask.sum(axis=0)
    candidate_kinase_positions = np.flatnonzero(substrate_counts >= min_substrates)

    if len(candidate_kinase_positions) == 0:
        empty_scores = pd.DataFrame(columns=list(phospho_matrix.columns), dtype=float)
        empty_counts = pd.Series(dtype=int, name="n_substrates")
        empty_counts.index.name = "kinase"
        return empty_scores, empty_counts

    matrix_values = aligned_matrix.to_numpy(dtype=float, copy=False)
    score_rows: list[np.ndarray] = []
    score_index: list[str] = []
    counts: dict[str, int] = {}

    for kinase_position in candidate_kinase_positions:
        selected_row_positions = np.flatnonzero(substrate_mask[:, kinase_position])
        kinase_values = matrix_values[selected_row_positions, :]
        kinase_scores = _nan_aware_mean_array(kinase_values)
        if np.isnan(kinase_scores).all():
            continue

        kinase_name = str(aligned_pred_mat.columns[kinase_position])
        counts[kinase_name] = int(selected_row_positions.size)
        score_index.append(kinase_name)
        score_rows.append(kinase_scores)

    if score_rows:
        score_frame = pd.DataFrame(
            score_rows,
            index=score_index,
            columns=aligned_matrix.columns,
            dtype=float,
        )
    else:
        score_frame = pd.DataFrame(columns=list(phospho_matrix.columns), dtype=float)

    count_series = pd.Series(counts, name="n_substrates").sort_values(ascending=False)
    count_series.index.name = "kinase"
    return score_frame, count_series


def _prediction_mask(pred_mat: pd.DataFrame, threshold: float) -> pd.DataFrame:
    return pred_mat.gt(threshold)


def _prediction_mask_array(pred_mat: pd.DataFrame, threshold: float) -> np.ndarray:
    return pred_mat.to_numpy(copy=False) > threshold


def _align_activity_inputs(
    *,
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    common_sites = pred_mat.index.intersection(phospho_matrix.index)
    return pred_mat.loc[common_sites], phospho_matrix.loc[common_sites]


def _build_site_position_lookup(index: pd.Index) -> dict[object, int]:
    return {site_id: position for position, site_id in enumerate(index)}


def _nan_aware_weighted_average(
    values: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    if values.ndim != 2:
        msg = "values must be a two-dimensional array"
        raise ValueError(msg)
    if weights.ndim != 1:
        msg = "weights must be a one-dimensional array"
        raise ValueError(msg)
    if values.shape[0] != weights.shape[0]:
        msg = "values and weights must align by substrate"
        raise ValueError(msg)

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
    return result


def _nan_aware_mean(values: pd.DataFrame) -> pd.Series:
    return pd.Series(
        _nan_aware_mean_array(values.to_numpy(dtype=float, copy=False)),
        index=values.columns,
        dtype=float,
    )


def _nan_aware_mean_array(values: np.ndarray) -> np.ndarray:
    if values.ndim != 2:
        msg = "values must be a two-dimensional array"
        raise ValueError(msg)

    valid_mask = ~np.isnan(values)
    valid_counts = valid_mask.sum(axis=0)
    sums = np.where(valid_mask, values, 0.0).sum(axis=0)
    result = np.full(values.shape[1], np.nan, dtype=float)
    np.divide(sums, valid_counts, out=result, where=valid_counts > 0)
    return result
