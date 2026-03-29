from __future__ import annotations

import numpy as np
import pandas as pd


def compute_weighted_kinase_activity(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    top_n_substrates: int = 20,
    min_substrates: int = 3,
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
        values = phospho_matrix.loc[substrates, samples].to_numpy(dtype=float)
        weighted_values = np.average(values, axis=0, weights=weights)
        kinase_mat.loc[kinase, :] = weighted_values

    return kinase_mat.dropna(how="all")


def build_kinase_target_table(
    pred_mat: pd.DataFrame,
    threshold: float = 0.6,
) -> pd.DataFrame:
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
    edges = build_kinase_target_table(pred_mat, threshold=threshold)
    counts = edges.groupby("kinase").size()
    counts = counts.reindex(pred_mat.columns, fill_value=0)
    counts.index.name = "kinase"
    return counts.rename("n_targets").sort_values(ascending=False)


def compute_ksea_scores(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    threshold: float = 0.6,
    min_substrates: int = 3,
) -> tuple[pd.DataFrame, pd.Series]:
    edges = build_kinase_target_table(pred_mat, threshold=threshold)

    substrate_map = {
        kinase: group["site_id"].tolist() for kinase, group in edges.groupby("kinase")
    }
    substrate_map = {k: v for k, v in substrate_map.items() if len(v) >= min_substrates}

    if not substrate_map:
        empty_scores = pd.DataFrame(columns=list(phospho_matrix.columns), dtype=float)
        empty_counts = pd.Series(dtype=int, name="n_substrates")
        empty_counts.index.name = "kinase"
        return empty_scores, empty_counts

    score_dict: dict[str, pd.Series] = {}
    counts: dict[str, int] = {}

    available_sites = set(phospho_matrix.index)

    for kinase, sites in substrate_map.items():
        present_sites = [site for site in sites if site in available_sites]
        counts[kinase] = len(present_sites)
        if len(present_sites) < min_substrates:
            continue
        score_dict[kinase] = phospho_matrix.loc[present_sites].mean(axis=0)

    score_frame = pd.DataFrame.from_dict(score_dict, orient="index")
    if score_frame.empty:
        score_frame = pd.DataFrame(columns=list(phospho_matrix.columns), dtype=float)
    count_series = pd.Series(counts, name="n_substrates").sort_values(ascending=False)
    count_series.index.name = "kinase"
    return score_frame, count_series
