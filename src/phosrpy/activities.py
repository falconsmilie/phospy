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

    for kinase in kinases:
        top_substrates = pred_mat[kinase].sort_values(ascending=False).head(top_n_substrates)
        substrates = [site for site in top_substrates.index if site in phospho_matrix.index]
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
    edges = (
        pred_mat.reset_index(names="site_id")
        .melt(id_vars="site_id", var_name="kinase", value_name="score")
    )
    return edges.loc[edges["score"] > threshold].sort_values(["kinase", "score"], ascending=[True, False])


def count_predicted_targets(
    pred_mat: pd.DataFrame,
    threshold: float = 0.6,
) -> pd.Series:
    edges = build_kinase_target_table(pred_mat, threshold=threshold)
    return edges.groupby("kinase").size().rename("n_targets").sort_values(ascending=False)


def compute_ksea_scores(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    threshold: float = 0.6,
    min_substrates: int = 3,
) -> tuple[pd.DataFrame, pd.Series]:
    edges = build_kinase_target_table(pred_mat, threshold=threshold)

    substrate_map = {
        kinase: group["site_id"].tolist()
        for kinase, group in edges.groupby("kinase")
    }
    substrate_map = {k: v for k, v in substrate_map.items() if len(v) >= min_substrates}

    score_dict: dict[str, pd.Series] = {}
    counts: dict[str, int] = {}

    for kinase, sites in substrate_map.items():
        present_sites = [site for site in sites if site in phospho_matrix.index]
        counts[kinase] = len(present_sites)
        if len(present_sites) < min_substrates:
            continue
        score_dict[kinase] = phospho_matrix.loc[present_sites].mean(axis=0)

    score_frame = pd.DataFrame.from_dict(score_dict, orient="index")
    count_series = pd.Series(counts, name="n_substrates").sort_values(ascending=False)
    return score_frame, count_series
