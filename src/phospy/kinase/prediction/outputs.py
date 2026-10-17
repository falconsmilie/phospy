"""Ranking and output construction for Kinase prediction."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.kinase.scoring.profile import _column_series


def rank_kinases_for_prediction(
    *,
    prediction_score_matrix: pd.DataFrame,
    candidate_substrates: dict[str, list[str]],
) -> pd.Series:
    """Rank kinases by the mean candidate-site score used for prediction."""

    ranking: dict[str, float] = {}
    for kinase, candidate_sites in candidate_substrates.items():
        if kinase not in prediction_score_matrix.columns:
            continue
        available_sites = [
            site for site in candidate_sites if site in prediction_score_matrix.index
        ]
        if not available_sites:
            continue
        values = (
            prediction_score_matrix.loc[available_sites, kinase].astype(float).dropna()
        )
        if values.empty:
            continue
        ranking[kinase] = float(values.mean(skipna=False))
    ranking_series = pd.Series(ranking, dtype=float, name="prediction_rank_score")
    return ranking_series.sort_values(ascending=False, kind="mergesort")


def build_prediction_outputs(
    *,
    prediction_score_matrix: pd.DataFrame,
    selected_kinases: pd.Index,
    candidate_substrates: dict[str, list[str]],
    top_k: int,
    retain_full_scores: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build the prediction matrix and substrate table for selected kinases."""

    pred_mat = pd.DataFrame(
        np.nan,
        index=prediction_score_matrix.index.copy(),
        columns=selected_kinases.copy(),
    )
    pred_mat.index.name = prediction_score_matrix.index.name
    pred_mat.columns.name = "kinase"

    score_index = prediction_score_matrix.index
    score_values = prediction_score_matrix.to_numpy(dtype=float, copy=False)
    site_ids = score_index.to_numpy(copy=False)
    use_vectorized_site_indexer = score_index.is_unique

    substrate_rows: list[dict[str, object]] = []
    for kinase_position, kinase in enumerate(selected_kinases):
        candidate_sites = candidate_substrates.get(str(kinase), [])
        if not candidate_sites:
            continue
        if retain_full_scores:
            full_scores = _column_series(prediction_score_matrix, str(kinase)).astype(
                float
            )
            pred_mat.iloc[:, kinase_position] = full_scores.to_numpy(
                dtype=float,
                copy=False,
            )
            available_sites = [site for site in candidate_sites if site in score_index]
            if not available_sites:
                continue
            ranked_sites = (
                full_scores.reindex(available_sites)
                .dropna()
                .nlargest(
                    top_k,
                    keep="first",
                )
            )
            for rank, (site_id, score) in enumerate(ranked_sites.items(), start=1):
                substrate_rows.append(
                    {
                        "kinase": str(kinase),
                        "substrate_site": site_id,
                        "score": float(score),
                        "rank": rank,
                    }
                )
            continue
        score_column_position = prediction_score_matrix.columns.get_loc(kinase)

        if use_vectorized_site_indexer:
            site_positions = score_index.get_indexer(pd.Index(candidate_sites))
            available_positions = site_positions[site_positions >= 0]
            if available_positions.size == 0:
                continue
            candidate_scores = score_values[available_positions, score_column_position]
            scored_mask = ~np.isnan(candidate_scores)
            if not scored_mask.any():
                continue
            available_positions = available_positions[scored_mask]
            candidate_scores = candidate_scores[scored_mask]
            retained_count = min(int(top_k), int(candidate_scores.size))
            if retained_count <= 0:
                continue
            if retained_count < int(candidate_scores.size):
                retained_positions = np.argpartition(
                    -candidate_scores,
                    retained_count - 1,
                )[:retained_count]
                selected_positions = available_positions[retained_positions]
                selected_scores = candidate_scores[retained_positions]
            else:
                selected_positions = available_positions
                selected_scores = candidate_scores

            ranked_positions = np.lexsort((selected_positions, -selected_scores))
            selected_positions = selected_positions[ranked_positions]
            selected_scores = selected_scores[ranked_positions]
            pred_mat.iloc[selected_positions, kinase_position] = selected_scores
            for rank, (site_position, score) in enumerate(
                zip(selected_positions, selected_scores, strict=True),
                start=1,
            ):
                substrate_rows.append(
                    {
                        "kinase": str(kinase),
                        "substrate_site": site_ids[site_position],
                        "score": float(score),
                        "rank": rank,
                    }
                )
            continue

        available_sites = [site for site in candidate_sites if site in score_index]
        if not available_sites:
            continue
        ranked_sites = (
            _column_series(prediction_score_matrix, str(kinase))
            .reindex(available_sites)
            .astype(float)
            .dropna()
            .nlargest(top_k, keep="first")
        )
        if ranked_sites.empty:
            continue
        pred_mat.loc[ranked_sites.index, kinase] = ranked_sites.to_numpy(
            dtype=float, copy=False
        )
        for rank, (site_id, score) in enumerate(ranked_sites.items(), start=1):
            substrate_rows.append(
                {
                    "kinase": str(kinase),
                    "substrate_site": site_id,
                    "score": float(score),
                    "rank": rank,
                }
            )

    substrate_list = pd.DataFrame(
        substrate_rows,
        columns=pd.Index(["kinase", "substrate_site", "score", "rank"]),
    )
    return pred_mat, substrate_list


__all__ = [
    "build_prediction_outputs",
    "rank_kinases_for_prediction",
]
