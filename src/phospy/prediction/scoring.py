"""Scoring helpers shared by kinase workflow execution."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd


def combine_profile_and_motif_scores(
    *,
    motif_scores: pd.DataFrame,
    profile_scores: pd.DataFrame,
    motif_sizes: pd.Series,
    profile_sizes: pd.Series,
    allow_profile_only_fallback: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Combine motif and profile scores using rank-derived kinase weights."""

    profile_kinases = set(profile_scores.columns)
    overlap = [kinase for kinase in motif_scores.columns if kinase in profile_kinases]

    if not overlap:
        if not allow_profile_only_fallback:
            raise ValueError(
                "no overlapping kinases between motif_scores and profile_scores"
            )
        weights = _profile_only_weights(profile_scores.columns)
        return profile_scores.copy(), weights

    _require_index_members(motif_sizes, overlap, name="motif_sizes")
    _require_index_members(profile_sizes, overlap, name="profile_sizes")

    motif_rank_weight = np.log(motif_sizes.loc[overlap].rank(method="average") + 1.0)
    profile_rank_weight = np.log(
        profile_sizes.loc[overlap].rank(method="average") + 1.0
    )
    total_weight = motif_rank_weight + profile_rank_weight
    combined_scores = (
        motif_scores.loc[:, overlap].multiply(motif_rank_weight, axis=1)
        + profile_scores.loc[:, overlap].multiply(profile_rank_weight, axis=1)
    ).divide(total_weight, axis=1)
    # Motif columns can become all-NaN after min-max scaling when no sequence
    # contrast exists. Preserve profile evidence in that case instead of
    # discarding the kinase column for downstream prediction.
    motif_missing_profile_present = (
        motif_scores.loc[:, overlap].isna() & profile_scores.loc[:, overlap].notna()
    )
    combined_scores = combined_scores.where(
        ~motif_missing_profile_present,
        profile_scores.loc[:, overlap],
    )
    weights = pd.DataFrame(
        {
            "motif_weight": motif_rank_weight / total_weight,
            "profile_weight": profile_rank_weight / total_weight,
            "motif_rank_weight": motif_rank_weight,
            "profile_rank_weight": profile_rank_weight,
        },
        index=overlap,
    )

    if allow_profile_only_fallback:
        profile_only = [
            kinase for kinase in profile_scores.columns if kinase not in set(overlap)
        ]
        if profile_only:
            combined_scores = pd.concat(
                [combined_scores, profile_scores.loc[:, profile_only]],
                axis=1,
            )
            weights = pd.concat(
                [weights, _profile_only_weights(profile_only)],
                axis=0,
            )

    weights.index.name = "kinase"
    return combined_scores, weights


def _require_index_members(
    series: pd.Series,
    members: Sequence[str],
    *,
    name: str,
) -> None:
    missing = [member for member in members if member not in series.index]
    if missing:
        raise ValueError(f"{name} is missing entries for: {', '.join(missing)}")


def _profile_only_weights(kinases: Sequence[str]) -> pd.DataFrame:
    weights = pd.DataFrame(
        {
            "motif_weight": 0.0,
            "profile_weight": 1.0,
            "motif_rank_weight": 0.0,
            "profile_rank_weight": 1.0,
        },
        index=list(kinases),
    )
    weights.index.name = "kinase"
    return weights


__all__ = ["combine_profile_and_motif_scores"]
