"""Scientific kernels for the supported kinase workflow route."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.api.configs import KINASE_SCORING_MIN_SUBSTRATES_FLOOR
from phospy.errors.workflows import WorkflowStageError


@dataclass(frozen=True, slots=True)
class KinaseProfileBuild:
    """Resolved profile inputs for kinase scoring and downstream prediction."""

    profile_matrix: pd.DataFrame
    quantified_substrates: dict[str, list[str]]
    substrate_counts: pd.Series


def build_kinase_profiles(
    *,
    phospho: pd.DataFrame,
    kinase_substrate_map: pd.DataFrame,
    min_substrates: int,
) -> KinaseProfileBuild:
    """Build kinase substrate profiles from quantified substrate rows."""

    if min_substrates < KINASE_SCORING_MIN_SUBSTRATES_FLOOR:
        raise WorkflowStageError(
            "kinase workflow internal invariant failed at seam="
            "kinase.science.min_substrate_floor; "
            f"min_substrates must be greater than or equal to "
            f"{KINASE_SCORING_MIN_SUBSTRATES_FLOOR}"
        )

    numeric_phospho = _require_finite_matrix(
        phospho,
        field_name="kinase.workflow.dataset.phospho",
    )
    observed_sites = set(numeric_phospho.index)

    profile_rows: dict[str, pd.Series] = {}
    quantified_substrates: dict[str, list[str]] = {}
    substrate_counts: dict[str, int] = {}

    for kinase, grouped in kinase_substrate_map.groupby("kinase", sort=False):
        unique_sites = list(dict.fromkeys(grouped.loc[:, "substrate_site"].tolist()))
        quantified_sites = [site for site in unique_sites if site in observed_sites]
        substrate_counts[str(kinase)] = len(quantified_sites)
        if len(quantified_sites) < min_substrates:
            continue
        quantified_matrix = numeric_phospho.loc[quantified_sites, :]
        if quantified_matrix.shape[0] == 1:
            profile = quantified_matrix.iloc[0].astype(float)
        else:
            profile = quantified_matrix.median(axis=0, skipna=False).astype(float)
        profile_rows[str(kinase)] = profile
        quantified_substrates[str(kinase)] = quantified_sites

    if profile_rows:
        profile_matrix = pd.DataFrame.from_dict(profile_rows, orient="index")
        profile_matrix = profile_matrix.loc[:, phospho.columns.copy()]
    else:
        profile_matrix = pd.DataFrame(columns=phospho.columns.copy(), dtype=float)

    profile_matrix.index.name = "kinase"
    substrate_counts_series = pd.Series(substrate_counts, dtype="int64", name="NumSub")
    substrate_counts_series.index.name = "kinase"

    return KinaseProfileBuild(
        profile_matrix=profile_matrix,
        quantified_substrates=quantified_substrates,
        substrate_counts=substrate_counts_series,
    )


def score_profile_correlations(
    *,
    phospho: pd.DataFrame,
    profile_matrix: pd.DataFrame,
) -> pd.DataFrame:
    """Score each site against kinase profiles using row-wise Pearson correlation."""

    if profile_matrix.empty:
        return pd.DataFrame(
            index=phospho.index.copy(),
            columns=profile_matrix.index.copy(),
            dtype=float,
        )
    if set(phospho.columns) != set(profile_matrix.columns):
        raise WorkflowStageError(
            "phospho sample columns must match kinase profile columns"
        )
    aligned_phospho = _require_finite_matrix(
        phospho.loc[:, profile_matrix.columns],
        field_name="kinase.workflow.scoring_phospho",
    )
    profile_matrix = _require_finite_matrix(
        profile_matrix,
        field_name="kinase.workflow.profile_matrix",
    )
    left = aligned_phospho.to_numpy(dtype=float)
    right = profile_matrix.to_numpy(dtype=float)

    left_centered = left - left.mean(axis=1, keepdims=True)
    right_centered = right - right.mean(axis=1, keepdims=True)
    left_scale = np.linalg.norm(left_centered, axis=1)
    right_scale = np.linalg.norm(right_centered, axis=1)
    denominator = np.outer(left_scale, right_scale)
    with np.errstate(divide="ignore", invalid="ignore"):
        correlation = left_centered @ right_centered.T / denominator
    correlation[denominator == 0.0] = np.nan

    scores = (correlation + 1.0) / 2.0
    valid = np.isfinite(scores)
    scores[valid] = np.clip(scores[valid], 0.0, 1.0)
    return pd.DataFrame(
        scores,
        index=phospho.index.copy(),
        columns=profile_matrix.index.copy(),
    )


def rank_kinases_for_prediction(
    *,
    score_matrix: pd.DataFrame,
    quantified_substrates: dict[str, list[str]],
) -> pd.Series:
    """Rank kinases by the mean candidate-site score used for prediction."""

    ranking: dict[str, float] = {}
    for kinase, candidate_sites in quantified_substrates.items():
        if kinase not in score_matrix.columns:
            continue
        available_sites = [
            site for site in candidate_sites if site in score_matrix.index
        ]
        if not available_sites:
            continue
        values = score_matrix.loc[available_sites, kinase].astype(float).dropna()
        if values.empty:
            continue
        ranking[kinase] = float(values.mean(skipna=False))
    ranking_series = pd.Series(ranking, dtype=float, name="prediction_rank_score")
    return ranking_series.sort_values(ascending=False)


def build_prediction_outputs(
    *,
    score_matrix: pd.DataFrame,
    selected_kinases: pd.Index,
    quantified_substrates: dict[str, list[str]],
    top_k: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build the prediction matrix and substrate table for selected kinases."""

    pred_mat = pd.DataFrame(
        0.0,
        index=score_matrix.index.copy(),
        columns=selected_kinases.copy(),
    )
    pred_mat.index.name = score_matrix.index.name
    pred_mat.columns.name = "kinase"

    substrate_rows: list[dict[str, object]] = []
    for kinase in selected_kinases:
        candidate_sites = quantified_substrates.get(str(kinase), [])
        available_sites = [
            site for site in candidate_sites if site in score_matrix.index
        ]
        if not available_sites:
            continue
        ranked_sites = (
            score_matrix.loc[available_sites, kinase]
            .astype(float)
            .dropna()
            .sort_values(ascending=False)
            .head(top_k)
        )
        if ranked_sites.empty:
            continue
        pred_mat.loc[ranked_sites.index, kinase] = ranked_sites.values
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
        columns=["kinase", "substrate_site", "score", "rank"],
    )
    return pred_mat, substrate_list


def _require_finite_matrix(
    value: pd.DataFrame,
    *,
    field_name: str,
) -> pd.DataFrame:
    numeric = value.astype(float)
    if not np.isfinite(numeric.to_numpy(dtype=float, copy=False)).all():
        raise WorkflowStageError(
            "kinase workflow internal invariant failed at "
            "seam=kinase.science.input_finite_values; "
            f"{field_name} must contain finite numeric values"
        )
    return numeric
