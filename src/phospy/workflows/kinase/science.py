"""Scientific kernels for the supported kinase workflow route."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.api.configs import (
    KINASE_PROFILE_MISSING_VALUE_STRATEGIES,
    KINASE_PROFILE_MISSING_VALUE_STRATEGY_MEDIAN_SKIPNA,
    KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT,
    KINASE_SCORING_MIN_SUBSTRATES_FLOOR,
    KinaseProfileMissingValueStrategy,
)
from phospy.errors.workflows import WorkflowStageError
from phospy.workflows.kinase.scoring_transforms import (
    shift_correlation_to_unit_support,
)


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
    allow_single_substrate_profiles: bool = False,
    profile_missing_value_strategy: KinaseProfileMissingValueStrategy = (
        KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT
    ),
) -> KinaseProfileBuild:
    """Build kinase substrate profiles from quantified substrate rows."""

    required_floor = (
        1 if allow_single_substrate_profiles else KINASE_SCORING_MIN_SUBSTRATES_FLOOR
    )
    if min_substrates < required_floor:
        raise WorkflowStageError(
            "kinase workflow internal invariant failed at seam="
            "kinase.science.min_substrate_floor; "
            f"min_substrates must be greater than or equal to "
            f"{required_floor}"
        )

    if profile_missing_value_strategy not in KINASE_PROFILE_MISSING_VALUE_STRATEGIES:
        allowed = ", ".join(sorted(KINASE_PROFILE_MISSING_VALUE_STRATEGIES))
        raise WorkflowStageError(
            "kinase workflow internal invariant failed at seam="
            "kinase.science.profile_missing_value_strategy; "
            "profile_missing_value_strategy must be one of: "
            f"{allowed}"
        )

    numeric_phospho = _require_numeric_matrix(
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
        elif (
            profile_missing_value_strategy
            == KINASE_PROFILE_MISSING_VALUE_STRATEGY_MEDIAN_SKIPNA
        ):
            profile = quantified_matrix.median(axis=0, skipna=True).astype(float)
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
    """Score each site against kinase profiles using row-wise Pearson correlation.

    Correlations with zero-variance denominators are represented as ``NaN`` so
    downstream stages can distinguish unsupported from low-confidence scores.
    """

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
    aligned_phospho = _require_numeric_matrix(
        phospho.loc[:, profile_matrix.columns],
        field_name="kinase.workflow.scoring_phospho",
    )
    profile_matrix = _require_numeric_matrix(
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

    scores = shift_correlation_to_unit_support(correlation)
    return pd.DataFrame(
        scores,
        index=phospho.index.copy(),
        columns=profile_matrix.index.copy(),
    )


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
    return ranking_series.sort_values(ascending=False)


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
            full_scores = prediction_score_matrix.loc[:, kinase].astype(float)
            pred_mat.iloc[:, kinase_position] = full_scores.to_numpy(
                dtype=float,
                copy=False,
            )
            available_sites = [site for site in candidate_sites if site in score_index]
            if not available_sites:
                continue
            ranked_sites = (
                full_scores.loc[available_sites]
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
            site_positions = score_index.get_indexer(candidate_sites)
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
            prediction_score_matrix.loc[available_sites, kinase]
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


def _require_numeric_matrix(
    value: pd.DataFrame,
    *,
    field_name: str,
) -> pd.DataFrame:
    boolean_columns = [
        str(column)
        for column in value.columns
        if pd.api.types.is_bool_dtype(value[column])
    ]
    if boolean_columns:
        raise WorkflowStageError(
            "kinase workflow internal invariant failed at "
            "seam=kinase.science.input_non_boolean_values; "
            f"{field_name} contains boolean columns: {', '.join(boolean_columns)}"
        )
    non_numeric_columns = [
        str(column)
        for column in value.columns
        if not pd.api.types.is_numeric_dtype(value[column])
    ]
    if non_numeric_columns:
        raise WorkflowStageError(
            "kinase workflow internal invariant failed at "
            "seam=kinase.science.input_numeric_values; "
            f"{field_name} contains non-numeric columns: "
            f"{', '.join(non_numeric_columns)}"
        )
    numeric = value.astype(float)
    if np.isinf(numeric.to_numpy(dtype=float, copy=False)).any():
        raise WorkflowStageError(
            "kinase workflow internal invariant failed at "
            "seam=kinase.science.input_non_infinite_values; "
            f"{field_name} must not contain infinite numeric values"
        )
    return numeric
