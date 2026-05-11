"""Simplified weighted substrate activity method implementation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.activities.models import (
    SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD,
    ActivityMethodSummary,
    KinaseActivityInputs,
    KinaseActivityResult,
)
from phospy.activities.threshold_membership import (
    build_activity_threshold_membership_diagnostics,
    threshold_membership_filtered_frame,
    threshold_membership_mask_array,
    threshold_membership_mask_frame,
)
from phospy.errors.workflows import WorkflowBoundaryError

_SITE_ID_COLUMN = "site_id"


@dataclass(frozen=True, slots=True)
class SimplifiedWeightedSubstrateActivityMethod:
    """Heuristic weighted activity method."""

    threshold: float
    min_substrates: int
    top_n_substrates: int

    def run(self, inputs: KinaseActivityInputs) -> KinaseActivityResult:
        weighted_activity = _compute_weighted_kinase_activity(
            pred_mat=inputs.pred_mat,
            phospho_matrix=inputs.phospho_matrix,
            top_n_substrates=self.top_n_substrates,
            min_substrates=self.min_substrates,
        )
        threshold_diagnostics = build_activity_threshold_membership_diagnostics(
            threshold_parameter="threshold",
            threshold_value=float(self.threshold),
        )
        (
            thresholded_substrate_mean_activity,
            thresholded_substrate_counts,
        ) = _compute_thresholded_substrate_mean_activity(
            pred_mat=inputs.pred_mat,
            phospho_matrix=inputs.phospho_matrix,
            threshold=self.threshold,
            min_substrates=self.min_substrates,
        )
        target_counts = _count_predicted_targets(
            pred_mat=inputs.pred_mat,
            threshold=self.threshold,
        )
        target_table = _build_kinase_target_table(
            pred_mat=inputs.pred_mat,
            threshold=self.threshold,
        )

        if weighted_activity.empty and thresholded_substrate_mean_activity.empty:
            _raise_boundary_error(
                seam="kinase.activity.valid_candidates",
                next_action=(
                    "lower activity_config.min_substrates, increase "
                    "activity_config.top_n_substrates, or lower "
                    "activity_config.threshold to retain kinase activity candidates"
                ),
                overlap_sites=int(inputs.overlap_summary.overlap_count),
                pred_mat_sites=int(inputs.overlap_summary.pred_mat_rows),
                phospho_sites=int(inputs.overlap_summary.phospho_rows),
                candidate_kinases=int(inputs.pred_mat.shape[1]),
                weighted_activity_kinases=0,
                thresholded_mean_activity_kinases=0,
                activity_config_threshold=float(self.threshold),
                activity_config_min_substrates=int(self.min_substrates),
                activity_config_top_n_substrates=int(self.top_n_substrates),
            )

        evaluated_pairs = int(inputs.pred_mat.shape[1] * inputs.phospho_matrix.shape[1])
        summary = ActivityMethodSummary(
            kinases_evaluated=int(inputs.pred_mat.shape[1]),
            kinase_condition_pairs_evaluated=evaluated_pairs,
            kinase_condition_pairs_computed=evaluated_pairs,
            kinase_condition_pairs_insufficient_substrates=0,
            kinase_condition_pairs_invalid_background_variance=0,
            kinase_condition_pairs_no_finite_background_values=0,
            kinase_condition_pairs_no_finite_substrate_values=0,
        )
        return KinaseActivityResult._from_owned(
            weighted_activity=weighted_activity,
            thresholded_substrate_mean_activity=thresholded_substrate_mean_activity,
            thresholded_substrate_counts=thresholded_substrate_counts,
            target_counts=target_counts,
            target_table=target_table,
            threshold_membership_diagnostics=threshold_diagnostics,
            statistics_table=None,
            method_summary=summary,
            activity_method=SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY_METHOD,
        )


def _compute_weighted_kinase_activity(
    *,
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    top_n_substrates: int,
    min_substrates: int,
) -> pd.DataFrame:
    kinases = pred_mat.columns.tolist()
    samples = phospho_matrix.columns.tolist()
    kinase_rows: list[np.ndarray] = []
    kinase_names: list[str] = []

    phospho_values = phospho_matrix.to_numpy(dtype=float, copy=False)
    phospho_site_positions = _build_site_position_lookup(phospho_matrix.index)
    pred_values = pred_mat.to_numpy(dtype=float, copy=False)
    pred_site_positions = np.asarray(
        [phospho_site_positions.get(site_id, -1) for site_id in pred_mat.index],
        dtype=int,
    )
    valid_prediction_mask = ~np.isnan(pred_values)
    non_nan_counts = valid_prediction_mask.sum(axis=0)
    sortable_pred_values = np.where(valid_prediction_mask, pred_values, -np.inf)
    sorted_site_positions = np.argsort(-sortable_pred_values, axis=0, kind="stable")

    for kinase_position, kinase_name in enumerate(kinases):
        top_count = min(top_n_substrates, int(non_nan_counts[kinase_position]))
        if top_count == 0:
            continue

        top_pred_positions = sorted_site_positions[:top_count, kinase_position]
        aligned_site_positions = pred_site_positions[top_pred_positions]
        aligned_mask = aligned_site_positions >= 0
        substrate_positions = aligned_site_positions[aligned_mask]

        if substrate_positions.size < min_substrates:
            continue

        weights_array = pred_values[top_pred_positions[aligned_mask], kinase_position]
        if float(weights_array.sum()) <= 0.0:
            continue

        substrate_values = phospho_values[substrate_positions, :]
        weighted_values = _nan_aware_weighted_average(substrate_values, weights_array)
        if np.isnan(weighted_values).all():
            continue

        kinase_names.append(str(kinase_name))
        kinase_rows.append(weighted_values)

    if not kinase_rows:
        return pd.DataFrame(columns=samples, dtype=float)

    result = pd.DataFrame(kinase_rows, index=kinase_names, columns=samples, dtype=float)
    result.index.name = "kinase"
    return result


def _compute_thresholded_substrate_mean_activity(
    *,
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    threshold: float,
    min_substrates: int,
) -> tuple[pd.DataFrame, pd.Series]:
    aligned_pred_mat, aligned_matrix = _align_activity_inputs(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
    )
    substrate_mask = _prediction_mask_array(aligned_pred_mat, threshold=threshold)
    substrate_counts = substrate_mask.sum(axis=0)
    candidate_kinase_positions = np.flatnonzero(substrate_counts >= min_substrates)

    if len(candidate_kinase_positions) == 0:
        empty_scores = pd.DataFrame(
            columns=list(phospho_matrix.columns),
            dtype=float,
        )
        empty_counts = pd.Series(dtype=int, name="n_substrates")
        empty_counts.index.name = "kinase"
        return empty_scores, empty_counts

    matrix_values = aligned_matrix.to_numpy(dtype=float, copy=False)
    score_rows: list[np.ndarray] = []
    score_index: list[str] = []
    counts: dict[str, int] = {}

    for kinase_position in candidate_kinase_positions:
        selected_row_positions = np.flatnonzero(substrate_mask[:, int(kinase_position)])
        kinase_values = matrix_values[selected_row_positions, :]
        kinase_scores = _nan_aware_mean_array(kinase_values)
        if np.isnan(kinase_scores).all():
            continue

        kinase_name = str(aligned_pred_mat.columns[int(kinase_position)])
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
        score_frame = pd.DataFrame(
            columns=list(phospho_matrix.columns),
            dtype=float,
        )

    score_frame.index.name = "kinase"
    count_series = pd.Series(counts, name="n_substrates").sort_values(ascending=False)
    count_series.index.name = "kinase"
    return score_frame, count_series


def _count_predicted_targets(
    *,
    pred_mat: pd.DataFrame,
    threshold: float,
) -> pd.Series:
    counts = _prediction_mask(pred_mat, threshold=threshold).sum(axis=0)
    counts = counts.reindex(pred_mat.columns, fill_value=0).astype(int)
    counts.index.name = "kinase"
    return counts.rename("n_targets").sort_values(ascending=False)


def _build_kinase_target_table(
    *,
    pred_mat: pd.DataFrame,
    threshold: float,
) -> pd.DataFrame:
    filtered = threshold_membership_filtered_frame(pred_mat, threshold=threshold)
    try:
        edges = filtered.stack(future_stack=True).rename("score").reset_index()
    except TypeError:
        edges = filtered.stack().rename("score").reset_index()
    edges = edges.loc[edges["score"].notna()]
    edges.columns = [_SITE_ID_COLUMN, "kinase", "score"]
    return edges.sort_values(["kinase", "score"], ascending=[True, False])


def _prediction_mask(pred_mat: pd.DataFrame, threshold: float) -> pd.DataFrame:
    return threshold_membership_mask_frame(pred_mat, threshold=threshold)


def _prediction_mask_array(pred_mat: pd.DataFrame, threshold: float) -> np.ndarray:
    return threshold_membership_mask_array(
        pred_mat.to_numpy(dtype=float, copy=False),
        threshold=threshold,
    )


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
        raise ValueError("values must be a two-dimensional array")
    if weights.ndim != 1:
        raise ValueError("weights must be a one-dimensional array")
    if values.shape[0] != weights.shape[0]:
        raise ValueError("values and weights must align by substrate")

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


def _nan_aware_mean_array(values: np.ndarray) -> np.ndarray:
    if values.ndim != 2:
        raise ValueError("values must be a two-dimensional array")
    valid_mask = ~np.isnan(values)
    valid_counts = valid_mask.sum(axis=0)
    sums = np.where(valid_mask, values, 0.0).sum(axis=0)
    result = np.full(values.shape[1], np.nan, dtype=float)
    np.divide(sums, valid_counts, out=result, where=valid_counts > 0)
    return result


def _raise_boundary_error(
    *,
    seam: str,
    next_action: str,
    **details: int | float,
) -> None:
    details_text = ", ".join(f"{key}={value}" for key, value in details.items())
    raise WorkflowBoundaryError(
        "kinase workflow boundary validation failed at "
        f"seam={seam}; {details_text}; next_action={next_action}"
    )


__all__ = ["SimplifiedWeightedSubstrateActivityMethod"]
