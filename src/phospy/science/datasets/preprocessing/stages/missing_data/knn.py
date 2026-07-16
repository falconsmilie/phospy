"""Deterministic KNN missing-data policy implementation.

The public `impute_knn` policy is a PhosPy-owned implementation, not a
delegation to scikit-learn. It preserves the original custom semantics:

- rows above `max_missing_fraction_per_row` are dropped before donor search;
- only `distance="nan_euclidean"` is accepted;
- donors for a missing cell must have an observed value in that cell's column;
- row-to-donor distances use only overlapping observed columns and are scaled
  by `n_columns / observed_overlap_count`;
- tied donors are ordered by `(str(row_id), original_position)`;
- selected donor values are averaged without distance weighting;
- if no donor has any overlapping observed value, the retained-column mean is
  used as the deterministic fallback.

The distance kernel is chunked over target rows with missing values. Guardrails
reject requests whose retained shape or estimated distance work is outside the
documented preprocessing performance budget.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.models import PreprocessingState

from .diagnostics import label_preview
from .models import KnnPolicyOutcome, RowImputationRecord

KNN_DISTANCE_METRIC_NAN_EUCLIDEAN = "nan_euclidean"
KNN_MAX_RETAINED_SITE_COUNT = 50_000
KNN_MAX_SAMPLE_COUNT = 64
KNN_MAX_DISTANCE_FEATURE_OPERATIONS = 2_000_000_000
KNN_DISTANCE_CHUNK_MATRIX_MIB = 96.0
KNN_PEAK_MEMORY_BUDGET_MIB = 384.0


def run_knn_policy(state: PreprocessingState) -> KnnPolicyOutcome:
    """Apply deterministic KNN imputation to the preprocessing state.

    Current semantics:

    - rows whose missing fraction is greater than
      `missing_data.max_missing_fraction_per_row` are dropped and reported as
      not imputable;
    - retained columns must each have at least one observed value;
    - each missing cell is imputed from up to `k` nearest retained donor rows
      that are observed in that cell's column;
    - donor distance is nan-euclidean over the target and donor row's shared
      observed columns;
    - exact distance ties are resolved by stringified row label, then original
      retained-row position;
    - when a missing cell has no donor with any shared observed column, the
      retained-column mean is used.

    Runtime is bounded by rejecting retained requests above
    `KNN_MAX_DISTANCE_FEATURE_OPERATIONS` estimated distance-feature operations
    and by processing accepted requests in target-row chunks capped by
    `KNN_DISTANCE_CHUNK_MATRIX_MIB` per pairwise matrix.
    """

    k = state.plan.missing_data_k
    distance = state.plan.missing_data_distance
    max_missing_fraction_per_row = state.plan.missing_data_max_missing_fraction_per_row
    if k is None or distance is None or max_missing_fraction_per_row is None:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.missing_data.policy="
            "'impute_knn' requires k, distance, and "
            "max_missing_fraction_per_row"
        )
    k_value = int(k)
    distance_value = str(distance).strip()
    max_missing_fraction_value = float(max_missing_fraction_per_row)
    if distance_value != KNN_DISTANCE_METRIC_NAN_EUCLIDEAN:
        raise PhosPyInputError(
            "dataset preprocessing stage 'missing_data' cannot apply "
            "missing_data.policy='impute_knn' because "
            "missing_data.distance must be 'nan_euclidean'."
        )

    missing_fraction = state.phospho.isna().mean(axis=1)
    retained_mask = missing_fraction <= max_missing_fraction_value
    dropped_missing_fraction = missing_fraction.loc[~retained_mask]
    filtered_phospho = state.phospho.loc[retained_mask].copy(deep=True)
    filtered_site_metadata = state.site_metadata.loc[filtered_phospho.index]

    if filtered_phospho.empty:
        imputed = filtered_phospho.copy(deep=True)
    else:
        all_missing_columns = filtered_phospho.columns[
            filtered_phospho.notna().sum(axis=0).to_numpy(dtype=int, copy=False) == 0
        ]
        if len(all_missing_columns) > 0:
            raise PhosPyInputError(
                "dataset preprocessing stage 'missing_data' cannot apply "
                "missing_data.policy='impute_knn' because one or more columns "
                "have no observed values after row filtering. "
                f"affected column labels (preview): {label_preview(all_missing_columns.tolist())}. "
                "adjust missing_data.max_missing_fraction_per_row or input data."
            )
        _validate_knn_scale_request(filtered_phospho)
        imputed_values = _deterministic_knn_impute(
            filtered_phospho,
            n_neighbors=k_value,
        )
        if imputed_values.shape[1] != filtered_phospho.shape[1]:
            raise PhosPyInputError(
                "dataset preprocessing stage 'missing_data' cannot apply "
                "missing_data.policy='impute_knn' because the imputer could not "
                "retain all matrix columns during imputation. "
                "ensure every retained column has at least one observed value."
            )
        imputed = pd.DataFrame(
            imputed_values,
            index=filtered_phospho.index.copy(),
            columns=filtered_phospho.columns.copy(),
        )

    if filtered_phospho.empty:
        imputed_mask = filtered_phospho.isna() & filtered_phospho.notna()
    else:
        imputed_mask = filtered_phospho.isna() & imputed.notna()

    imputed_cell_count = int(imputed_mask.to_numpy().sum())
    imputed_row_ids = (
        tuple(
            str(row_id)
            for row_id in imputed.index[
                imputed_mask.any(axis=1).to_numpy(dtype=bool, copy=False)
            ].tolist()
        )
        if not imputed.empty
        else ()
    )
    imputed_column_ids = (
        tuple(
            str(column_name)
            for column_name in imputed.columns[
                imputed_mask.any(axis=0).to_numpy(dtype=bool, copy=False)
            ].tolist()
        )
        if not imputed.empty
        else ()
    )
    unresolved_mask = imputed.isna() & filtered_phospho.isna()
    unresolved_row_ids = (
        tuple(
            str(row_id)
            for row_id in imputed.index[
                unresolved_mask.any(axis=1).to_numpy(dtype=bool, copy=False)
            ].tolist()
        )
        if not imputed.empty
        else ()
    )
    unresolved_column_ids = (
        tuple(
            str(column_name)
            for column_name in imputed.columns[
                unresolved_mask.any(axis=0).to_numpy(dtype=bool, copy=False)
            ].tolist()
        )
        if not imputed.empty
        else ()
    )
    output_missing_cell_count = int(imputed.isna().to_numpy().sum())
    if output_missing_cell_count > 0:
        raise PhosPyInputError(
            "dataset preprocessing stage 'missing_data' could not complete "
            "missing_data.policy='impute_knn' because missing values remain after "
            "imputation. "
            f"remaining rows (preview): {label_preview(list(unresolved_row_ids))}. "
            f"remaining columns (preview): {label_preview(list(unresolved_column_ids))}. "
            "adjust missing_data.max_missing_fraction_per_row, k, or input data."
        )

    dropped_row_ids = tuple(
        str(row_id) for row_id in dropped_missing_fraction.index.tolist()
    )
    imputed_rows = (
        tuple(
            RowImputationRecord(
                row_id=str(row_id),
                imputed_columns=tuple(
                    str(column_name)
                    for column_name in imputed.columns[
                        imputed_mask.loc[row_id]
                    ].tolist()
                ),
                imputed_cell_count=int(imputed_mask.loc[row_id].sum()),
            )
            for row_id in imputed.index[
                imputed_mask.any(axis=1).to_numpy(dtype=bool, copy=False)
            ]
        )
        if not imputed.empty
        else ()
    )
    dropped_rows_missing_fraction = tuple(
        (str(row_id), float(missing_fraction_value))
        for row_id, missing_fraction_value in dropped_missing_fraction.items()
    )
    return KnnPolicyOutcome(
        phospho=imputed,
        site_metadata=filtered_site_metadata,
        imputed_mask=imputed_mask,
        k=k_value,
        distance=distance_value,
        max_missing_fraction_per_row=max_missing_fraction_value,
        dropped_row_ids=dropped_row_ids,
        dropped_rows_missing_fraction=dropped_rows_missing_fraction,
        imputed_cell_count=imputed_cell_count,
        imputed_row_ids=imputed_row_ids,
        imputed_column_ids=imputed_column_ids,
        output_missing_cell_count=output_missing_cell_count,
        rows_not_imputable=dropped_row_ids,
        imputed_rows=imputed_rows,
    )


def _deterministic_knn_impute(
    phospho: pd.DataFrame,
    *,
    n_neighbors: int,
) -> np.ndarray:
    """Impute with chunked nan-euclidean KNN and deterministic donor ties."""

    if n_neighbors < 1:
        raise PhosPyInputError(
            "dataset preprocessing stage 'missing_data' cannot apply "
            "missing_data.policy='impute_knn' because k must be at least 1."
        )
    values = phospho.to_numpy(dtype=float, copy=True)
    imputed_values = values.copy()
    missing_mask = np.isnan(values)
    if not bool(missing_mask.any()):
        return imputed_values

    column_means = np.nanmean(values, axis=0)
    row_sort_order = np.asarray(
        sorted(
            range(int(values.shape[0])),
            key=lambda position: (str(phospho.index[position]), int(position)),
        ),
        dtype=np.int64,
    )
    sorted_values = values[row_sort_order, :]
    sorted_observed = ~np.isnan(sorted_values)
    sorted_values_filled = np.where(sorted_observed, sorted_values, 0.0)
    sorted_values_squared = sorted_values_filled * sorted_values_filled
    sorted_observed_float = sorted_observed.astype(float, copy=False)
    eligible_donor_positions_by_column = tuple(
        np.flatnonzero(sorted_observed[:, column_position])
        for column_position in range(int(values.shape[1]))
    )

    target_positions = np.flatnonzero(missing_mask.any(axis=1))
    chunk_size = _knn_distance_chunk_size(retained_row_count=int(values.shape[0]))
    for chunk_start in range(0, int(target_positions.size), chunk_size):
        chunk_positions = target_positions[chunk_start : chunk_start + chunk_size]
        distances = _nan_euclidean_distance_block(
            target_values=values[chunk_positions, :],
            sorted_values_filled=sorted_values_filled,
            sorted_values_squared=sorted_values_squared,
            sorted_observed_float=sorted_observed_float,
        )
        for local_position, row_position_raw in enumerate(chunk_positions):
            row_position = int(row_position_raw)
            row_distances = distances[local_position, :]
            for column_position_raw in np.flatnonzero(missing_mask[row_position, :]):
                column_position = int(column_position_raw)
                donor_positions = eligible_donor_positions_by_column[column_position]
                eligible_distances = row_distances[donor_positions]
                finite_mask = np.isfinite(eligible_distances)
                if bool(finite_mask.any()):
                    selected_positions = _select_knn_donor_positions(
                        donor_positions=donor_positions[finite_mask],
                        distances=eligible_distances[finite_mask],
                        n_neighbors=n_neighbors,
                    )
                    imputed_values[row_position, column_position] = float(
                        np.mean(sorted_values[selected_positions, column_position])
                    )
                else:
                    imputed_values[row_position, column_position] = float(
                        column_means[column_position]
                    )

    return imputed_values


def _validate_knn_scale_request(phospho: pd.DataFrame) -> None:
    """Reject retained KNN requests outside the documented execution envelope."""

    missing_target_row_count = int(phospho.isna().any(axis=1).sum())
    if missing_target_row_count == 0:
        return

    retained_row_count = int(phospho.shape[0])
    sample_count = int(phospho.shape[1])
    if retained_row_count > KNN_MAX_RETAINED_SITE_COUNT:
        raise PhosPyInputError(
            "dataset preprocessing stage 'missing_data' cannot apply "
            "missing_data.policy='impute_knn' because the retained matrix has "
            f"{retained_row_count} site rows after missing-fraction filtering; "
            f"the custom deterministic KNN guardrail is {KNN_MAX_RETAINED_SITE_COUNT}. "
            "Reduce the number of retained sites, lower "
            "missing_data.max_missing_fraction_per_row, or choose a simpler "
            "missing-data policy such as 'impute_row_median'."
        )
    if sample_count > KNN_MAX_SAMPLE_COUNT:
        raise PhosPyInputError(
            "dataset preprocessing stage 'missing_data' cannot apply "
            "missing_data.policy='impute_knn' because the retained matrix has "
            f"{sample_count} sample columns; the custom deterministic KNN "
            f"guardrail is {KNN_MAX_SAMPLE_COUNT}. Reduce sample count before "
            "KNN imputation or choose a simpler missing-data policy such as "
            "'impute_row_median'."
        )

    estimated_feature_operations = _estimate_knn_feature_operations(phospho)
    if estimated_feature_operations > KNN_MAX_DISTANCE_FEATURE_OPERATIONS:
        raise PhosPyInputError(
            "dataset preprocessing stage 'missing_data' cannot apply "
            "missing_data.policy='impute_knn' because the retained request would "
            "require an impractical number of deterministic nan-euclidean "
            "distance operations. "
            f"retained_rows={retained_row_count}, "
            f"rows_with_missing_values={missing_target_row_count}, "
            f"sample_columns={sample_count}, "
            f"estimated_distance_feature_operations={estimated_feature_operations}, "
            f"budget={KNN_MAX_DISTANCE_FEATURE_OPERATIONS}. Reduce the number "
            "of retained missing rows, lower "
            "missing_data.max_missing_fraction_per_row, pre-filter low-value "
            "features, or choose 'impute_row_median'."
        )


def _estimate_knn_feature_operations(phospho: pd.DataFrame) -> int:
    missing_target_row_count = int(phospho.isna().any(axis=1).sum())
    return int(missing_target_row_count * int(phospho.shape[0]) * int(phospho.shape[1]))


def _knn_distance_chunk_size(*, retained_row_count: int) -> int:
    bytes_per_pairwise_matrix = max(
        1.0,
        float(KNN_DISTANCE_CHUNK_MATRIX_MIB) * 1024.0 * 1024.0,
    )
    bytes_per_target_row = max(1, int(retained_row_count)) * 8.0
    return max(1, int(bytes_per_pairwise_matrix // bytes_per_target_row))


def _nan_euclidean_distance_block(
    *,
    target_values: np.ndarray,
    sorted_values_filled: np.ndarray,
    sorted_values_squared: np.ndarray,
    sorted_observed_float: np.ndarray,
) -> np.ndarray:
    """Return target-by-donor nan-euclidean distances in donor tie order."""

    target_observed = ~np.isnan(target_values)
    target_observed_float = target_observed.astype(float, copy=False)
    target_values_filled = np.where(target_observed, target_values, 0.0)
    target_values_squared = target_values_filled * target_values_filled

    overlap_counts = target_observed_float @ sorted_observed_float.T
    squared_distances = target_values_squared @ sorted_observed_float.T
    squared_distances += target_observed_float @ sorted_values_squared.T
    squared_distances -= 2.0 * (target_values_filled @ sorted_values_filled.T)
    np.maximum(squared_distances, 0.0, out=squared_distances)

    with np.errstate(divide="ignore", invalid="ignore"):
        squared_distances *= float(target_values.shape[1])
        squared_distances /= overlap_counts
    squared_distances[overlap_counts <= 0.0] = np.inf
    return np.sqrt(squared_distances, out=squared_distances)


def _select_knn_donor_positions(
    *,
    donor_positions: np.ndarray,
    distances: np.ndarray,
    n_neighbors: int,
) -> np.ndarray:
    """Select donor positions from arrays already ordered by deterministic tie key."""

    if n_neighbors == 1:
        return donor_positions[[int(np.argmin(distances))]]
    if int(donor_positions.size) <= int(n_neighbors):
        selected_order = np.argsort(distances, kind="mergesort")
        return donor_positions[selected_order]
    selected_order = np.argsort(distances, kind="mergesort")[:n_neighbors]
    return donor_positions[selected_order]
