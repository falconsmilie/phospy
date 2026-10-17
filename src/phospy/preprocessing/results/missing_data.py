"""Passive schema metadata for missing-data stage diagnostics."""

from __future__ import annotations

from typing import TypeAlias

JsonPrimitive: TypeAlias = None | str | bool | int | float
JsonValue: TypeAlias = JsonPrimitive | list["JsonValue"] | dict[str, "JsonValue"]

MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1 = 1
MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V2 = 2

V1_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS = frozenset(
    (
        "diagnostics_schema_version",
        "missing_data_policy",
        "imputation_method_id",
        "imputation_method_family",
        "input_missing_cell_count",
        "output_missing_cell_count",
        "imputed_cell_count",
        "affected_row_count",
        "affected_column_count",
        "affected_row_ids",
        "affected_column_ids",
        "imputed_row_ids",
        "imputed_column_ids",
        "dropped_row_ids",
        "imputed_row_count",
        "imputed_column_count",
        "dropped_row_count",
        "random_seed",
        "method_parameters",
        "matrix_scale_requirement",
        "imputation_input_scale",
        "imputation_input_scale_source",
        "imputation_operation_order",
        "stage_order",
        "missingness_mask_hash",
        "imputation_mask_hash",
        "left_censored_assumption",
        "rows_not_imputable",
        "row_medians_used",
        "per_column_distribution_parameters",
        "dropped_rows_above_max_missing_fraction",
        "neighbour_count",
        "distance_metric",
        "knn_no_overlap_policy",
        "knn_no_overlap_policy_version",
        "knn_nearest_neighbour_imputed_cell_count",
        "knn_nearest_neighbour_imputed_row_ids",
        "knn_nearest_neighbour_imputed_column_ids",
        "knn_column_mean_fallback_imputed_cell_count",
        "knn_column_mean_fallback_row_ids",
        "knn_column_mean_fallback_column_ids",
        "knn_nearest_neighbour_imputation_mask_hash",
        "knn_column_mean_fallback_imputation_mask_hash",
        "knn_fully_column_mean_fallback_row_ids",
        "diagnostic_caveat_codes",
    )
)

V2_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS = frozenset(
    (*V1_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS, "group_aware")
)

__all__ = [
    "MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1",
    "MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V2",
    "V1_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS",
    "V2_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS",
]
