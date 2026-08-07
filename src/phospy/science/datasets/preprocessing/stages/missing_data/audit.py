"""Row-audit construction for missing-data stage."""

from __future__ import annotations

import pandas as pd

from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    PreprocessingPlan,
)
from phospy.science.datasets.preprocessing.policy_models import MissingDataPolicy
from phospy.science.datasets.preprocessing.report_schema import PreprocessingRowAuditRow
from phospy.science.datasets.processing_state import JsonValue

from .models import (
    KnnPolicyOutcome,
    MinProbPolicyOutcome,
    MissingDataInputProfile,
    RowMedianPolicyOutcome,
)


def build_row_median_audit_records(
    *,
    plan: PreprocessingPlan,
    input_profile: MissingDataInputProfile,
    outcome: RowMedianPolicyOutcome,
) -> list[PreprocessingRowAuditRow]:
    snapshot_base = _build_row_median_snapshot_base(
        input_profile=input_profile,
        output_missing_cell_count=outcome.output_missing_cell_count,
        imputed_cell_count=outcome.imputed_cell_count,
        stage_order=plan.stage_order,
        min_observed_values=outcome.min_observed_values,
    )
    records: list[PreprocessingRowAuditRow] = []
    for source_row_id, observed_value_count in outcome.dropped_row_observed_values:
        records.append(
            PreprocessingRowAuditRow(
                stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
                action="dropped",
                reason=(
                    "dropped because observed value count is below "
                    "missing_data.min_observed_values"
                ),
                source_row_id=source_row_id,
                site_id=source_row_id,
                retained=False,
                retained_row_id=pd.NA,
                source_rows=(source_row_id,),
                retained_row=pd.NA,
                parameter_snapshot={
                    **snapshot_base,
                    "observed_values": int(observed_value_count),
                },
            )
        )
    for row in outcome.imputed_rows:
        records.append(
            PreprocessingRowAuditRow(
                stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
                action="imputed",
                reason="missing values imputed with row median",
                source_row_id=row.row_id,
                site_id=row.row_id,
                retained=True,
                retained_row_id=row.row_id,
                source_rows=(row.row_id,),
                retained_row=row.row_id,
                parameter_snapshot={
                    **snapshot_base,
                    "imputed_columns": row.imputed_columns,
                    "imputed_cell_count": int(row.imputed_cell_count),
                    "row_median": float(outcome.row_medians_used[row.row_id]),
                },
            )
        )
    return records


def build_knn_audit_records(
    *,
    plan: PreprocessingPlan,
    input_profile: MissingDataInputProfile,
    outcome: KnnPolicyOutcome,
) -> list[PreprocessingRowAuditRow]:
    snapshot_base = _build_knn_snapshot_base(
        input_profile=input_profile,
        output_missing_cell_count=outcome.output_missing_cell_count,
        imputed_cell_count=outcome.imputed_cell_count,
        stage_order=plan.stage_order,
        k=outcome.k,
        distance=outcome.distance,
        max_missing_fraction_per_row=outcome.max_missing_fraction_per_row,
        no_overlap_policy=outcome.no_overlap_policy,
        no_overlap_policy_version=outcome.no_overlap_policy_version,
        nearest_neighbour_imputed_cell_count=(
            outcome.nearest_neighbour_imputed_cell_count
        ),
        column_mean_fallback_imputed_cell_count=(
            outcome.column_mean_fallback_imputed_cell_count
        ),
        nearest_neighbour_imputation_mask_hash=(
            outcome.nearest_neighbour_imputation_mask_hash
        ),
        column_mean_fallback_imputation_mask_hash=(
            outcome.column_mean_fallback_imputation_mask_hash
        ),
    )
    records: list[PreprocessingRowAuditRow] = []
    for source_row_id, missing_fraction_value in outcome.dropped_rows_missing_fraction:
        records.append(
            PreprocessingRowAuditRow(
                stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
                action="dropped",
                reason=(
                    "dropped because missing fraction exceeds "
                    "missing_data.max_missing_fraction_per_row"
                ),
                source_row_id=source_row_id,
                site_id=source_row_id,
                retained=False,
                retained_row_id=pd.NA,
                source_rows=(source_row_id,),
                retained_row=pd.NA,
                parameter_snapshot={
                    **snapshot_base,
                    "row_missing_fraction": float(missing_fraction_value),
                },
            )
        )
    for row in outcome.imputed_rows:
        records.append(
            PreprocessingRowAuditRow(
                stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
                action="imputed",
                reason="missing values imputed with knn",
                source_row_id=row.row_id,
                site_id=row.row_id,
                retained=True,
                retained_row_id=row.row_id,
                source_rows=(row.row_id,),
                retained_row=row.row_id,
                parameter_snapshot={
                    **snapshot_base,
                    "imputed_columns": row.imputed_columns,
                    "imputed_cell_count": int(row.imputed_cell_count),
                    "nearest_neighbour_imputed_columns": (
                        row.nearest_neighbour_imputed_columns
                    ),
                    "nearest_neighbour_imputed_cell_count": int(
                        len(row.nearest_neighbour_imputed_columns)
                    ),
                    "column_mean_fallback_columns": (row.column_mean_fallback_columns),
                    "column_mean_fallback_imputed_cell_count": int(
                        len(row.column_mean_fallback_columns)
                    ),
                    "fully_column_mean_fallback_imputed": bool(
                        row.imputed_cell_count > 0
                        and len(row.nearest_neighbour_imputed_columns) == 0
                        and len(row.column_mean_fallback_columns)
                        == int(row.imputed_cell_count)
                    ),
                },
            )
        )
    return records


def build_minprob_audit_records(
    *,
    plan: PreprocessingPlan,
    input_profile: MissingDataInputProfile,
    outcome: MinProbPolicyOutcome,
) -> list[PreprocessingRowAuditRow]:
    snapshot_base = _build_minprob_snapshot_base(
        input_profile=input_profile,
        output_missing_cell_count=outcome.output_missing_cell_count,
        imputed_cell_count=outcome.imputed_cell_count,
        stage_order=plan.stage_order,
        q=outcome.q,
        width=outcome.width,
        seed=outcome.seed,
        max_missing_fraction_per_row=outcome.max_missing_fraction_per_row,
    )
    records: list[PreprocessingRowAuditRow] = []
    for source_row_id, missing_fraction_value in outcome.dropped_rows_missing_fraction:
        records.append(
            PreprocessingRowAuditRow(
                stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
                action="dropped",
                reason=(
                    "dropped because missing fraction exceeds "
                    "missing_data.max_missing_fraction_per_row"
                ),
                source_row_id=source_row_id,
                site_id=source_row_id,
                retained=False,
                retained_row_id=pd.NA,
                source_rows=(source_row_id,),
                retained_row=pd.NA,
                parameter_snapshot={
                    **snapshot_base,
                    "row_missing_fraction": float(missing_fraction_value),
                },
            )
        )
    for row in outcome.imputed_rows:
        row_column_distributions = {
            column_name: outcome.per_column_distribution_parameters[column_name]
            for column_name in row.imputed_columns
            if column_name in outcome.per_column_distribution_parameters
        }
        records.append(
            PreprocessingRowAuditRow(
                stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
                action="imputed",
                reason="missing values imputed with minprob",
                source_row_id=row.row_id,
                site_id=row.row_id,
                retained=True,
                retained_row_id=row.row_id,
                source_rows=(row.row_id,),
                retained_row=row.row_id,
                parameter_snapshot={
                    **snapshot_base,
                    "imputed_columns": row.imputed_columns,
                    "imputed_cell_count": int(row.imputed_cell_count),
                    "column_distribution_parameters": row_column_distributions,
                },
            )
        )
    return records


def _build_row_audit_snapshot_common(
    *,
    input_profile: MissingDataInputProfile,
    output_missing_cell_count: int,
    imputed_cell_count: int,
    stage_order: tuple[str, ...],
) -> dict[str, JsonValue]:
    return {
        "input_missing_cell_count": int(input_profile.input_missing_cell_count),
        "output_missing_cell_count": int(output_missing_cell_count),
        "imputed_cell_count": int(imputed_cell_count),
        "affected_row_count": int(input_profile.affected_row_count),
        "affected_column_count": int(input_profile.affected_column_count),
        "missingness_mask_hash": str(input_profile.missingness_mask_hash),
        "stage_order": [str(stage) for stage in stage_order],
    }


def _build_row_median_snapshot_base(
    *,
    input_profile: MissingDataInputProfile,
    output_missing_cell_count: int,
    imputed_cell_count: int,
    stage_order: tuple[str, ...],
    min_observed_values: int,
) -> dict[str, JsonValue]:
    return {
        "missing_data_policy": MissingDataPolicy.IMPUTE_ROW_MEDIAN.value,
        "missing_data_min_observed_values": int(min_observed_values),
        **_build_row_audit_snapshot_common(
            input_profile=input_profile,
            output_missing_cell_count=output_missing_cell_count,
            imputed_cell_count=imputed_cell_count,
            stage_order=stage_order,
        ),
    }


def _build_knn_snapshot_base(
    *,
    input_profile: MissingDataInputProfile,
    output_missing_cell_count: int,
    imputed_cell_count: int,
    stage_order: tuple[str, ...],
    k: int,
    distance: str,
    max_missing_fraction_per_row: float,
    no_overlap_policy: str,
    no_overlap_policy_version: int,
    nearest_neighbour_imputed_cell_count: int,
    column_mean_fallback_imputed_cell_count: int,
    nearest_neighbour_imputation_mask_hash: str,
    column_mean_fallback_imputation_mask_hash: str,
) -> dict[str, JsonValue]:
    return {
        "missing_data_policy": MissingDataPolicy.IMPUTE_KNN.value,
        **_build_row_audit_snapshot_common(
            input_profile=input_profile,
            output_missing_cell_count=output_missing_cell_count,
            imputed_cell_count=imputed_cell_count,
            stage_order=stage_order,
        ),
        "k": int(k),
        "distance": str(distance),
        "max_missing_fraction_per_row": float(max_missing_fraction_per_row),
        "no_overlap_policy": str(no_overlap_policy),
        "no_overlap_policy_version": int(no_overlap_policy_version),
        "nearest_neighbour_imputed_cell_count": int(
            nearest_neighbour_imputed_cell_count
        ),
        "column_mean_fallback_imputed_cell_count": int(
            column_mean_fallback_imputed_cell_count
        ),
        "nearest_neighbour_imputation_mask_hash": str(
            nearest_neighbour_imputation_mask_hash
        ),
        "column_mean_fallback_imputation_mask_hash": str(
            column_mean_fallback_imputation_mask_hash
        ),
    }


def _build_minprob_snapshot_base(
    *,
    input_profile: MissingDataInputProfile,
    output_missing_cell_count: int,
    imputed_cell_count: int,
    stage_order: tuple[str, ...],
    q: float,
    width: float,
    seed: int,
    max_missing_fraction_per_row: float,
) -> dict[str, JsonValue]:
    return {
        "missing_data_policy": MissingDataPolicy.IMPUTE_MINPROB.value,
        **_build_row_audit_snapshot_common(
            input_profile=input_profile,
            output_missing_cell_count=output_missing_cell_count,
            imputed_cell_count=imputed_cell_count,
            stage_order=stage_order,
        ),
        "q": float(q),
        "width": float(width),
        "seed": int(seed),
        "max_missing_fraction_per_row": float(max_missing_fraction_per_row),
    }
