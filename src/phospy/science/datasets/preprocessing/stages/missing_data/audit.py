"""Row-audit construction for missing-data stage."""

from __future__ import annotations

from typing import cast

import pandas as pd

from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    PreprocessingPlan,
)
from phospy.science.datasets.preprocessing.policy_models import MissingDataPolicy
from phospy.science.datasets.preprocessing.report_schema import PreprocessingRowAuditRow
from phospy.science.datasets.processing_state import JsonValue

from .models import (
    GroupAwarePolicyOutcome,
    GroupMissingnessClassification,
    KnnPolicyOutcome,
    MinProbPolicyOutcome,
    MissingDataInputProfile,
    RowMedianPolicyOutcome,
)


def build_group_aware_audit_records(
    *,
    plan: PreprocessingPlan,
    input_profile: MissingDataInputProfile,
    outcome: GroupAwarePolicyOutcome,
) -> list[PreprocessingRowAuditRow]:
    """Build row-level routing and mechanism audit records."""

    resolved_group_samples: dict[str, JsonValue] = {
        group: cast(JsonValue, list(samples))
        for group, samples in outcome.routing.resolved_groups.original_sample_order_by_group.items()
    }
    sample_group_by_column = outcome.routing.resolved_groups.group_by_original_sample
    routing_facts_by_row = {
        row_id: tuple(
            fact for fact in outcome.routing.group_facts if fact.row_id == row_id
        )
        for row_id in (
            *outcome.routing.retained_row_ids,
            *outcome.routing.dropped_row_ids,
        )
    }
    snapshot_base: dict[str, JsonValue] = {
        "missing_data_policy": MissingDataPolicy.IMPUTE_GROUP_AWARE.value,
        **_build_row_audit_snapshot_common(
            input_profile=input_profile,
            output_missing_cell_count=outcome.output_missing_cell_count,
            imputed_cell_count=outcome.imputed_cell_count,
            stage_order=plan.stage_order,
        ),
        "group_column": str(plan.missing_data_group_column),
        "min_partial_observed_fraction": float(
            outcome.routing.min_partial_observed_fraction
        ),
        "min_reference_observed_fraction": float(
            outcome.routing.min_reference_observed_fraction
        ),
        "k": int(outcome.k),
        "distance": str(outcome.distance),
        "no_overlap_policy": str(outcome.no_overlap_policy),
        "q": float(outcome.q),
        "width": float(outcome.width),
        "seed": int(outcome.seed),
        "observed_group_sizes": {
            group: len(samples)
            for group, samples in outcome.routing.resolved_groups.original_sample_order_by_group.items()
        },
        "resolved_group_samples": resolved_group_samples,
        "retained_row_count": len(outcome.routing.retained_row_ids),
        "dropped_unsupported_row_count": len(outcome.routing.dropped_row_ids),
        "knn_routed_cell_count": outcome.knn_target_cell_count,
        "minprob_routed_cell_count": outcome.minprob_target_cell_count,
        "knn_imputed_cell_count": outcome.knn_imputed_cell_count,
        "minprob_imputed_cell_count": outcome.minprob_imputed_cell_count,
        "knn_target_mask_hash": outcome.knn_target_mask_hash,
        "minprob_target_mask_hash": outcome.minprob_target_mask_hash,
        "knn_imputation_mask_hash": outcome.knn_imputation_mask_hash,
        "minprob_imputation_mask_hash": outcome.minprob_imputation_mask_hash,
        "imputation_mask_hash": outcome.imputation_mask_hash,
        "minprob_left_censored_assumption": True,
    }
    records: list[PreprocessingRowAuditRow] = []
    for dropped in outcome.routing.dropped_row_reasons:
        partial_groups = tuple(
            group
            for group, classification in dropped.reasons_by_group.items()
            if classification is GroupMissingnessClassification.UNSUPPORTED_PARTIAL
        )
        absence_groups = tuple(
            group
            for group, classification in dropped.reasons_by_group.items()
            if classification
            is GroupMissingnessClassification.UNSUPPORTED_FULLY_MISSING
        )
        if partial_groups and absence_groups:
            reason = "unsupported partial coverage and unsupported fully-missing group"
        elif partial_groups:
            reason = "unsupported partial coverage"
        else:
            reason = "unsupported fully-missing group"
        records.append(
            PreprocessingRowAuditRow(
                stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
                action="dropped",
                reason=reason,
                source_row_id=dropped.row_id,
                site_id=dropped.row_id,
                retained=False,
                retained_row_id=pd.NA,
                source_rows=(dropped.row_id,),
                retained_row=pd.NA,
                parameter_snapshot={
                    **snapshot_base,
                    "unsupported_groups": {
                        str(group): category.value
                        for group, category in dropped.route_categories_by_group.items()
                    },
                    "unsupported_partial_groups": partial_groups,
                    "unsupported_absence_groups": absence_groups,
                    "observed_finite_count_by_group": {
                        fact.group_label: fact.observed_finite_count
                        for fact in routing_facts_by_row[dropped.row_id]
                    },
                    "observed_fraction_by_group": {
                        fact.group_label: fact.observed_fraction
                        for fact in routing_facts_by_row[dropped.row_id]
                    },
                },
            )
        )
    for row in outcome.imputed_rows:
        minprob_row_mask = cast(
            pd.Series,
            outcome.minprob_imputed_mask.loc[row.row_id],
        )
        minprob_columns = tuple(
            str(column)
            for column in outcome.phospho.columns[
                minprob_row_mask.to_numpy(dtype=bool, copy=False)
            ].tolist()
        )
        records.append(
            PreprocessingRowAuditRow(
                stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
                action="imputed",
                reason="group-aware row retained/imputed",
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
                    "knn_imputed_columns": row.nearest_neighbour_imputed_columns,
                    "minprob_imputed_columns": minprob_columns,
                    "knn_affected_groups": list(
                        dict.fromkeys(
                            sample_group_by_column[column]
                            for column in row.nearest_neighbour_imputed_columns
                        )
                    ),
                    "minprob_affected_groups": list(
                        dict.fromkeys(
                            sample_group_by_column[column] for column in minprob_columns
                        )
                    ),
                    "mechanisms": [
                        mechanism
                        for mechanism, columns in (
                            (
                                "partial_observation_knn",
                                row.nearest_neighbour_imputed_columns,
                            ),
                            ("asymmetric_absence_minprob", minprob_columns),
                        )
                        if columns
                    ],
                    "knn_imputed_cell_count_for_row": len(
                        row.nearest_neighbour_imputed_columns
                    ),
                    "minprob_imputed_cell_count_for_row": len(minprob_columns),
                    "minprob_column_distribution_parameters": {
                        column: outcome.per_column_distribution_parameters[column]
                        for column in minprob_columns
                        if column in outcome.per_column_distribution_parameters
                    },
                },
            )
        )
    return records


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
