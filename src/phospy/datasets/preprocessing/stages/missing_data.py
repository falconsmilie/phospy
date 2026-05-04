"""Missing-data preprocessing stage for dataset building."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer

from phospy.api.configs import (
    DATASET_MISSING_DATA_POLICY_FORBID,
    DATASET_MISSING_DATA_POLICY_IMPUTE_KNN,
    DATASET_MISSING_DATA_POLICY_IMPUTE_MINPROB,
    DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
)
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    PreprocessingStageResult,
    PreprocessingState,
    append_row_audit_records,
)
from phospy.datasets.preprocessing.report_rows import report_rows_from_row_audit_rows
from phospy.datasets.preprocessing.report_schema import PreprocessingRowAuditRow
from phospy.datasets.processing_state import (
    MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1,
    JsonValue,
)
from phospy.errors.input import PhosPyInputError
from phospy.provenance.hashing import hash_table


class MissingDataStage:
    """Apply the configured missing-data policy to phospho/site tables."""

    stage_key = DATASET_PREPROCESSING_STAGE_MISSING_DATA

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        input_missing_mask = state.phospho.isna()
        input_missing_cell_count = int(input_missing_mask.to_numpy().sum())
        input_affected_row_ids = tuple(
            str(row_id)
            for row_id in state.phospho.index[input_missing_mask.any(axis=1)]
        )
        input_affected_column_ids = tuple(
            str(column_id)
            for column_id in state.phospho.columns[input_missing_mask.any(axis=0)]
        )
        affected_row_count = int(len(input_affected_row_ids))
        affected_column_count = int(len(input_affected_column_ids))
        missingness_mask_hash = _hash_missingness_mask(input_missing_mask)
        if state.plan.missing_data_policy == DATASET_MISSING_DATA_POLICY_FORBID:
            _fail_if_forbid_policy_has_missing_values(state.phospho)
            diagnostics = _build_missing_data_diagnostics(
                missing_data_policy=state.plan.missing_data_policy,
                imputation_method_id=None,
                imputation_method_family=None,
                input_missing_cell_count=input_missing_cell_count,
                output_missing_cell_count=input_missing_cell_count,
                imputed_cell_count=0,
                affected_row_ids=input_affected_row_ids,
                affected_column_ids=input_affected_column_ids,
                imputed_row_ids=(),
                imputed_column_ids=(),
                dropped_row_ids=(),
                random_seed=None,
                method_parameters={},
                matrix_scale_requirement=None,
                stage_order=state.plan.stage_order,
                missingness_mask_hash=missingness_mask_hash,
                left_censored_assumption=None,
                rows_not_imputable=(),
                row_medians_used={},
                per_column_distribution_parameters=None,
                dropped_rows_above_max_missing_fraction=(),
                neighbour_count=None,
                distance_metric=None,
            )
            return PreprocessingStageResult(
                state=state,
                diagnostics={
                    "dropped_row_ids": (),
                    "dropped_row_count": 0,
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "stage executed",
                    "diagnostics": diagnostics,
                },
            )

        if state.plan.missing_data_policy == DATASET_MISSING_DATA_POLICY_IMPUTE_MINPROB:
            return _run_minprob_missing_data_stage(
                state=state,
                input_missing_cell_count=input_missing_cell_count,
                input_affected_row_ids=input_affected_row_ids,
                input_affected_column_ids=input_affected_column_ids,
                affected_row_count=affected_row_count,
                affected_column_count=affected_column_count,
                missingness_mask_hash=missingness_mask_hash,
            )

        if state.plan.missing_data_policy == DATASET_MISSING_DATA_POLICY_IMPUTE_KNN:
            return _run_knn_missing_data_stage(
                state=state,
                input_missing_cell_count=input_missing_cell_count,
                input_affected_row_ids=input_affected_row_ids,
                input_affected_column_ids=input_affected_column_ids,
                affected_row_count=affected_row_count,
                affected_column_count=affected_column_count,
                missingness_mask_hash=missingness_mask_hash,
            )

        if (
            state.plan.missing_data_policy
            != DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config contains an unsupported "
                "missing_data.policy"
            )

        min_observed_values = state.plan.missing_data_min_observed_values
        if not isinstance(min_observed_values, int):
            raise PhosPyInputError(
                "dataset build request "
                "preprocessing_config.missing_data.min_observed_values must be an "
                "int when missing_data.policy='impute_row_median'"
            )
        if min_observed_values > state.phospho.shape[1]:
            raise PhosPyInputError(
                "dataset build request "
                "preprocessing_config.missing_data.min_observed_values "
                f"({min_observed_values}) cannot exceed the number of phospho "
                f"samples ({state.phospho.shape[1]})"
            )

        observed_counts = state.phospho.notna().sum(axis=1)
        retained_mask = observed_counts >= min_observed_values
        dropped_observed_counts = observed_counts.loc[~retained_mask]
        filtered_phospho = state.phospho.loc[retained_mask]
        filtered_site_metadata = state.site_metadata.loc[filtered_phospho.index]

        row_medians = filtered_phospho.median(axis=1, skipna=True)
        imputed = filtered_phospho.T.fillna(row_medians).T
        imputed_mask = filtered_phospho.isna() & imputed.notna()
        row_audit_records: list[PreprocessingRowAuditRow] = []
        dropped_row_details = tuple(
            (str(row_id), int(observed_value_count))
            for row_id, observed_value_count in dropped_observed_counts.items()
        )

        imputed_row_flags = imputed_mask.any(axis=1)
        imputed_rows = filtered_phospho.index[imputed_row_flags]
        row_medians_used = {
            str(row_id): float(row_medians.loc[row_id]) for row_id in imputed_rows
        }
        imputed_cell_count = int(imputed_mask.to_numpy().sum())
        imputed_row_ids = tuple(str(row_id) for row_id in imputed_rows.tolist())
        imputed_column_ids = tuple(
            str(column_name)
            for column_name in filtered_phospho.columns[imputed_mask.any(axis=0)]
        )
        rows_not_imputable = tuple(
            str(row_id)
            for row_id in imputed.index[
                imputed.isna().any(axis=1) & filtered_phospho.isna().any(axis=1)
            ].tolist()
        )
        dropped_row_ids = tuple(
            str(row_id) for row_id in dropped_observed_counts.index.tolist()
        )
        output_missing_cell_count = int(imputed.isna().to_numpy().sum())
        row_audit_snapshot_base = _build_row_audit_snapshot_base(
            missing_data_policy=DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
            missing_data_min_observed_values=int(min_observed_values),
            input_missing_cell_count=input_missing_cell_count,
            output_missing_cell_count=output_missing_cell_count,
            imputed_cell_count=imputed_cell_count,
            affected_row_count=affected_row_count,
            affected_column_count=affected_column_count,
            missingness_mask_hash=missingness_mask_hash,
            stage_order=state.plan.stage_order,
        )
        for source_row_id, observed_value_count in dropped_row_details:
            row_audit_records.append(
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
                        **row_audit_snapshot_base,
                        "observed_values": int(observed_value_count),
                    },
                )
            )
        for row_id in imputed_rows:
            source_row_id = str(row_id)
            row_imputed_mask = imputed_mask.loc[row_id]
            imputed_columns = tuple(
                str(column_name)
                for column_name in filtered_phospho.columns[row_imputed_mask].tolist()
            )
            row_audit_records.append(
                PreprocessingRowAuditRow(
                    stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
                    action="imputed",
                    reason="missing values imputed with row median",
                    source_row_id=source_row_id,
                    site_id=source_row_id,
                    retained=True,
                    retained_row_id=source_row_id,
                    source_rows=(source_row_id,),
                    retained_row=source_row_id,
                    parameter_snapshot={
                        **row_audit_snapshot_base,
                        "imputed_columns": imputed_columns,
                        "imputed_cell_count": int(row_imputed_mask.sum()),
                        "row_median": float(row_medians.loc[row_id]),
                    },
                )
            )
        next_state = append_row_audit_records(state, row_audit_records)
        diagnostics = _build_missing_data_diagnostics(
            missing_data_policy=state.plan.missing_data_policy,
            imputation_method_id="row_median",
            imputation_method_family="deterministic_row_statistic",
            input_missing_cell_count=input_missing_cell_count,
            output_missing_cell_count=output_missing_cell_count,
            imputed_cell_count=imputed_cell_count,
            affected_row_ids=input_affected_row_ids,
            affected_column_ids=input_affected_column_ids,
            imputed_row_ids=imputed_row_ids,
            imputed_column_ids=imputed_column_ids,
            dropped_row_ids=dropped_row_ids,
            random_seed=None,
            method_parameters={"min_observed_values": int(min_observed_values)},
            matrix_scale_requirement=None,
            stage_order=state.plan.stage_order,
            missingness_mask_hash=missingness_mask_hash,
            left_censored_assumption=False,
            rows_not_imputable=rows_not_imputable,
            row_medians_used=row_medians_used,
            per_column_distribution_parameters=None,
            dropped_rows_above_max_missing_fraction=(),
            neighbour_count=None,
            distance_metric=None,
        )
        return PreprocessingStageResult(
            state=replace(
                next_state,
                phospho=imputed,
                site_metadata=filtered_site_metadata,
            ),
            report_rows=report_rows_from_row_audit_rows(row_audit_records),
            diagnostics={
                "dropped_row_ids": dropped_row_ids,
                "dropped_row_count": int(len(dropped_row_ids)),
                "imputed_cell_count": imputed_cell_count,
                "imputed_row_ids": imputed_row_ids,
                "notes": "stage executed",
                "diagnostics": diagnostics,
            },
        )


def _run_knn_missing_data_stage(
    *,
    state: PreprocessingState,
    input_missing_cell_count: int,
    input_affected_row_ids: tuple[str, ...],
    input_affected_column_ids: tuple[str, ...],
    affected_row_count: int,
    affected_column_count: int,
    missingness_mask_hash: str,
) -> PreprocessingStageResult:
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
                f"affected column labels (preview): {_label_preview(all_missing_columns.tolist())}. "
                "adjust missing_data.max_missing_fraction_per_row or input data."
            )
        imputer = KNNImputer(
            n_neighbors=k_value,
            metric="nan_euclidean",
        )
        imputed_values = imputer.fit_transform(filtered_phospho)
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
            f"remaining rows (preview): {_label_preview(list(unresolved_row_ids))}. "
            f"remaining columns (preview): {_label_preview(list(unresolved_column_ids))}. "
            "adjust missing_data.max_missing_fraction_per_row, k, or input data."
        )

    dropped_row_ids = tuple(
        str(row_id) for row_id in dropped_missing_fraction.index.tolist()
    )
    rows_not_imputable = dropped_row_ids
    row_audit_records: list[PreprocessingRowAuditRow] = []
    row_audit_snapshot_base = _build_knn_row_audit_snapshot_base(
        input_missing_cell_count=input_missing_cell_count,
        output_missing_cell_count=output_missing_cell_count,
        imputed_cell_count=imputed_cell_count,
        affected_row_count=affected_row_count,
        affected_column_count=affected_column_count,
        missingness_mask_hash=missingness_mask_hash,
        stage_order=state.plan.stage_order,
        k=k_value,
        distance=distance_value,
        max_missing_fraction_per_row=max_missing_fraction_value,
    )
    for row_id, missing_fraction_value in dropped_missing_fraction.items():
        source_row_id = str(row_id)
        row_audit_records.append(
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
                    **row_audit_snapshot_base,
                    "row_missing_fraction": float(missing_fraction_value),
                },
            )
        )
    if not imputed.empty:
        for row_id in imputed.index[
            imputed_mask.any(axis=1).to_numpy(dtype=bool, copy=False)
        ]:
            source_row_id = str(row_id)
            row_imputed_mask = imputed_mask.loc[row_id]
            imputed_columns = tuple(
                str(column_name)
                for column_name in imputed.columns[row_imputed_mask].tolist()
            )
            row_audit_records.append(
                PreprocessingRowAuditRow(
                    stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
                    action="imputed",
                    reason="missing values imputed with knn",
                    source_row_id=source_row_id,
                    site_id=source_row_id,
                    retained=True,
                    retained_row_id=source_row_id,
                    source_rows=(source_row_id,),
                    retained_row=source_row_id,
                    parameter_snapshot={
                        **row_audit_snapshot_base,
                        "imputed_columns": imputed_columns,
                        "imputed_cell_count": int(row_imputed_mask.sum()),
                    },
                )
            )

    next_state = append_row_audit_records(state, row_audit_records)
    diagnostics = _build_missing_data_diagnostics(
        missing_data_policy=state.plan.missing_data_policy,
        imputation_method_id="knn",
        imputation_method_family="nearest_neighbour",
        input_missing_cell_count=input_missing_cell_count,
        output_missing_cell_count=output_missing_cell_count,
        imputed_cell_count=imputed_cell_count,
        affected_row_ids=input_affected_row_ids,
        affected_column_ids=input_affected_column_ids,
        imputed_row_ids=imputed_row_ids,
        imputed_column_ids=imputed_column_ids,
        dropped_row_ids=dropped_row_ids,
        random_seed=None,
        method_parameters={
            "k": int(k_value),
            "distance": distance_value,
            "max_missing_fraction_per_row": float(max_missing_fraction_value),
        },
        matrix_scale_requirement=None,
        stage_order=state.plan.stage_order,
        missingness_mask_hash=missingness_mask_hash,
        left_censored_assumption=False,
        rows_not_imputable=rows_not_imputable,
        row_medians_used={},
        per_column_distribution_parameters=None,
        dropped_rows_above_max_missing_fraction=dropped_row_ids,
        neighbour_count=int(k_value),
        distance_metric=distance_value,
    )
    return PreprocessingStageResult(
        state=replace(
            next_state,
            phospho=imputed,
            site_metadata=filtered_site_metadata,
        ),
        report_rows=report_rows_from_row_audit_rows(row_audit_records),
        diagnostics={
            "dropped_row_ids": dropped_row_ids,
            "dropped_row_count": int(len(dropped_row_ids)),
            "imputed_cell_count": imputed_cell_count,
            "imputed_row_ids": imputed_row_ids,
            "notes": "stage executed",
            "diagnostics": diagnostics,
        },
    )


def _run_minprob_missing_data_stage(
    *,
    state: PreprocessingState,
    input_missing_cell_count: int,
    input_affected_row_ids: tuple[str, ...],
    input_affected_column_ids: tuple[str, ...],
    affected_row_count: int,
    affected_column_count: int,
    missingness_mask_hash: str,
) -> PreprocessingStageResult:
    if state.plan.intensity_transform_policy != "log2":
        raise PhosPyInputError(
            "dataset build request preprocessing_config.missing_data.policy="
            "'impute_minprob' requires log2-scale values. Configure "
            "preprocessing_config.intensity_transform.policy='log2'."
        )

    q = state.plan.missing_data_q
    width = state.plan.missing_data_width
    seed = state.plan.missing_data_seed
    max_missing_fraction_per_row = state.plan.missing_data_max_missing_fraction_per_row
    if (
        q is None
        or width is None
        or seed is None
        or max_missing_fraction_per_row is None
    ):
        raise PhosPyInputError(
            "dataset build request preprocessing_config.missing_data.policy="
            "'impute_minprob' requires q, width, seed, and "
            "max_missing_fraction_per_row"
        )
    q_value = float(q)
    width_value = float(width)
    seed_value = int(seed)
    max_missing_fraction_value = float(max_missing_fraction_per_row)

    missing_fraction = state.phospho.isna().mean(axis=1)
    retained_mask = missing_fraction <= max_missing_fraction_value
    dropped_missing_fraction = missing_fraction.loc[~retained_mask]
    filtered_phospho = state.phospho.loc[retained_mask].copy(deep=True)
    filtered_site_metadata = state.site_metadata.loc[filtered_phospho.index]

    eps = float(np.finfo(float).eps)
    per_column_distribution_parameters: dict[str, dict[str, JsonValue]] = {}
    for column_index, column_name in enumerate(filtered_phospho.columns):
        column_label = str(column_name)
        column = filtered_phospho.loc[:, column_name]
        missing_mask = column.isna().to_numpy(dtype=bool, copy=False)
        missing_count = int(missing_mask.sum())
        observed_values = column.dropna().to_numpy(dtype=float, copy=False)
        observed_count = int(observed_values.size)
        if observed_count == 0 and missing_count > 0:
            raise PhosPyInputError(
                "dataset preprocessing stage 'missing_data' cannot apply "
                "missing_data.policy='impute_minprob' because column "
                f"{column_label!r} has no observed values after row filtering; "
                "adjust missing_data.max_missing_fraction_per_row or input data."
            )

        quantile_value = (
            float(np.quantile(observed_values, q_value)) if observed_count > 0 else 0.0
        )
        if observed_count > 1:
            observed_sd = float(np.std(observed_values, ddof=1))
        elif observed_count == 1:
            observed_sd = 0.0
        else:
            observed_sd = 0.0
        if not np.isfinite(observed_sd) or observed_sd <= 0.0:
            observed_sd = (
                float(np.std(observed_values, ddof=0)) if observed_count > 0 else 0.0
            )
        if not np.isfinite(observed_sd) or observed_sd <= 0.0:
            observed_sd = eps

        imputation_sd = max(observed_sd * width_value, eps)
        imputation_mean = float(quantile_value - (1.8 * imputation_sd))
        lower_tail = observed_values[observed_values <= quantile_value]
        lower_tail_mean = (
            float(np.mean(lower_tail))
            if int(lower_tail.size) > 0
            else float(quantile_value)
        )

        per_column_distribution_parameters[column_label] = {
            "observed_count": int(observed_count),
            "missing_count": int(missing_count),
            "q": float(q_value),
            "width": float(width_value),
            "lower_q_quantile": float(quantile_value),
            "lower_tail_mean": float(lower_tail_mean),
            "observed_sd": float(observed_sd),
            "imputation_mean": float(imputation_mean),
            "imputation_sd": float(imputation_sd),
        }

        if missing_count == 0:
            continue
        column_rng = np.random.default_rng(seed_value + column_index)
        draws = column_rng.normal(
            loc=imputation_mean,
            scale=imputation_sd,
            size=missing_count,
        )
        missing_index = filtered_phospho.index[missing_mask]
        filtered_phospho.loc[missing_index, column_name] = draws

    if filtered_phospho.empty:
        imputed_mask = filtered_phospho.isna() & filtered_phospho.notna()
    else:
        imputed_mask = (
            state.phospho.loc[retained_mask].isna() & filtered_phospho.notna()
        )
    imputed_cell_count = int(imputed_mask.to_numpy().sum())
    imputed_row_ids = (
        tuple(
            str(row_id)
            for row_id in filtered_phospho.index[
                imputed_mask.any(axis=1).to_numpy(dtype=bool, copy=False)
            ].tolist()
        )
        if not filtered_phospho.empty
        else ()
    )
    imputed_column_ids = (
        tuple(
            str(column_name)
            for column_name in filtered_phospho.columns[
                imputed_mask.any(axis=0).to_numpy(dtype=bool, copy=False)
            ].tolist()
        )
        if not filtered_phospho.empty
        else ()
    )
    rows_not_imputable = (
        tuple(
            str(row_id)
            for row_id in filtered_phospho.index[
                filtered_phospho.isna().any(axis=1).to_numpy(dtype=bool, copy=False)
            ].tolist()
        )
        if not filtered_phospho.empty
        else ()
    )
    dropped_row_ids = tuple(
        str(row_id) for row_id in dropped_missing_fraction.index.tolist()
    )
    output_missing_cell_count = int(filtered_phospho.isna().to_numpy().sum())
    row_audit_records: list[PreprocessingRowAuditRow] = []
    row_audit_snapshot_base = _build_minprob_row_audit_snapshot_base(
        input_missing_cell_count=input_missing_cell_count,
        output_missing_cell_count=output_missing_cell_count,
        imputed_cell_count=imputed_cell_count,
        affected_row_count=affected_row_count,
        affected_column_count=affected_column_count,
        missingness_mask_hash=missingness_mask_hash,
        stage_order=state.plan.stage_order,
        q=q_value,
        width=width_value,
        seed=seed_value,
        max_missing_fraction_per_row=max_missing_fraction_value,
    )
    for row_id, missing_fraction_value in dropped_missing_fraction.items():
        source_row_id = str(row_id)
        row_audit_records.append(
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
                    **row_audit_snapshot_base,
                    "row_missing_fraction": float(missing_fraction_value),
                },
            )
        )

    if not filtered_phospho.empty:
        for row_id in filtered_phospho.index[
            imputed_mask.any(axis=1).to_numpy(dtype=bool, copy=False)
        ]:
            source_row_id = str(row_id)
            row_imputed_mask = imputed_mask.loc[row_id]
            imputed_columns = tuple(
                str(column_name)
                for column_name in filtered_phospho.columns[row_imputed_mask].tolist()
            )
            row_column_distributions = {
                column_name: per_column_distribution_parameters[column_name]
                for column_name in imputed_columns
                if column_name in per_column_distribution_parameters
            }
            row_audit_records.append(
                PreprocessingRowAuditRow(
                    stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
                    action="imputed",
                    reason="missing values imputed with minprob",
                    source_row_id=source_row_id,
                    site_id=source_row_id,
                    retained=True,
                    retained_row_id=source_row_id,
                    source_rows=(source_row_id,),
                    retained_row=source_row_id,
                    parameter_snapshot={
                        **row_audit_snapshot_base,
                        "imputed_columns": imputed_columns,
                        "imputed_cell_count": int(row_imputed_mask.sum()),
                        "column_distribution_parameters": row_column_distributions,
                    },
                )
            )

    next_state = append_row_audit_records(state, row_audit_records)
    diagnostics = _build_missing_data_diagnostics(
        missing_data_policy=state.plan.missing_data_policy,
        imputation_method_id="minprob",
        imputation_method_family="left_censored_random",
        input_missing_cell_count=input_missing_cell_count,
        output_missing_cell_count=output_missing_cell_count,
        imputed_cell_count=imputed_cell_count,
        affected_row_ids=input_affected_row_ids,
        affected_column_ids=input_affected_column_ids,
        imputed_row_ids=imputed_row_ids,
        imputed_column_ids=imputed_column_ids,
        dropped_row_ids=dropped_row_ids,
        random_seed=seed_value,
        method_parameters={
            "q": float(q_value),
            "width": float(width_value),
            "seed": int(seed_value),
            "max_missing_fraction_per_row": float(max_missing_fraction_value),
        },
        matrix_scale_requirement="log2",
        stage_order=state.plan.stage_order,
        missingness_mask_hash=missingness_mask_hash,
        left_censored_assumption=True,
        rows_not_imputable=rows_not_imputable,
        row_medians_used={},
        per_column_distribution_parameters=per_column_distribution_parameters,
        dropped_rows_above_max_missing_fraction=dropped_row_ids,
        neighbour_count=None,
        distance_metric=None,
    )
    return PreprocessingStageResult(
        state=replace(
            next_state,
            phospho=filtered_phospho,
            site_metadata=filtered_site_metadata,
        ),
        report_rows=report_rows_from_row_audit_rows(row_audit_records),
        diagnostics={
            "dropped_row_ids": dropped_row_ids,
            "dropped_row_count": int(len(dropped_row_ids)),
            "imputed_cell_count": imputed_cell_count,
            "imputed_row_ids": imputed_row_ids,
            "notes": "stage executed",
            "diagnostics": diagnostics,
        },
    )


def _fail_if_forbid_policy_has_missing_values(phospho: pd.DataFrame) -> None:
    missing_mask = phospho.isna()
    missing_cell_count = int(missing_mask.to_numpy().sum())
    if missing_cell_count == 0:
        return

    affected_row_count = int(missing_mask.any(axis=1).sum())
    affected_column_count = int(missing_mask.any(axis=0).sum())
    affected_row_preview = _label_preview(
        phospho.index[missing_mask.any(axis=1)].tolist()
    )
    affected_column_preview = _label_preview(
        phospho.columns[missing_mask.any(axis=0)].tolist()
    )
    raise PhosPyInputError(
        "dataset preprocessing stage 'missing_data' rejected missing values because "
        "missing_data.policy='forbid'; "
        f"found {missing_cell_count} missing values across "
        f"{affected_row_count} rows and {affected_column_count} columns. "
        f"affected row labels (preview): {affected_row_preview}. "
        f"affected column labels (preview): {affected_column_preview}. "
        "choose missing_data.policy='impute_row_median' or clean the input data."
    )


def _label_preview(values: list[object], *, max_items: int = 3) -> str:
    if not values:
        return "none"
    rendered = [repr(str(value)) for value in values[:max_items]]
    remaining_count = len(values) - len(rendered)
    if remaining_count > 0:
        rendered.append(f"+{remaining_count} more")
    return ", ".join(rendered)


def _hash_missingness_mask(mask: pd.DataFrame) -> str:
    """Return stable fingerprint for input missingness structure."""

    return hash_table(
        mask.astype("int8"),
        name="missing_data.input_missingness_mask",
    )


def _build_missing_data_diagnostics(
    *,
    missing_data_policy: str,
    imputation_method_id: str | None,
    imputation_method_family: str | None,
    input_missing_cell_count: int,
    output_missing_cell_count: int,
    imputed_cell_count: int,
    affected_row_ids: tuple[str, ...],
    affected_column_ids: tuple[str, ...],
    imputed_row_ids: tuple[str, ...],
    imputed_column_ids: tuple[str, ...],
    dropped_row_ids: tuple[str, ...],
    random_seed: int | None,
    method_parameters: dict[str, JsonValue],
    matrix_scale_requirement: str | None,
    stage_order: tuple[str, ...],
    missingness_mask_hash: str,
    left_censored_assumption: bool | None,
    rows_not_imputable: tuple[str, ...],
    row_medians_used: dict[str, float],
    per_column_distribution_parameters: dict[str, dict[str, JsonValue]] | None,
    dropped_rows_above_max_missing_fraction: tuple[str, ...],
    neighbour_count: int | None,
    distance_metric: str | None,
) -> dict[str, JsonValue]:
    diagnostics: dict[str, JsonValue] = {
        "diagnostics_schema_version": MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1,
        "missing_data_policy": missing_data_policy,
        "imputation_method_id": imputation_method_id,
        "imputation_method_family": imputation_method_family,
        "input_missing_cell_count": int(input_missing_cell_count),
        "output_missing_cell_count": int(output_missing_cell_count),
        "imputed_cell_count": int(imputed_cell_count),
        "affected_row_count": int(len(affected_row_ids)),
        "affected_column_count": int(len(affected_column_ids)),
        "affected_row_ids": list(affected_row_ids),
        "affected_column_ids": list(affected_column_ids),
        "imputed_row_ids": list(imputed_row_ids),
        "imputed_column_ids": list(imputed_column_ids),
        "dropped_row_ids": list(dropped_row_ids),
        "random_seed": random_seed,
        "method_parameters": dict(method_parameters),
        "matrix_scale_requirement": matrix_scale_requirement,
        "stage_order": list(stage_order),
        "missingness_mask_hash": missingness_mask_hash,
        "left_censored_assumption": left_censored_assumption,
        "rows_not_imputable": list(rows_not_imputable),
        "row_medians_used": {
            str(row_id): float(row_median)
            for row_id, row_median in row_medians_used.items()
        },
        "dropped_rows_above_max_missing_fraction": list(
            dropped_rows_above_max_missing_fraction
        ),
        "neighbour_count": neighbour_count,
        "distance_metric": distance_metric,
    }
    if per_column_distribution_parameters is not None:
        per_column_distribution_payload: dict[str, JsonValue] = {
            str(column_name): dict(parameters)
            for column_name, parameters in per_column_distribution_parameters.items()
        }
        diagnostics["per_column_distribution_parameters"] = (
            per_column_distribution_payload
        )
    return diagnostics


def _build_row_audit_snapshot_base(
    *,
    missing_data_policy: str,
    missing_data_min_observed_values: int,
    input_missing_cell_count: int,
    output_missing_cell_count: int,
    imputed_cell_count: int,
    affected_row_count: int,
    affected_column_count: int,
    missingness_mask_hash: str,
    stage_order: tuple[str, ...],
) -> dict[str, JsonValue]:
    return {
        "missing_data_policy": missing_data_policy,
        "missing_data_min_observed_values": int(missing_data_min_observed_values),
        "input_missing_cell_count": int(input_missing_cell_count),
        "output_missing_cell_count": int(output_missing_cell_count),
        "imputed_cell_count": int(imputed_cell_count),
        "affected_row_count": int(affected_row_count),
        "affected_column_count": int(affected_column_count),
        "missingness_mask_hash": str(missingness_mask_hash),
        "stage_order": [str(stage) for stage in stage_order],
    }


def _build_knn_row_audit_snapshot_base(
    *,
    input_missing_cell_count: int,
    output_missing_cell_count: int,
    imputed_cell_count: int,
    affected_row_count: int,
    affected_column_count: int,
    missingness_mask_hash: str,
    stage_order: tuple[str, ...],
    k: int,
    distance: str,
    max_missing_fraction_per_row: float,
) -> dict[str, JsonValue]:
    return {
        "missing_data_policy": DATASET_MISSING_DATA_POLICY_IMPUTE_KNN,
        "input_missing_cell_count": int(input_missing_cell_count),
        "output_missing_cell_count": int(output_missing_cell_count),
        "imputed_cell_count": int(imputed_cell_count),
        "affected_row_count": int(affected_row_count),
        "affected_column_count": int(affected_column_count),
        "missingness_mask_hash": str(missingness_mask_hash),
        "stage_order": [str(stage) for stage in stage_order],
        "k": int(k),
        "distance": str(distance),
        "max_missing_fraction_per_row": float(max_missing_fraction_per_row),
    }


def _build_minprob_row_audit_snapshot_base(
    *,
    input_missing_cell_count: int,
    output_missing_cell_count: int,
    imputed_cell_count: int,
    affected_row_count: int,
    affected_column_count: int,
    missingness_mask_hash: str,
    stage_order: tuple[str, ...],
    q: float,
    width: float,
    seed: int,
    max_missing_fraction_per_row: float,
) -> dict[str, JsonValue]:
    return {
        "missing_data_policy": DATASET_MISSING_DATA_POLICY_IMPUTE_MINPROB,
        "input_missing_cell_count": int(input_missing_cell_count),
        "output_missing_cell_count": int(output_missing_cell_count),
        "imputed_cell_count": int(imputed_cell_count),
        "affected_row_count": int(affected_row_count),
        "affected_column_count": int(affected_column_count),
        "missingness_mask_hash": str(missingness_mask_hash),
        "stage_order": [str(stage) for stage in stage_order],
        "q": float(q),
        "width": float(width),
        "seed": int(seed),
        "max_missing_fraction_per_row": float(max_missing_fraction_per_row),
    }


__all__ = ["MissingDataStage"]
