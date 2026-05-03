"""Missing-data preprocessing stage for dataset building."""

from __future__ import annotations

from dataclasses import replace

import pandas as pd

from phospy.api.configs import (
    DATASET_MISSING_DATA_POLICY_FORBID,
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
        for row_id, observed_value_count in dropped_observed_counts.items():
            source_row_id = str(row_id)
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
                        "missing_data_policy": DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
                        "missing_data_min_observed_values": int(min_observed_values),
                        "observed_values": int(observed_value_count),
                    },
                )
            )

        imputed_row_flags = imputed_mask.any(axis=1)
        imputed_rows = filtered_phospho.index[imputed_row_flags]
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
                        "missing_data_policy": DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
                        "missing_data_min_observed_values": int(min_observed_values),
                        "imputed_columns": imputed_columns,
                        "imputed_cell_count": int(row_imputed_mask.sum()),
                        "row_median": float(row_medians.loc[row_id]),
                    },
                )
            )
        next_state = append_row_audit_records(state, row_audit_records)
        dropped_row_ids = tuple(
            str(row_id) for row_id in dropped_observed_counts.index.tolist()
        )
        output_missing_cell_count = int(imputed.isna().to_numpy().sum())
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
) -> dict[str, JsonValue]:
    return {
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
    }


__all__ = ["MissingDataStage"]
