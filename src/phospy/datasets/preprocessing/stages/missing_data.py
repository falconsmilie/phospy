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
from phospy.errors.input import PhosPyInputError


class MissingDataStage:
    """Apply the configured missing-data policy to phospho/site tables."""

    stage_key = DATASET_PREPROCESSING_STAGE_MISSING_DATA

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        if state.plan.missing_data_policy == DATASET_MISSING_DATA_POLICY_FORBID:
            return PreprocessingStageResult(
                state=state,
                diagnostics={
                    "dropped_row_ids": (),
                    "dropped_row_count": 0,
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "stage executed",
                    "diagnostics": {
                        "min_observed_values": state.plan.missing_data_min_observed_values,
                        "imputed_row_ids": [],
                    },
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
        row_audit_records: list[dict[str, object]] = []
        for row_id, observed_value_count in dropped_observed_counts.items():
            source_row_id = str(row_id)
            row_audit_records.append(
                {
                    "stage": DATASET_PREPROCESSING_STAGE_MISSING_DATA,
                    "action": "dropped",
                    "reason": (
                        "dropped because observed value count is below "
                        "missing_data.min_observed_values"
                    ),
                    "source_row_id": source_row_id,
                    "site_id": source_row_id,
                    "retained": False,
                    "retained_row_id": pd.NA,
                    "source_rows": (source_row_id,),
                    "retained_row": pd.NA,
                    "parameter_snapshot": {
                        "missing_data_policy": DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
                        "missing_data_min_observed_values": int(min_observed_values),
                        "observed_values": int(observed_value_count),
                    },
                }
            )

        imputed_row_flags = imputed_mask.any(axis=1)
        imputed_rows = filtered_phospho.index[imputed_row_flags]
        imputed_cell_count = int(imputed_mask.to_numpy().sum())
        imputed_row_ids = tuple(str(row_id) for row_id in imputed_rows.tolist())
        for row_id in imputed_rows:
            source_row_id = str(row_id)
            row_imputed_mask = imputed_mask.loc[row_id]
            imputed_columns = tuple(
                str(column_name)
                for column_name in filtered_phospho.columns[row_imputed_mask].tolist()
            )
            row_audit_records.append(
                {
                    "stage": DATASET_PREPROCESSING_STAGE_MISSING_DATA,
                    "action": "imputed",
                    "reason": "missing values imputed with row median",
                    "source_row_id": source_row_id,
                    "site_id": source_row_id,
                    "retained": True,
                    "retained_row_id": source_row_id,
                    "source_rows": (source_row_id,),
                    "retained_row": source_row_id,
                    "parameter_snapshot": {
                        "missing_data_policy": DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
                        "missing_data_min_observed_values": int(min_observed_values),
                        "imputed_columns": imputed_columns,
                        "imputed_cell_count": int(row_imputed_mask.sum()),
                        "row_median": float(row_medians.loc[row_id]),
                    },
                }
            )
        next_state = append_row_audit_records(state, row_audit_records)
        dropped_row_ids = tuple(
            str(row_id) for row_id in dropped_observed_counts.index.tolist()
        )
        return PreprocessingStageResult(
            state=replace(
                next_state,
                phospho=imputed,
                site_metadata=filtered_site_metadata,
            ),
            diagnostics={
                "dropped_row_ids": dropped_row_ids,
                "dropped_row_count": int(len(dropped_row_ids)),
                "imputed_cell_count": imputed_cell_count,
                "imputed_row_ids": imputed_row_ids,
                "notes": "stage executed",
                "diagnostics": {
                    "min_observed_values": state.plan.missing_data_min_observed_values,
                    "imputed_row_ids": [row_id for row_id in imputed_row_ids],
                },
            },
        )


__all__ = ["MissingDataStage"]
