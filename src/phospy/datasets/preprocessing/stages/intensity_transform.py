"""Intensity-transform stage for dataset preprocessing."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd

from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.errors.input import PhosPyInputError
from phospy.policy_models import IntensityTransformPolicy
from phospy.provenance.hashing import hash_table


class IntensityTransformStage:
    """Apply configured quantitative intensity transform to phospho values."""

    stage_key = DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        policy = state.plan.intensity_transform_policy
        if policy is IntensityTransformPolicy.IDENTITY:
            identity_diagnostics: dict[str, object] = {
                "policy": policy.value,
                "pseudocount": float(state.plan.intensity_transform_pseudocount),
                "affected_matrices": ["phospho"],
                "input_phospho_hash": hash_table(
                    state.phospho,
                    name="intensity_transform.input.phospho",
                ),
                "output_phospho_hash": hash_table(
                    state.phospho,
                    name="intensity_transform.output.phospho",
                ),
            }
            if state.total is not None:
                identity_diagnostics["affected_matrices"] = ["phospho", "total"]
                identity_diagnostics["input_total_hash"] = hash_table(
                    state.total,
                    name="intensity_transform.input.total",
                )
                identity_diagnostics["output_total_hash"] = hash_table(
                    state.total,
                    name="intensity_transform.output.total",
                )
            return PreprocessingStageResult(
                state=state,
                diagnostics={
                    "dropped_row_ids": (),
                    "dropped_row_count": 0,
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "stage executed",
                    "diagnostics": identity_diagnostics,
                },
            )
        if policy is not IntensityTransformPolicy.LOG2:
            raise PhosPyInputError(
                "dataset build request preprocessing_config contains an unsupported "
                "intensity_transform.policy"
            )

        pseudocount = float(state.plan.intensity_transform_pseudocount)
        if not np.isfinite(pseudocount):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.intensity_transform."
                "pseudocount must be finite"
            )
        if pseudocount < 0:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.intensity_transform."
                "pseudocount must be greater than or equal to 0"
            )

        _require_numeric_columns(
            state.phospho,
            field_name="phospho",
            operation_name="log2",
        )
        transformed_phospho = _apply_log2(
            state.phospho,
            pseudocount=pseudocount,
            field_name="phospho",
        )
        transformed_total = state.total
        affected_matrices = ["phospho"]
        if transformed_total is not None:
            _require_numeric_columns(
                transformed_total,
                field_name="total",
                operation_name="log2",
            )
            transformed_total = _apply_log2(
                transformed_total,
                pseudocount=pseudocount,
                field_name="total",
            )
            affected_matrices.append("total")
        next_state = replace(
            state,
            phospho=transformed_phospho,
            total=transformed_total,
        )
        diagnostics: dict[str, object] = {
            "policy": policy.value,
            "pseudocount": pseudocount,
            "affected_matrices": affected_matrices,
            "input_phospho_hash": hash_table(
                state.phospho,
                name="intensity_transform.input.phospho",
            ),
            "output_phospho_hash": hash_table(
                transformed_phospho,
                name="intensity_transform.output.phospho",
            ),
        }
        if state.total is not None and transformed_total is not None:
            diagnostics["input_total_hash"] = hash_table(
                state.total,
                name="intensity_transform.input.total",
            )
            diagnostics["output_total_hash"] = hash_table(
                transformed_total,
                name="intensity_transform.output.total",
            )
        return PreprocessingStageResult(
            state=next_state,
            diagnostics={
                "dropped_row_ids": (),
                "dropped_row_count": 0,
                "imputed_cell_count": 0,
                "imputed_row_ids": (),
                "notes": "stage executed",
                "diagnostics": diagnostics,
            },
        )


def _apply_log2(
    matrix: pd.DataFrame,
    *,
    pseudocount: float,
    field_name: str,
) -> pd.DataFrame:
    adjusted = matrix + pseudocount
    invalid_mask = matrix.notna() & (adjusted <= 0)
    if bool(invalid_mask.to_numpy().any()):
        invalid_positions = np.argwhere(invalid_mask.to_numpy())
        first_row, first_col = invalid_positions[0]
        row_label = str(matrix.index[first_row])
        column_label = str(matrix.columns[first_col])
        offending_value = matrix.iat[int(first_row), int(first_col)]
        raise PhosPyInputError(
            "dataset build request preprocessing intensity_transform.policy='log2' "
            "requires all non-missing values plus pseudocount to be greater than 0. "
            f"First invalid value at {field_name}[{row_label!r}, {column_label!r}]="
            f"{offending_value!r} with pseudocount={pseudocount}. Increase pseudocount "
            "or adjust/remove non-positive values before applying log2."
        )
    transformed = np.log2(adjusted)
    transformed.index = matrix.index.copy()
    transformed.columns = matrix.columns.copy()
    return transformed


def _require_numeric_columns(
    matrix: pd.DataFrame,
    *,
    field_name: str,
    operation_name: str,
) -> None:
    non_numeric = [
        str(column)
        for column in matrix.columns
        if (
            not pd.api.types.is_numeric_dtype(matrix[column])
            or pd.api.types.is_bool_dtype(matrix[column])
        )
    ]
    if non_numeric:
        raise PhosPyInputError(
            "dataset build request preprocessing "
            f"intensity_transform.policy='{operation_name}' requires numeric "
            f"{field_name} columns. Non-numeric columns: {', '.join(non_numeric)}"
        )


__all__ = ["IntensityTransformStage"]
