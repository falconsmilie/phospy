"""Intensity-transform stage for dataset preprocessing."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd

from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
    PreprocessingStateTableKey,
)
from phospy.datasets.preprocessing.stage_contract import PreprocessingStageContract
from phospy.errors.input import PhosPyInputError
from phospy.policy_models import IntensityTransformPolicy
from phospy.provenance.hashing import hash_table
from phospy.transformations.models import IntensityScaleKind


class IntensityTransformStage:
    """Apply configured quantitative intensity transform to phospho values."""

    stage_key = DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        policy = state.plan.intensity_transform_policy
        if policy is IntensityTransformPolicy.IDENTITY:
            identity_diagnostics: dict[str, object] = {
                "policy": policy.value,
                "pseudocount": float(state.plan.intensity_transform_pseudocount),
                "output_intensity_scale_kind": IntensityScaleKind.LINEAR.value,
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
            "output_intensity_scale_kind": IntensityScaleKind.LOG2.value,
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
    transformed_values = np.log2(adjusted.to_numpy(dtype=float, copy=False))
    return pd.DataFrame(
        transformed_values,
        index=matrix.index.copy(),
        columns=matrix.columns.copy(),
    )


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


def _resolve_operation(plan: PreprocessingPlan) -> str:
    return plan.intensity_transform_policy.value


def _resolve_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    return {"pseudocount": float(plan.intensity_transform_pseudocount)}


INTENSITY_TRANSFORM_STAGE_CONTRACT = PreprocessingStageContract(
    stage_key=DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    display_label=DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    provenance_stage=DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    operation_name=_resolve_operation,
    serialize_parameters=_resolve_parameters,
    consumed_input_tables=(
        PreprocessingStateTableKey.DATASET_PHOSPHO,
        PreprocessingStateTableKey.DATASET_TOTAL,
    ),
    produced_output_tables=(
        PreprocessingStateTableKey.DATASET_PHOSPHO,
        PreprocessingStateTableKey.DATASET_TOTAL,
    ),
    stage_factory=IntensityTransformStage,
    backend="numpy",
    diagnostics_metadata={
        "known_diagnostics_fields": (
            "policy",
            "pseudocount",
            "output_intensity_scale_kind",
            "affected_matrices",
            "input_phospho_hash",
            "output_phospho_hash",
            "input_total_hash",
            "output_total_hash",
        )
    },
)


__all__ = ["INTENSITY_TRANSFORM_STAGE_CONTRACT", "IntensityTransformStage"]
