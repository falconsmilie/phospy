"""Log2 transformer for dataset quantitative matrices."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.transformations.contracts import TransformationResult
from phospy.science.transformations.models import (
    IntensityScaleKind,
    IntensityScaleState,
    MatrixIntensityScaleState,
)

_LOG2_ESTABLISHED_BY = "phospy.science.transformations.transformers.log2"


class Log2Transformer:
    """Apply validated log2 transformation with an explicit pseudocount."""

    preserves_input_scale_state = False
    changes_numeric_values = True
    requires_established_input_state = False

    def __init__(self, *, pseudocount: float) -> None:
        self._pseudocount = _validate_pseudocount(pseudocount)

    @property
    def pseudocount(self) -> float:
        return self._pseudocount

    def run(
        self,
        phospho: pd.DataFrame,
        total: pd.DataFrame | None = None,
    ) -> TransformationResult:
        _require_numeric_columns(
            phospho,
            field_name="phospho",
            operation_name="log2",
        )
        transformed_phospho = _apply_log2(
            phospho,
            pseudocount=self._pseudocount,
            field_name="phospho",
        )
        transformed_total = total
        affected_matrices = ["phospho"]
        if transformed_total is not None:
            _require_numeric_columns(
                transformed_total,
                field_name="total",
                operation_name="log2",
            )
            transformed_total = _apply_log2(
                transformed_total,
                pseudocount=self._pseudocount,
                field_name="total",
            )
            affected_matrices.append("total")

        state = IntensityScaleState(
            phospho=MatrixIntensityScaleState.log2(established_by=_LOG2_ESTABLISHED_BY),
            total=(
                MatrixIntensityScaleState.log2(established_by=_LOG2_ESTABLISHED_BY)
                if transformed_total is not None
                else None
            ),
        )
        transformer_name = f"{self.__class__.__module__}.{self.__class__.__qualname__}"
        return TransformationResult(
            phospho=transformed_phospho,
            total=transformed_total,
            state=state,
            provenance={
                "policy": "log2",
                "operation": "log2",
                "pseudocount": self._pseudocount,
                "output_intensity_scale_kind": IntensityScaleKind.LOG2.value,
                "affected_matrices": affected_matrices,
                "transformer_name": transformer_name,
                "transformer_state": _serialize_intensity_scale_state(state),
            },
        )


def _validate_pseudocount(pseudocount: float) -> float:
    resolved = float(pseudocount)
    if not np.isfinite(resolved):
        raise PhosPyInputError(
            "dataset build request preprocessing_config.intensity_transform."
            "pseudocount must be finite"
        )
    if resolved < 0:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.intensity_transform."
            "pseudocount must be greater than or equal to 0"
        )
    return resolved


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


def _serialize_intensity_scale_state(state: IntensityScaleState) -> dict[str, object]:
    return {
        "phospho": {
            "kind": state.phospho.kind.value,
            "transformed": bool(state.phospho.transformed),
            "established_by": state.phospho.established_by,
        },
        "total": (
            None
            if state.total is None
            else {
                "kind": state.total.kind.value,
                "transformed": bool(state.total.transformed),
                "established_by": state.total.established_by,
            }
        ),
    }
