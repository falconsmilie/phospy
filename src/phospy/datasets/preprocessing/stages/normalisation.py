"""Normalisation stage for dataset preprocessing."""

from __future__ import annotations

from dataclasses import replace
from numbers import Real
from typing import TypedDict

import numpy as np
import pandas as pd

from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
    PreprocessingStateTableKey,
)
from phospy.datasets.preprocessing.policy_models import NormalisationPolicy
from phospy.datasets.preprocessing.stage_contract import PreprocessingStageContract
from phospy.errors.input import PhosPyInputError
from phospy.provenance.hashing import hash_table


class _ShapePayload(TypedDict):
    rows: int
    columns: int


class _SampleSummary(TypedDict):
    non_missing_count: int
    missing_count: int
    mean: float | None
    median: float | None
    std: float | None
    min: float | None
    max: float | None


class _NormalisationDiagnostics(TypedDict):
    method: str
    parameters: dict[str, object]
    policy: str
    affected_columns: list[str]
    input_matrix_shape: _ShapePayload
    output_matrix_shape: _ShapePayload
    per_sample_summary_before: dict[str, _SampleSummary]
    per_sample_summary_after: dict[str, _SampleSummary]
    rows_dropped: bool
    columns_dropped: bool
    dropped_row_ids: tuple[str, ...]
    dropped_column_ids: tuple[str, ...]
    dropped_row_count: int
    dropped_column_count: int
    input_phospho_hash: str
    output_phospho_hash: str
    note: str | None


class NormalisationStage:
    """Apply configured normalisation policy to phosphosite sample columns."""

    stage_key = DATASET_PREPROCESSING_STAGE_NORMALISATION

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        policy = state.plan.normalisation_policy
        method_parameters = _resolve_method_parameters(policy)
        if policy is NormalisationPolicy.NONE:
            diagnostics = _build_diagnostics(
                policy=policy,
                method_parameters=method_parameters,
                before=state.phospho,
                after=state.phospho,
            )
            return PreprocessingStageResult(
                state=state,
                diagnostics={
                    "dropped_row_ids": diagnostics["dropped_row_ids"],
                    "dropped_row_count": diagnostics["dropped_row_count"],
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "stage executed",
                    "diagnostics": diagnostics,
                },
            )

        _require_non_empty_matrix(state.phospho, policy_name=policy.value)
        _require_numeric_columns(state.phospho, policy_name=policy.value)

        if policy is NormalisationPolicy.MEDIAN_CENTER:
            normalised = _median_center(state.phospho)
            next_state = replace(state, phospho=normalised)
            diagnostics = _build_diagnostics(
                policy=policy,
                method_parameters=method_parameters,
                before=state.phospho,
                after=normalised,
            )
            return PreprocessingStageResult(
                state=next_state,
                diagnostics={
                    "dropped_row_ids": diagnostics["dropped_row_ids"],
                    "dropped_row_count": diagnostics["dropped_row_count"],
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "stage executed",
                    "diagnostics": diagnostics,
                },
            )
        if policy is NormalisationPolicy.QUANTILE:
            normalised = _quantile_normalise(state.phospho)
            next_state = replace(state, phospho=normalised)
            diagnostics = _build_diagnostics(
                policy=policy,
                method_parameters=method_parameters,
                before=state.phospho,
                after=normalised,
            )
            return PreprocessingStageResult(
                state=next_state,
                diagnostics={
                    "dropped_row_ids": diagnostics["dropped_row_ids"],
                    "dropped_row_count": diagnostics["dropped_row_count"],
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "stage executed",
                    "diagnostics": diagnostics,
                },
            )
        raise PhosPyInputError(
            "dataset build request preprocessing_config contains an unsupported "
            "normalisation.policy"
        )


def _build_diagnostics(
    *,
    policy: NormalisationPolicy,
    method_parameters: dict[str, object],
    before: pd.DataFrame,
    after: pd.DataFrame,
) -> _NormalisationDiagnostics:
    dropped_row_ids = _resolve_dropped_labels(before.index, after.index)
    dropped_column_ids = _resolve_dropped_labels(before.columns, after.columns)
    diagnostics: _NormalisationDiagnostics = {
        "method": policy.value,
        "parameters": dict(method_parameters),
        "policy": policy.value,
        "affected_columns": [str(column) for column in after.columns.tolist()],
        "input_matrix_shape": _shape_payload(before),
        "output_matrix_shape": _shape_payload(after),
        "per_sample_summary_before": _per_sample_summary(before),
        "per_sample_summary_after": _per_sample_summary(after),
        "rows_dropped": bool(dropped_row_ids),
        "columns_dropped": bool(dropped_column_ids),
        "dropped_row_ids": dropped_row_ids,
        "dropped_column_ids": dropped_column_ids,
        "dropped_row_count": len(dropped_row_ids),
        "dropped_column_count": len(dropped_column_ids),
        "input_phospho_hash": hash_table(
            before,
            name="normalisation.input.phospho",
        ),
        "output_phospho_hash": hash_table(
            after,
            name="normalisation.output.phospho",
        ),
        "note": None,
    }
    if policy is not NormalisationPolicy.NONE:
        diagnostics["note"] = (
            "quantile normalisation used"
            if policy is NormalisationPolicy.QUANTILE
            else "median centering used"
        )
    return diagnostics


def _shape_payload(matrix: pd.DataFrame) -> _ShapePayload:
    return {
        "rows": int(matrix.shape[0]),
        "columns": int(matrix.shape[1]),
    }


def _resolve_dropped_labels(before: pd.Index, after: pd.Index) -> tuple[str, ...]:
    after_values = {str(value) for value in after.tolist()}
    return tuple(
        str(value) for value in before.tolist() if str(value) not in after_values
    )


def _safe_float(value: object) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool):
        return None
    if not isinstance(value, Real):
        return None
    resolved = float(value)
    if not np.isfinite(resolved):
        return None
    return resolved


def _per_sample_summary(matrix: pd.DataFrame) -> dict[str, _SampleSummary]:
    summary: dict[str, _SampleSummary] = {}
    for column in matrix.columns.tolist():
        series = matrix.loc[:, column]
        non_missing_count = int(series.notna().sum())
        missing_count = int(series.isna().sum())
        observed = series.dropna().astype("float64")
        summary[str(column)] = {
            "non_missing_count": non_missing_count,
            "missing_count": missing_count,
            "mean": _safe_float(observed.mean()) if non_missing_count > 0 else None,
            "median": (
                _safe_float(observed.median()) if non_missing_count > 0 else None
            ),
            "std": _safe_float(observed.std(ddof=1)) if non_missing_count > 0 else None,
            "min": _safe_float(observed.min()) if non_missing_count > 0 else None,
            "max": _safe_float(observed.max()) if non_missing_count > 0 else None,
        }
    return summary


def _require_non_empty_matrix(matrix: pd.DataFrame, *, policy_name: str) -> None:
    if matrix.empty:
        raise PhosPyInputError(
            "dataset build request preprocessing "
            f"normalisation.policy='{policy_name}' requires a non-empty phospho "
            "matrix"
        )


def _require_numeric_columns(matrix: pd.DataFrame, *, policy_name: str) -> None:
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
            f"normalisation.policy='{policy_name}' requires numeric phospho columns. "
            "Non-numeric columns: " + ", ".join(non_numeric)
        )


def _median_center(matrix: pd.DataFrame) -> pd.DataFrame:
    medians = matrix.median(axis=0, skipna=True)
    centered = matrix.subtract(medians, axis=1)
    centered.index = matrix.index.copy()
    centered.columns = matrix.columns.copy()
    return centered


def _quantile_normalise(matrix: pd.DataFrame) -> pd.DataFrame:
    """Quantile-normalise sample columns with deterministic tie handling.

    Performance note: this path is dense and sort-heavy. Runtime scales roughly
    with O(n_sites * n_samples * log(n_sites)) and requires additional float64
    matrix copies.
    """

    as_float = matrix.astype("float64")
    sorted_values = np.sort(as_float.to_numpy(copy=True), axis=0)
    rank_means = np.nanmean(sorted_values, axis=1)

    normalised_columns: dict[object, pd.Series] = {}
    for column in as_float.columns:
        series = as_float.loc[:, column]
        observed = series.dropna()
        if observed.empty:
            normalised_columns[column] = pd.Series(
                np.nan,
                index=series.index.copy(),
                dtype="float64",
            )
            continue
        sorted_observed = observed.sort_values(kind="mergesort")
        quantile_values = _assign_quantile_values_with_deterministic_ties(
            sorted_values=sorted_observed.to_numpy(copy=True),
            rank_means=rank_means[: len(sorted_observed)],
        )
        mapped = pd.Series(np.nan, index=series.index.copy(), dtype="float64")
        mapped.loc[sorted_observed.index] = quantile_values
        normalised_columns[column] = mapped

    normalised = pd.DataFrame(
        {column: normalised_columns[column] for column in matrix.columns},
        index=matrix.index.copy(),
    )
    normalised.columns = matrix.columns.copy()
    return normalised


def _assign_quantile_values_with_deterministic_ties(
    *,
    sorted_values: np.ndarray,
    rank_means: np.ndarray,
) -> np.ndarray:
    assigned = np.empty_like(sorted_values, dtype="float64")
    start = 0
    while start < sorted_values.shape[0]:
        stop = start + 1
        while (
            stop < sorted_values.shape[0]
            and sorted_values[stop] == sorted_values[start]
        ):
            stop += 1
        assigned[start:stop] = float(np.nanmean(rank_means[start:stop]))
        start = stop
    return assigned


def _resolve_operation(plan: PreprocessingPlan) -> str:
    return plan.normalisation_policy.value


def _resolve_method_parameters(policy: NormalisationPolicy) -> dict[str, object]:
    if policy is NormalisationPolicy.NONE:
        return {"applied": False}
    if policy is NormalisationPolicy.MEDIAN_CENTER:
        return {
            "applied": True,
            "centering_statistic": "median",
            "axis": "columns",
            "skipna": True,
        }
    if policy is NormalisationPolicy.QUANTILE:
        return {
            "applied": True,
            "target_distribution": "mean_rank_distribution",
            "tie_strategy": "deterministic_rank_average",
            "dtype": "float64",
        }
    raise PhosPyInputError(
        "dataset build request preprocessing_config contains an unsupported "
        "normalisation.policy"
    )


def _resolve_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    return _resolve_method_parameters(plan.normalisation_policy)


NORMALISATION_STAGE_CONTRACT = PreprocessingStageContract(
    stage_key=DATASET_PREPROCESSING_STAGE_NORMALISATION,
    display_label=DATASET_PREPROCESSING_STAGE_NORMALISATION,
    provenance_stage=DATASET_PREPROCESSING_STAGE_NORMALISATION,
    operation_name=_resolve_operation,
    serialize_parameters=_resolve_parameters,
    consumed_input_tables=(PreprocessingStateTableKey.DATASET_PHOSPHO,),
    produced_output_tables=(PreprocessingStateTableKey.DATASET_PHOSPHO,),
    stage_factory=NormalisationStage,
    backend="numpy",
    diagnostics_metadata={
        "known_diagnostics_fields": (
            "method",
            "parameters",
            "policy",
            "affected_columns",
            "input_matrix_shape",
            "output_matrix_shape",
            "per_sample_summary_before",
            "per_sample_summary_after",
            "rows_dropped",
            "columns_dropped",
            "dropped_row_ids",
            "dropped_column_ids",
            "dropped_row_count",
            "dropped_column_count",
            "input_phospho_hash",
            "output_phospho_hash",
            "note",
        )
    },
)


__all__ = ["NORMALISATION_STAGE_CONTRACT", "NormalisationStage"]
