"""Normalisation stage for dataset preprocessing."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd

from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
    PreprocessingStateTableKey,
)
from phospy.datasets.preprocessing.stage_contract import PreprocessingStageContract
from phospy.errors.input import PhosPyInputError
from phospy.policy_models import NormalisationPolicy
from phospy.provenance.hashing import hash_table


class NormalisationStage:
    """Apply configured normalisation policy to phosphosite sample columns."""

    stage_key = DATASET_PREPROCESSING_STAGE_NORMALISATION

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        policy = state.plan.normalisation_policy
        if policy is NormalisationPolicy.NONE:
            return PreprocessingStageResult(
                state=state,
                diagnostics={
                    "dropped_row_ids": (),
                    "dropped_row_count": 0,
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "stage executed",
                    "diagnostics": {
                        "policy": policy.value,
                        "affected_columns": [
                            str(column) for column in state.phospho.columns.tolist()
                        ],
                        "input_phospho_hash": hash_table(
                            state.phospho,
                            name="normalisation.input.phospho",
                        ),
                        "output_phospho_hash": hash_table(
                            state.phospho,
                            name="normalisation.output.phospho",
                        ),
                    },
                },
            )

        _require_non_empty_matrix(state.phospho, policy_name=policy.value)
        _require_numeric_columns(state.phospho, policy_name=policy.value)

        if policy is NormalisationPolicy.MEDIAN_CENTER:
            normalised = _median_center(state.phospho)
            next_state = replace(state, phospho=normalised)
            return PreprocessingStageResult(
                state=next_state,
                diagnostics={
                    "dropped_row_ids": (),
                    "dropped_row_count": 0,
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "stage executed",
                    "diagnostics": _build_diagnostics(
                        policy=policy,
                        before=state.phospho,
                        after=normalised,
                    ),
                },
            )
        if policy is NormalisationPolicy.QUANTILE:
            normalised = _quantile_normalise(state.phospho)
            next_state = replace(state, phospho=normalised)
            return PreprocessingStageResult(
                state=next_state,
                diagnostics={
                    "dropped_row_ids": (),
                    "dropped_row_count": 0,
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "stage executed",
                    "diagnostics": _build_diagnostics(
                        policy=policy,
                        before=state.phospho,
                        after=normalised,
                    ),
                },
            )
        raise PhosPyInputError(
            "dataset build request preprocessing_config contains an unsupported "
            "normalisation.policy"
        )


def _build_diagnostics(
    *,
    policy: NormalisationPolicy,
    before: pd.DataFrame,
    after: pd.DataFrame,
) -> dict[str, object]:
    diagnostics: dict[str, object] = {
        "policy": policy.value,
        "affected_columns": [str(column) for column in after.columns.tolist()],
        "input_phospho_hash": hash_table(
            before,
            name="normalisation.input.phospho",
        ),
        "output_phospho_hash": hash_table(
            after,
            name="normalisation.output.phospho",
        ),
    }
    if policy is not NormalisationPolicy.NONE:
        diagnostics["note"] = (
            "quantile normalisation used"
            if policy is NormalisationPolicy.QUANTILE
            else "median centering used"
        )
    return diagnostics


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


def _resolve_parameters(_plan: PreprocessingPlan) -> dict[str, object]:
    return {}


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
            "policy",
            "affected_columns",
            "input_phospho_hash",
            "output_phospho_hash",
            "note",
        )
    },
)


__all__ = ["NORMALISATION_STAGE_CONTRACT", "NormalisationStage"]
