"""Normalisation stage for dataset preprocessing."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd

from phospy.api.configs import (
    DATASET_NORMALISATION_POLICY_MEDIAN_CENTER,
    DATASET_NORMALISATION_POLICY_NONE,
    DATASET_NORMALISATION_POLICY_QUANTILE,
)
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    PreprocessingState,
)
from phospy.errors.input import PhosPyInputError


class NormalisationStage:
    """Apply configured normalisation policy to phosphosite sample columns."""

    stage_key = DATASET_PREPROCESSING_STAGE_NORMALISATION

    def run(self, state: PreprocessingState) -> PreprocessingState:
        policy = state.plan.normalisation_policy
        if policy == DATASET_NORMALISATION_POLICY_NONE:
            return state

        _require_non_empty_matrix(state.phospho, policy_name=policy)
        _require_numeric_columns(state.phospho, policy_name=policy)

        if policy == DATASET_NORMALISATION_POLICY_MEDIAN_CENTER:
            normalised = _median_center(state.phospho)
            return replace(state, phospho=normalised)
        if policy == DATASET_NORMALISATION_POLICY_QUANTILE:
            normalised = _quantile_normalise(state.phospho)
            return replace(state, phospho=normalised)
        raise PhosPyInputError(
            "dataset build request preprocessing_config contains an unsupported "
            "normalisation.policy"
        )


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
    rank_means = pd.DataFrame(sorted_values).mean(axis=1, skipna=True).to_numpy()

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


__all__ = ["NormalisationStage"]
