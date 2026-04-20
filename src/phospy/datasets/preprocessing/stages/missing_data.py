"""Missing-data preprocessing stage for dataset building."""

from __future__ import annotations

from dataclasses import replace

from phospy.api.configs import (
    DATASET_MISSING_DATA_POLICY_FORBID,
    DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
)
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    PreprocessingState,
)
from phospy.errors.input import PhosPyInputError


class MissingDataStage:
    """Apply the configured missing-data policy to phospho/site tables."""

    stage_key = DATASET_PREPROCESSING_STAGE_MISSING_DATA

    def run(self, state: PreprocessingState) -> PreprocessingState:
        if state.plan.missing_data_policy == DATASET_MISSING_DATA_POLICY_FORBID:
            return state

        if (
            state.plan.missing_data_policy
            != DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config contains an unsupported "
                "missing_data_policy"
            )

        min_observed_values = state.plan.min_observed_values
        if not isinstance(min_observed_values, int):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.min_observed_values "
                "must be an int when missing_data_policy='impute_row_median'"
            )
        if min_observed_values > state.phospho.shape[1]:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.min_observed_values "
                f"({min_observed_values}) cannot exceed the number of phospho "
                f"samples ({state.phospho.shape[1]})"
            )

        observed_counts = state.phospho.notna().sum(axis=1)
        retained_mask = observed_counts >= min_observed_values
        filtered_phospho = state.phospho.loc[retained_mask]
        filtered_site_metadata = state.site_metadata.loc[filtered_phospho.index]

        row_medians = filtered_phospho.median(axis=1, skipna=True)
        imputed = filtered_phospho.T.fillna(row_medians).T
        return replace(
            state,
            phospho=imputed,
            site_metadata=filtered_site_metadata,
        )


__all__ = ["MissingDataStage"]
