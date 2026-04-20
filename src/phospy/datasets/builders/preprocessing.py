"""Dataset-builder preprocessing application."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.api.configs import (
    DATASET_MISSING_DATA_POLICY_FORBID,
    DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
    DatasetPreprocessingConfig,
)
from phospy.errors.input import PhosPyInputError


@dataclass(frozen=True, slots=True)
class PreprocessedDatasetTables:
    """Preprocessed builder tables prior to transformation-state establishment."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame


class DatasetPreprocessor:
    """Apply supported public preprocessing policies for dataset building."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        config: DatasetPreprocessingConfig,
    ) -> PreprocessedDatasetTables:
        if config.missing_data_policy == DATASET_MISSING_DATA_POLICY_FORBID:
            return PreprocessedDatasetTables(
                phospho=phospho,
                site_metadata=site_metadata,
            )
        if config.missing_data_policy == DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN:
            return self._run_row_median_imputation(
                phospho=phospho,
                site_metadata=site_metadata,
                min_observed_values=config.min_observed_values,
            )
        raise PhosPyInputError(
            "dataset build request preprocessing_config contains an unsupported "
            "missing_data_policy"
        )

    def _run_row_median_imputation(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        min_observed_values: int | None,
    ) -> PreprocessedDatasetTables:
        if not isinstance(min_observed_values, int):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.min_observed_values "
                "must be an int when missing_data_policy='impute_row_median'"
            )
        if min_observed_values > phospho.shape[1]:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.min_observed_values "
                f"({min_observed_values}) cannot exceed the number of phospho "
                f"samples ({phospho.shape[1]})"
            )

        observed_counts = phospho.notna().sum(axis=1)
        retained_mask = observed_counts >= min_observed_values
        filtered_phospho = phospho.loc[retained_mask]
        filtered_site_metadata = site_metadata.loc[filtered_phospho.index]

        row_medians = filtered_phospho.median(axis=1, skipna=True)
        imputed = filtered_phospho.T.fillna(row_medians).T
        return PreprocessedDatasetTables(
            phospho=imputed,
            site_metadata=filtered_site_metadata,
        )
