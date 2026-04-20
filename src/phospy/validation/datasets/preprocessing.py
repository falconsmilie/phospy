"""Dataset-build preprocessing config validation."""

from __future__ import annotations

from phospy.api.configs import (
    DATASET_MISSING_DATA_POLICIES,
    DATASET_MISSING_DATA_POLICY_FORBID,
    DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
    DatasetPreprocessingConfig,
)
from phospy.errors.input import PhosPyInputError


class DatasetPreprocessingConfigValidator:
    """Validate the supported public dataset preprocessing config."""

    def run(self, config: DatasetPreprocessingConfig) -> DatasetPreprocessingConfig:
        if not isinstance(config, DatasetPreprocessingConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config must be a "
                "DatasetPreprocessingConfig"
            )

        policy = config.missing_data_policy
        if policy not in DATASET_MISSING_DATA_POLICIES:
            supported = ", ".join(sorted(DATASET_MISSING_DATA_POLICIES))
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data_policy "
                f"must be one of: {supported}"
            )

        min_observed_values = config.min_observed_values
        if policy == DATASET_MISSING_DATA_POLICY_FORBID:
            if min_observed_values is not None:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.min_observed_values "
                    "must be None when missing_data_policy='forbid'"
                )
            return config

        if policy == DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN:
            if not isinstance(min_observed_values, int):
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.min_observed_values "
                    "must be an int when missing_data_policy='impute_row_median'"
                )
            if min_observed_values < 1:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.min_observed_values "
                    "must be greater than or equal to 1 when "
                    "missing_data_policy='impute_row_median'"
                )
            return config

        raise PhosPyInputError(
            "dataset build request preprocessing_config contains an unsupported "
            "missing_data_policy"
        )
