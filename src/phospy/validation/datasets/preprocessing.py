"""Dataset-build preprocessing config validation."""

from __future__ import annotations

from phospy.api.configs import (
    DATASET_INTENSITY_TRANSFORM_POLICY_LOG2,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL,
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DatasetTotalProteinCorrectionConfig,
)
from phospy.errors.input import PhosPyInputError


class DatasetPreprocessingConfigValidator:
    """Validate the public config interpreted into the internal preprocessing plan."""

    def run(self, config: DatasetPreprocessingConfig) -> DatasetPreprocessingConfig:
        if not isinstance(config, DatasetPreprocessingConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config must be a "
                "DatasetPreprocessingConfig"
            )

        self._validate_intensity_transform(config.intensity_transform)
        self._validate_normalisation(config.normalisation)
        self._validate_missing_data(config.missing_data)
        self._validate_total_protein_correction(config.total_protein_correction)
        self._validate_site_matrix(config.site_matrix)
        self._validate_comparisons(config.comparisons)
        self._validate_total_protein_correction_scale_contract(config)
        return config

    def _validate_intensity_transform(
        self, config: DatasetIntensityTransformConfig
    ) -> None:
        if not isinstance(config, DatasetIntensityTransformConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.intensity_transform "
                "must be a DatasetIntensityTransformConfig"
            )

    def _validate_normalisation(self, config: DatasetNormalisationConfig) -> None:
        if not isinstance(config, DatasetNormalisationConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.normalisation must be a "
                "DatasetNormalisationConfig"
            )

    def _validate_missing_data(self, config: DatasetMissingDataConfig) -> None:
        if not isinstance(config, DatasetMissingDataConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data must be a "
                "DatasetMissingDataConfig"
            )

    def _validate_total_protein_correction(
        self, config: DatasetTotalProteinCorrectionConfig
    ) -> None:
        if not isinstance(config, DatasetTotalProteinCorrectionConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction "
                "must be a DatasetTotalProteinCorrectionConfig"
            )

    def _validate_site_matrix(self, config: DatasetSiteMatrixConfig) -> None:
        if not isinstance(config, DatasetSiteMatrixConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.site_matrix must be a "
                "DatasetSiteMatrixConfig"
            )

    def _validate_comparisons(self, config: DatasetComparisonBuildingConfig) -> None:
        if not isinstance(config, DatasetComparisonBuildingConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons must be a "
                "DatasetComparisonBuildingConfig"
            )

    def _validate_total_protein_correction_scale_contract(
        self,
        config: DatasetPreprocessingConfig,
    ) -> None:
        requested_policy = config.total_protein_correction.policy
        resolved_policy = requested_policy
        if resolved_policy == DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE:
            return
        if (
            resolved_policy
            == DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL
            and config.intensity_transform.policy
            != DATASET_INTENSITY_TRANSFORM_POLICY_LOG2
        ):
            raise PhosPyInputError(
                "dataset build request "
                f"preprocessing_config.total_protein_correction.policy={requested_policy!r} "
                "requires log2-scale phospho and total values. Configure "
                "preprocessing_config.intensity_transform.policy='log2', or disable "
                "total-protein correction."
            )
