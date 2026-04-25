"""Dataset-build preprocessing config validation."""

from __future__ import annotations

from phospy.api.configs import (
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
