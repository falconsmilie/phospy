"""Dataset-build preprocessing config validation."""

from __future__ import annotations

from phospy.contracts.configs import (
    DATASET_INTENSITY_TRANSFORM_POLICY_LOG2,
    DATASET_MISSING_DATA_POLICY_IMPUTE_MINPROB,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL,
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetLocalisationConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DatasetTotalProteinCorrectionConfig,
)
from phospy.errors.input import PhosPyInputError
from phospy.validation.configs.preprocessing import (
    validate_preprocessing_section_type,
)


class DatasetPreprocessingConfigValidator:
    """Validate the public config interpreted into the internal preprocessing plan."""

    def run(self, config: DatasetPreprocessingConfig) -> DatasetPreprocessingConfig:
        self._validate_intensity_transform(config.intensity_transform)
        self._validate_normalisation(config.normalisation)
        self._validate_missing_data(config.missing_data)
        self._validate_total_protein_correction(config.total_protein_correction)
        self._validate_site_matrix(config.site_matrix)
        self._validate_comparisons(config.comparisons)
        self._validate_localisation(config.localisation)
        self._validate_total_protein_correction_scale_contract(config)
        self._validate_minprob_scale_contract(config)
        return config

    def _validate_intensity_transform(
        self, config: DatasetIntensityTransformConfig
    ) -> None:
        validate_preprocessing_section_type(
            config,
            field_name="dataset build request preprocessing_config.intensity_transform",
            expected_type=DatasetIntensityTransformConfig,
        )

    def _validate_normalisation(self, config: DatasetNormalisationConfig) -> None:
        validate_preprocessing_section_type(
            config,
            field_name="dataset build request preprocessing_config.normalisation",
            expected_type=DatasetNormalisationConfig,
        )

    def _validate_missing_data(self, config: DatasetMissingDataConfig) -> None:
        validate_preprocessing_section_type(
            config,
            field_name="dataset build request preprocessing_config.missing_data",
            expected_type=DatasetMissingDataConfig,
        )

    def _validate_total_protein_correction(
        self, config: DatasetTotalProteinCorrectionConfig
    ) -> None:
        validate_preprocessing_section_type(
            config,
            field_name=(
                "dataset build request preprocessing_config.total_protein_correction"
            ),
            expected_type=DatasetTotalProteinCorrectionConfig,
        )

    def _validate_site_matrix(self, config: DatasetSiteMatrixConfig) -> None:
        validate_preprocessing_section_type(
            config,
            field_name="dataset build request preprocessing_config.site_matrix",
            expected_type=DatasetSiteMatrixConfig,
        )

    def _validate_comparisons(self, config: DatasetComparisonBuildingConfig) -> None:
        validate_preprocessing_section_type(
            config,
            field_name="dataset build request preprocessing_config.comparisons",
            expected_type=DatasetComparisonBuildingConfig,
        )

    def _validate_localisation(self, config: DatasetLocalisationConfig) -> None:
        validate_preprocessing_section_type(
            config,
            field_name="dataset build request preprocessing_config.localisation",
            expected_type=DatasetLocalisationConfig,
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

    def _validate_minprob_scale_contract(
        self,
        config: DatasetPreprocessingConfig,
    ) -> None:
        if config.missing_data.policy != DATASET_MISSING_DATA_POLICY_IMPUTE_MINPROB:
            return
        if config.intensity_transform.policy != DATASET_INTENSITY_TRANSFORM_POLICY_LOG2:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data.policy="
                "'impute_minprob' requires "
                "preprocessing_config.intensity_transform.policy='log2'. "
                "Set intensity_transform.policy='log2' or choose a different "
                "missing_data policy."
            )
