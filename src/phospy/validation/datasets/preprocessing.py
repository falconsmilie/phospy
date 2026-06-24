"""Dataset-build preprocessing config and plan validation."""

from __future__ import annotations

from phospy.contracts.configs import (
    DATASET_INTENSITY_TRANSFORM_POLICY_LOG2,
    DATASET_MISSING_DATA_POLICY_IMPUTE_MINPROB,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL,
    DatasetBatchCorrectionConfig,
    DatasetComparisonBuildingConfig,
    DatasetGroupCoverageFilterConfig,
    DatasetIntensityTransformConfig,
    DatasetLocalisationConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DatasetProteinAwarePreparationConfig,
    DatasetSiteMatrixConfig,
    DatasetTotalProteinCorrectionConfig,
)
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
)
from phospy.validation.configs.preprocessing import (
    validate_preprocessing_section_type,
)

_BATCH_CORRECTION_DOWNSTREAM_BOUNDARY_STAGES = (
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
)


class DatasetPreprocessingConfigValidator:
    """Validate the public config interpreted into the internal preprocessing plan."""

    def run(self, config: DatasetPreprocessingConfig) -> DatasetPreprocessingConfig:
        self._validate_intensity_transform(config.intensity_transform)
        self._validate_normalisation(config.normalisation)
        self._validate_missing_data(config.missing_data)
        self._validate_group_coverage_filter(config.group_coverage_filter)
        self._validate_total_protein_correction(config.total_protein_correction)
        self._validate_site_matrix(config.site_matrix)
        self._validate_comparisons(config.comparisons)
        self._validate_localisation(config.localisation)
        self._validate_batch_correction(config.batch_correction)
        self._validate_protein_aware_preparation(config.protein_aware_preparation)
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

    def _validate_group_coverage_filter(
        self, config: DatasetGroupCoverageFilterConfig
    ) -> None:
        validate_preprocessing_section_type(
            config,
            field_name=(
                "dataset build request preprocessing_config.group_coverage_filter"
            ),
            expected_type=DatasetGroupCoverageFilterConfig,
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

    def _validate_batch_correction(self, config: DatasetBatchCorrectionConfig) -> None:
        validate_preprocessing_section_type(
            config,
            field_name="dataset build request preprocessing_config.batch_correction",
            expected_type=DatasetBatchCorrectionConfig,
        )

    def _validate_protein_aware_preparation(
        self, config: DatasetProteinAwarePreparationConfig
    ) -> None:
        validate_preprocessing_section_type(
            config,
            field_name=(
                "dataset build request preprocessing_config.protein_aware_preparation"
            ),
            expected_type=DatasetProteinAwarePreparationConfig,
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


class PreprocessingStageOrderValidator:
    """Validate scientific preprocessing stage-order constraints."""

    def run(
        self,
        *,
        stage_order: tuple[str, ...],
        batch_correction_requested: bool,
    ) -> None:
        stages = tuple(str(stage).strip() for stage in stage_order)
        blank_positions = [
            position for position, stage in enumerate(stages) if stage == ""
        ]
        if blank_positions:
            raise PhosPyInputError(
                "dataset preprocessing plan stage_order contains blank stage "
                f"entries at positions {blank_positions}"
            )
        duplicates = [
            stage for stage in dict.fromkeys(stages) if stages.count(stage) > 1
        ]
        if duplicates:
            raise PhosPyInputError(
                "dataset preprocessing plan stage_order contains duplicate stages: "
                + ", ".join(repr(stage) for stage in duplicates)
            )
        if not batch_correction_requested:
            return
        if DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION not in stages:
            raise PhosPyInputError(
                "dataset preprocessing plan requests batch correction but "
                "stage_order does not include 'batch_correction'. Build plans "
                "from DatasetPreprocessingConfig or include the batch_correction "
                "stage explicitly."
            )

        batch_position = stages.index(DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION)
        if (
            DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM in stages
            and batch_position
            < stages.index(DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM)
        ):
            raise PhosPyInputError(
                "dataset preprocessing plan has unsupported stage_order: "
                "batch_correction must run after intensity_transform when both "
                "stages are configured"
            )

        downstream_before_batch = [
            stage
            for stage in _BATCH_CORRECTION_DOWNSTREAM_BOUNDARY_STAGES
            if stage in stages and stages.index(stage) < batch_position
        ]
        if downstream_before_batch:
            raise PhosPyInputError(
                "dataset preprocessing plan has unsupported stage_order: "
                "batch_correction cannot run after downstream stages have consumed "
                "the matrix because that would weaken the analysis-ready dataset "
                "boundary; downstream stages before batch_correction: "
                + ", ".join(downstream_before_batch)
            )
