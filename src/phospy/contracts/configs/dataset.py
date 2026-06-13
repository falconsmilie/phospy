"""Dataset-level configuration aggregations."""
# pyright: reportUnnecessaryIsInstance=false
# Runtime boundary guards are intentionally retained for untyped external callers.

from __future__ import annotations

from dataclasses import dataclass, field

from phospy.contracts.configs.preprocessing import (
    DATASET_COMPARISON_BUILDING_POLICY_NONE,
    DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY,
    DATASET_INTENSITY_TRANSFORM_POLICY_LOG2,
    DATASET_MISSING_DATA_POLICY_FORBID,
    DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
    DATASET_NORMALISATION_POLICY_MEDIAN_CENTER,
    DATASET_NORMALISATION_POLICY_NONE,
    DATASET_SITE_MATRIX_POLICY_AS_INPUT,
    DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
    DatasetBatchCorrectionConfig,
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetLocalisationConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetProteinAwarePreparationConfig,
    DatasetRuvReadinessConfig,
    DatasetSiteMatrixConfig,
    DatasetSiteSequenceResolutionConfig,
    DatasetTotalProteinCorrectionConfig,
)
from phospy.validation.configs.preprocessing import (
    validate_preprocessing_section_type,
)


@dataclass(frozen=True, slots=True)
class DatasetPreprocessingConfig:
    """Public grouped preprocessing policy for dataset building.

    The builder owns this policy surface. Groups are intentionally separated so
    supported preprocessing science remains user-visible:

    - `intensity_transform`: quantitative transform policy.
    - `normalisation`: sample-wise normalisation policy.
    - `missing_data`: missing-value handling policy.
    - `total_protein_correction`: total/protein correction policy.
    - `protein_aware_preparation`: config-only intent for future aligned
      phosphosite/protein model input preparation. It does not run preparation
      or modelling during dataset build.
    - `site_matrix`: site-matrix construction policy.
    - `site_sequence_resolution`: optional local FASTA-backed site-sequence
      resolution policy.
    - `comparisons`: comparison-building policy.
    - `localisation`: phosphosite-localisation eligibility policy.
    - `batch_correction`: batch-correction intent declaration. Currently
      config-only; no correction is executed by the dataset builder.
    - `ruv_readiness`: readiness reporting contract for future RUV-compatible
      preprocessing (report-only; no correction).
    """

    intensity_transform: DatasetIntensityTransformConfig = field(
        default_factory=DatasetIntensityTransformConfig
    )
    normalisation: DatasetNormalisationConfig = field(
        default_factory=DatasetNormalisationConfig
    )
    missing_data: DatasetMissingDataConfig = field(
        default_factory=DatasetMissingDataConfig
    )
    total_protein_correction: DatasetTotalProteinCorrectionConfig = field(
        default_factory=DatasetTotalProteinCorrectionConfig
    )
    site_matrix: DatasetSiteMatrixConfig = field(
        default_factory=DatasetSiteMatrixConfig
    )
    site_sequence_resolution: DatasetSiteSequenceResolutionConfig = field(
        default_factory=DatasetSiteSequenceResolutionConfig
    )
    comparisons: DatasetComparisonBuildingConfig = field(
        default_factory=DatasetComparisonBuildingConfig
    )
    localisation: DatasetLocalisationConfig = field(
        default_factory=DatasetLocalisationConfig
    )
    batch_correction: DatasetBatchCorrectionConfig = field(
        default_factory=DatasetBatchCorrectionConfig
    )
    ruv_readiness: DatasetRuvReadinessConfig = field(
        default_factory=DatasetRuvReadinessConfig
    )
    protein_aware_preparation: DatasetProteinAwarePreparationConfig = field(
        default_factory=DatasetProteinAwarePreparationConfig
    )

    def __post_init__(self) -> None:
        validate_preprocessing_section_type(
            self.intensity_transform,
            field_name=(
                "dataset build request preprocessing_config.intensity_transform"
            ),
            expected_type=DatasetIntensityTransformConfig,
        )
        validate_preprocessing_section_type(
            self.normalisation,
            field_name="dataset build request preprocessing_config.normalisation",
            expected_type=DatasetNormalisationConfig,
        )
        validate_preprocessing_section_type(
            self.missing_data,
            field_name="dataset build request preprocessing_config.missing_data",
            expected_type=DatasetMissingDataConfig,
        )
        validate_preprocessing_section_type(
            self.total_protein_correction,
            field_name=(
                "dataset build request preprocessing_config.total_protein_correction"
            ),
            expected_type=DatasetTotalProteinCorrectionConfig,
        )
        validate_preprocessing_section_type(
            self.site_matrix,
            field_name="dataset build request preprocessing_config.site_matrix",
            expected_type=DatasetSiteMatrixConfig,
        )
        validate_preprocessing_section_type(
            self.site_sequence_resolution,
            field_name=(
                "dataset build request preprocessing_config.site_sequence_resolution"
            ),
            expected_type=DatasetSiteSequenceResolutionConfig,
        )
        validate_preprocessing_section_type(
            self.comparisons,
            field_name="dataset build request preprocessing_config.comparisons",
            expected_type=DatasetComparisonBuildingConfig,
        )
        validate_preprocessing_section_type(
            self.localisation,
            field_name="dataset build request preprocessing_config.localisation",
            expected_type=DatasetLocalisationConfig,
        )
        validate_preprocessing_section_type(
            self.batch_correction,
            field_name="dataset build request preprocessing_config.batch_correction",
            expected_type=DatasetBatchCorrectionConfig,
        )
        validate_preprocessing_section_type(
            self.ruv_readiness,
            field_name="dataset build request preprocessing_config.ruv_readiness",
            expected_type=DatasetRuvReadinessConfig,
        )
        validate_preprocessing_section_type(
            self.protein_aware_preparation,
            field_name=(
                "dataset build request preprocessing_config.protein_aware_preparation"
            ),
            expected_type=DatasetProteinAwarePreparationConfig,
        )

    @classmethod
    def default(cls) -> DatasetPreprocessingConfig:
        """Return the package default preprocessing profile."""
        return cls()

    @classmethod
    def strict(cls) -> DatasetPreprocessingConfig:
        """Return strict preprocessing with no imputation or implicit transforms."""
        return cls(
            intensity_transform=DatasetIntensityTransformConfig(
                policy=DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY
            ),
            normalisation=DatasetNormalisationConfig(
                policy=DATASET_NORMALISATION_POLICY_NONE
            ),
            missing_data=DatasetMissingDataConfig(
                policy=DATASET_MISSING_DATA_POLICY_FORBID
            ),
            total_protein_correction=DatasetTotalProteinCorrectionConfig(
                policy=DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE
            ),
            site_matrix=DatasetSiteMatrixConfig(
                policy=DATASET_SITE_MATRIX_POLICY_AS_INPUT
            ),
            comparisons=DatasetComparisonBuildingConfig(
                policy=DATASET_COMPARISON_BUILDING_POLICY_NONE
            ),
        )

    @classmethod
    def from_raw_phosphosite_table(cls) -> DatasetPreprocessingConfig:
        """Return a practical profile for common raw phosphosite matrices."""
        return cls(
            intensity_transform=DatasetIntensityTransformConfig(
                policy=DATASET_INTENSITY_TRANSFORM_POLICY_LOG2
            ),
            normalisation=DatasetNormalisationConfig(
                policy=DATASET_NORMALISATION_POLICY_MEDIAN_CENTER
            ),
            missing_data=DatasetMissingDataConfig(
                policy=DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
                min_observed_values=1,
            ),
            total_protein_correction=DatasetTotalProteinCorrectionConfig(
                policy=DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE
            ),
            site_matrix=DatasetSiteMatrixConfig(
                policy=DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA
            ),
            comparisons=DatasetComparisonBuildingConfig(
                policy=DATASET_COMPARISON_BUILDING_POLICY_NONE
            ),
        )


__all__ = [
    "DatasetPreprocessingConfig",
]
