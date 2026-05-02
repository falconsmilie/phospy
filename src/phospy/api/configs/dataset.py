"""Dataset-level configuration aggregations."""
# pyright: reportUnnecessaryIsInstance=false
# Runtime boundary guards are intentionally retained for untyped external callers.

from __future__ import annotations

from dataclasses import dataclass, field

from phospy.api.configs.preprocessing import (
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
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetSiteMatrixConfig,
    DatasetTotalProteinCorrectionConfig,
)
from phospy.errors.input import PhosPyInputError


@dataclass(frozen=True, slots=True)
class DatasetPreprocessingConfig:
    """Public grouped preprocessing policy for dataset building.

    The builder owns this policy surface. Groups are intentionally separated so
    supported preprocessing science remains user-visible:

    - `intensity_transform`: quantitative transform policy.
    - `normalisation`: sample-wise normalisation policy.
    - `missing_data`: missing-value handling policy.
    - `total_protein_correction`: total/protein correction policy.
    - `site_matrix`: site-matrix construction policy.
    - `comparisons`: comparison-building policy.
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
    comparisons: DatasetComparisonBuildingConfig = field(
        default_factory=DatasetComparisonBuildingConfig
    )

    def __post_init__(self) -> None:
        if not isinstance(self.intensity_transform, DatasetIntensityTransformConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.intensity_transform "
                "must be a DatasetIntensityTransformConfig"
            )
        if not isinstance(self.normalisation, DatasetNormalisationConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.normalisation must be a "
                "DatasetNormalisationConfig"
            )
        if not isinstance(self.missing_data, DatasetMissingDataConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data must be a "
                "DatasetMissingDataConfig"
            )
        if not isinstance(
            self.total_protein_correction, DatasetTotalProteinCorrectionConfig
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction "
                "must be a DatasetTotalProteinCorrectionConfig"
            )
        if not isinstance(self.site_matrix, DatasetSiteMatrixConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.site_matrix must be a "
                "DatasetSiteMatrixConfig"
            )
        if not isinstance(self.comparisons, DatasetComparisonBuildingConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons must be a "
                "DatasetComparisonBuildingConfig"
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
