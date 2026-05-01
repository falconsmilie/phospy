"""Dataset-level configuration aggregations."""
# pyright: reportUnnecessaryIsInstance=false
# Runtime boundary guards are intentionally retained for untyped external callers.

from __future__ import annotations

from dataclasses import dataclass, field

from phospy.api.configs.preprocessing import (
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


__all__ = [
    "DatasetPreprocessingConfig",
]
