"""Internal dataset preprocessing planning and state models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import pandas as pd

from phospy.api.configs import (
    DATASET_COMPARISON_BUILDING_DEFAULT_SAMPLE_GROUP_COLUMN,
    DATASET_COMPARISON_BUILDING_POLICY_NONE,
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_MAX_MEAN_SIGNAL,
    DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
    DATASET_SITE_MATRIX_POLICY_AS_INPUT,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
    DatasetComparisonBuildingPolicy,
    DatasetComparisonPair,
    DatasetMissingDataPolicy,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixDuplicateSiteStrategy,
    DatasetSiteMatrixMissingDataPolicy,
    DatasetSiteMatrixPolicy,
    DatasetTotalProteinCorrectionPolicy,
)

DATASET_PREPROCESSING_STAGE_MISSING_DATA = "missing_data"
DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION = "total_protein_correction"
DATASET_PREPROCESSING_STAGE_SITE_MATRIX = "site_matrix"
DATASET_PREPROCESSING_STAGE_COMPARISONS = "comparisons"
DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT = (DATASET_PREPROCESSING_STAGE_MISSING_DATA,)


@dataclass(frozen=True, slots=True)
class PreprocessingPlan:
    """Execution-ready internal preprocessing plan derived from public config."""

    missing_data_policy: DatasetMissingDataPolicy
    missing_data_min_observed_values: int | None
    total_protein_correction_policy: DatasetTotalProteinCorrectionPolicy
    site_matrix_policy: DatasetSiteMatrixPolicy
    comparison_building_policy: DatasetComparisonBuildingPolicy = (
        DATASET_COMPARISON_BUILDING_POLICY_NONE
    )
    site_matrix_duplicate_site_strategy: DatasetSiteMatrixDuplicateSiteStrategy = (
        DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_MAX_MEAN_SIGNAL
    )
    site_matrix_missing_data_policy: DatasetSiteMatrixMissingDataPolicy = (
        DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING
    )
    site_matrix_minimum_observed_values: int | None = None
    comparison_sample_group_column: str = (
        DATASET_COMPARISON_BUILDING_DEFAULT_SAMPLE_GROUP_COLUMN
    )
    comparison_pairs: tuple[DatasetComparisonPair, ...] | None = None
    stage_order: tuple[str, ...] = DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT

    @classmethod
    def from_config(cls, config: DatasetPreprocessingConfig) -> PreprocessingPlan:
        stage_order: list[str] = [DATASET_PREPROCESSING_STAGE_MISSING_DATA]
        if (
            config.total_protein_correction.policy
            != DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE
        ):
            stage_order.append(DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION)
        if config.site_matrix.policy != DATASET_SITE_MATRIX_POLICY_AS_INPUT:
            stage_order.append(DATASET_PREPROCESSING_STAGE_SITE_MATRIX)
        if config.comparisons.policy != DATASET_COMPARISON_BUILDING_POLICY_NONE:
            stage_order.append(DATASET_PREPROCESSING_STAGE_COMPARISONS)
        return cls(
            missing_data_policy=config.missing_data.policy,
            missing_data_min_observed_values=config.missing_data.min_observed_values,
            total_protein_correction_policy=config.total_protein_correction.policy,
            site_matrix_policy=config.site_matrix.policy,
            site_matrix_duplicate_site_strategy=config.site_matrix.duplicate_site_strategy,
            site_matrix_missing_data_policy=config.site_matrix.missing_data_policy,
            site_matrix_minimum_observed_values=config.site_matrix.minimum_observed_values,
            comparison_building_policy=config.comparisons.policy,
            comparison_sample_group_column=config.comparisons.sample_group_column,
            comparison_pairs=(
                None
                if config.comparisons.pairs is None
                else tuple(config.comparisons.pairs)
            ),
            stage_order=tuple(stage_order),
        )

    @classmethod
    def default(cls) -> PreprocessingPlan:
        return cls.from_config(DatasetPreprocessingConfig())


@dataclass(frozen=True, slots=True)
class PreprocessingState:
    """Internal preprocessing state carried between ordered stages."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None
    total: pd.DataFrame | None
    plan: PreprocessingPlan
    comparisons: pd.DataFrame | None = None


class PreprocessingStage(Protocol):
    """Single internal preprocessing stage contract."""

    stage_key: str

    def run(self, state: PreprocessingState) -> PreprocessingState:
        """Apply a preprocessing stage and return the next state."""


__all__ = [
    "DATASET_PREPROCESSING_STAGE_COMPARISONS",
    "DATASET_PREPROCESSING_STAGE_MISSING_DATA",
    "DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT",
    "DATASET_PREPROCESSING_STAGE_SITE_MATRIX",
    "DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION",
    "PreprocessingPlan",
    "PreprocessingStage",
    "PreprocessingState",
]
