"""Internal dataset preprocessing planning and state models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import pandas as pd

from phospy.api.configs import DatasetMissingDataPolicy, DatasetPreprocessingConfig

DATASET_PREPROCESSING_STAGE_MISSING_DATA = "missing_data"
DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION = "total_protein_correction"
DATASET_PREPROCESSING_STAGE_SITE_MATRIX = "site_matrix"
DATASET_PREPROCESSING_STAGE_COMPARISONS = "comparisons"
DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT = (DATASET_PREPROCESSING_STAGE_MISSING_DATA,)


@dataclass(frozen=True, slots=True)
class PreprocessingPlan:
    """Execution-ready internal preprocessing plan derived from public config."""

    missing_data_policy: DatasetMissingDataPolicy
    min_observed_values: int | None
    stage_order: tuple[str, ...] = DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT

    @classmethod
    def from_config(cls, config: DatasetPreprocessingConfig) -> PreprocessingPlan:
        return cls(
            missing_data_policy=config.missing_data_policy,
            min_observed_values=config.min_observed_values,
        )


@dataclass(frozen=True, slots=True)
class PreprocessingState:
    """Internal preprocessing state carried between ordered stages."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None
    total: pd.DataFrame | None
    config: DatasetPreprocessingConfig
    plan: PreprocessingPlan


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
