"""Normalisation preprocessing policy configuration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from phospy.contracts.configs.preprocessing._validation import (
    validate_normalisation_config,
)

DATASET_NORMALISATION_POLICY_NONE = "none"
DATASET_NORMALISATION_POLICY_MEDIAN_CENTER = "median_center"
DATASET_NORMALISATION_POLICY_QUANTILE = "quantile"
DatasetNormalisationPolicy = Literal["none", "median_center", "quantile"]
DATASET_NORMALISATION_POLICIES = frozenset(
    {
        DATASET_NORMALISATION_POLICY_NONE,
        DATASET_NORMALISATION_POLICY_MEDIAN_CENTER,
        DATASET_NORMALISATION_POLICY_QUANTILE,
    }
)


@dataclass(frozen=True, slots=True)
class DatasetNormalisationConfig:
    """Public normalisation policy options for dataset building.

    - `"none"`: no normalisation (strict default).
    - `"median_center"`: subtract sample-wise medians.
    - `"quantile"`: force sample columns to share one empirical distribution.
    """

    policy: DatasetNormalisationPolicy = DATASET_NORMALISATION_POLICY_NONE

    def __post_init__(self) -> None:
        validate_normalisation_config(
            policy=self.policy,
            supported_policies=DATASET_NORMALISATION_POLICIES,
        )


__all__ = [
    "DATASET_NORMALISATION_POLICIES",
    "DATASET_NORMALISATION_POLICY_MEDIAN_CENTER",
    "DATASET_NORMALISATION_POLICY_NONE",
    "DATASET_NORMALISATION_POLICY_QUANTILE",
    "DatasetNormalisationConfig",
    "DatasetNormalisationPolicy",
]
