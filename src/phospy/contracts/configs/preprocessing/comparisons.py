"""Comparison-building preprocessing policy configuration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from phospy.contracts.configs.preprocessing._validation import (
    validate_comparison_building_config,
)

DATASET_COMPARISON_BUILDING_POLICY_NONE = "none"
DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS = "sample_metadata_pairs"
DatasetComparisonBuildingPolicy = Literal["none", "sample_metadata_pairs"]
DatasetComparisonPair = tuple[str, str]
DATASET_COMPARISON_BUILDING_DEFAULT_SAMPLE_GROUP_COLUMN = "comparison_group"
DATASET_COMPARISON_BUILDING_POLICIES = frozenset(
    {
        DATASET_COMPARISON_BUILDING_POLICY_NONE,
        DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS,
    }
)


@dataclass(frozen=True, slots=True)
class DatasetComparisonBuildingConfig:
    """Public comparison-building policy options for dataset building.

    - `"none"`: do not build dataset-level pairwise comparisons.
    - `"sample_metadata_pairs"`: build comparison columns from grouped sample
      metadata.

    For `"sample_metadata_pairs"`:

    - `sample_group_column` must exist in `sample_metadata` and define one
      non-empty group label per sample.
    - `pairs` supports explicit pass-through comparisons as `(left, right)`
      tuples. If omitted, comparisons are inferred from all observed groups.
    """

    policy: DatasetComparisonBuildingPolicy = DATASET_COMPARISON_BUILDING_POLICY_NONE
    sample_group_column: str = DATASET_COMPARISON_BUILDING_DEFAULT_SAMPLE_GROUP_COLUMN
    pairs: tuple[DatasetComparisonPair, ...] | None = None

    def __post_init__(self) -> None:
        validate_comparison_building_config(
            policy=self.policy,
            sample_group_column=self.sample_group_column,
            pairs=self.pairs,
            supported_policies=DATASET_COMPARISON_BUILDING_POLICIES,
            policy_none=DATASET_COMPARISON_BUILDING_POLICY_NONE,
            policy_sample_metadata_pairs=(
                DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS
            ),
        )


__all__ = [
    "DATASET_COMPARISON_BUILDING_DEFAULT_SAMPLE_GROUP_COLUMN",
    "DATASET_COMPARISON_BUILDING_POLICIES",
    "DATASET_COMPARISON_BUILDING_POLICY_NONE",
    "DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS",
    "DatasetComparisonBuildingConfig",
    "DatasetComparisonBuildingPolicy",
    "DatasetComparisonPair",
]
