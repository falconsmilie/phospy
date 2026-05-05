"""Comparison-building preprocessing policy configuration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from phospy.errors.input import PhosPyInputError

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
        policy = self.policy
        if policy not in DATASET_COMPARISON_BUILDING_POLICIES:
            supported = ", ".join(sorted(DATASET_COMPARISON_BUILDING_POLICIES))
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons.policy "
                f"must be one of: {supported}"
            )
        sample_group_column = self.sample_group_column
        if not isinstance(sample_group_column, str) or not sample_group_column.strip():
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons."
                "sample_group_column must be a non-empty string"
            )
        pairs = self.pairs
        if policy == DATASET_COMPARISON_BUILDING_POLICY_NONE:
            if pairs is not None:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.comparisons.pairs "
                    "must be None when comparisons.policy='none'"
                )
            return
        if policy != DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS:
            raise PhosPyInputError(
                "dataset build request preprocessing_config contains an unsupported "
                "comparisons.policy"
            )
        if pairs is None:
            return
        if not isinstance(pairs, (tuple, list)):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons.pairs must be "
                "a sequence of (left_group, right_group) pairs when provided"
            )
        resolved_pairs = tuple(pairs)
        if not resolved_pairs:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons.pairs must "
                "contain at least one pair when provided"
            )
        seen_pairs: set[tuple[str, str]] = set()
        for pair in resolved_pairs:
            if not isinstance(pair, tuple) or len(pair) != 2:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.comparisons.pairs "
                    "must contain only (left_group, right_group) tuples"
                )
            left_group, right_group = pair
            if not isinstance(left_group, str) or not left_group.strip():
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.comparisons.pairs "
                    "must contain non-empty left_group strings"
                )
            if not isinstance(right_group, str) or not right_group.strip():
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.comparisons.pairs "
                    "must contain non-empty right_group strings"
                )
            left = left_group.strip()
            right = right_group.strip()
            if left == right:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.comparisons.pairs "
                    "cannot contain self-comparison pairs"
                )
            canonical_pair = (left, right) if left <= right else (right, left)
            if canonical_pair in seen_pairs:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.comparisons.pairs "
                    "contains duplicate pairs regardless of direction"
                )
            seen_pairs.add(canonical_pair)


__all__ = [
    "DATASET_COMPARISON_BUILDING_DEFAULT_SAMPLE_GROUP_COLUMN",
    "DATASET_COMPARISON_BUILDING_POLICIES",
    "DATASET_COMPARISON_BUILDING_POLICY_NONE",
    "DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS",
    "DatasetComparisonBuildingConfig",
    "DatasetComparisonBuildingPolicy",
    "DatasetComparisonPair",
]
