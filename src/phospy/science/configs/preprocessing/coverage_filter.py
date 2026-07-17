"""Group-aware phosphosite coverage filter configuration."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.science.configs.preprocessing._validation import (
    validate_group_coverage_filter_config,
)


@dataclass(frozen=True, slots=True)
class DatasetGroupCoverageFilterConfig:
    """Public config for condition/replicate-aware coverage filtering.

    This object describes the requested filtering rule. When enabled in
    `DatasetPreprocessingConfig`, the dataset preprocessing pipeline applies it
    before missing-data handling and analysis-ready dataset creation.

    When `enabled=True`, provide `group_column`, exactly one threshold, and the
    minimum number of groups that must pass the threshold:

    - `min_finite_observations_per_group`: minimum finite sample values in a
      group.
    - `min_finite_fraction_per_group`: minimum finite-value fraction in a group.
    - `min_groups_passing_threshold`: number of groups that must satisfy the
      selected threshold.
    """

    enabled: bool = False
    group_column: str | None = None
    min_finite_observations_per_group: int | None = None
    min_finite_fraction_per_group: float | None = None
    min_groups_passing_threshold: int = 1

    def __post_init__(self) -> None:
        validate_group_coverage_filter_config(
            enabled=self.enabled,
            group_column=self.group_column,
            min_finite_observations_per_group=(self.min_finite_observations_per_group),
            min_finite_fraction_per_group=self.min_finite_fraction_per_group,
            min_groups_passing_threshold=self.min_groups_passing_threshold,
        )


__all__ = [
    "DatasetGroupCoverageFilterConfig",
]
