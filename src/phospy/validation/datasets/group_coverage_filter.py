"""Compatibility facade for preprocessing-owned group coverage validation."""

from phospy.science.datasets.preprocessing.group_coverage_metadata import (
    GroupCoverageFilterMetadataValidator,
    ResolvedGroupCoverageFilterMetadata,
    require_numeric_group_coverage_matrix,
)

__all__ = [
    "GroupCoverageFilterMetadataValidator",
    "ResolvedGroupCoverageFilterMetadata",
    "require_numeric_group_coverage_matrix",
]
