from __future__ import annotations

from phospy.policies import PolicyEnum


class ImputationInputScale(PolicyEnum):
    LINEAR = "linear"
    LOG2 = "log2"


class IntensityTransformPolicy(PolicyEnum):
    IDENTITY = "identity"
    LOG2 = "log2"


class SiteMatrixPolicy(PolicyEnum):
    AS_INPUT = "as_input"
    BUILD_FROM_METADATA = "build_from_metadata"


class SiteMatrixDuplicateSitePolicy(PolicyEnum):
    MAX_MEAN_SIGNAL = "max_mean_signal"
    FIRST = "first"
    AGGREGATE_MEAN = "aggregate_mean"
    AGGREGATE_MEDIAN = "aggregate_median"
    ERROR = "error"


class SiteMatrixMissingDataPolicy(PolicyEnum):
    DROP_ANY_MISSING = "drop_any_missing"
    RETAIN_MISSING = "retain_missing"
    REQUIRE_MIN_OBSERVED_VALUES = "require_min_observed_values"


class ComparisonBuildingPolicy(PolicyEnum):
    NONE = "none"
    SAMPLE_METADATA_PAIRS = "sample_metadata_pairs"


class LocalisationEligibilityMode(PolicyEnum):
    REQUIRE_THRESHOLD = "require_threshold"
    ALLOW_MISSING_WITH_WAIVER = "allow_missing_with_waiver"
    IGNORE = "ignore"


__all__ = [
    "ComparisonBuildingPolicy",
    "ImputationInputScale",
    "IntensityTransformPolicy",
    "LocalisationEligibilityMode",
    "SiteMatrixDuplicateSitePolicy",
    "SiteMatrixMissingDataPolicy",
    "SiteMatrixPolicy",
]
