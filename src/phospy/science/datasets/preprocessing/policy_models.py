from __future__ import annotations

from phospy.policies import PolicyEnum


class MissingDataPolicy(PolicyEnum):
    FORBID = "forbid"
    IMPUTE_ROW_MEDIAN = "impute_row_median"
    IMPUTE_MINPROB = "impute_minprob"
    IMPUTE_KNN = "impute_knn"


class ImputationInputScale(PolicyEnum):
    LINEAR = "linear"
    LOG2 = "log2"


class TotalProteinCorrectionPolicy(PolicyEnum):
    NONE = "none"
    SUBTRACT_LOG_TOTAL = "subtract_log_total"


class TotalProteinCorrectionIdentityMatchingPolicy(PolicyEnum):
    STRICT = "strict"
    GENE_SYMBOL_NORMALISED = "gene_symbol_normalised"


class IntensityTransformPolicy(PolicyEnum):
    IDENTITY = "identity"
    LOG2 = "log2"


class NormalisationPolicy(PolicyEnum):
    NONE = "none"
    MEDIAN_CENTER = "median_center"
    QUANTILE = "quantile"


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


class SiteSequenceResolutionMode(PolicyEnum):
    VALIDATE_EXISTING_AND_FILL_MISSING = "validate_existing_and_fill_missing"
    FILL_MISSING_ONLY = "fill_missing_only"
    VALIDATE_EXISTING_ONLY = "validate_existing_only"
    REPLACE_EXISTING = "replace_existing"


class SiteSequenceConflictPolicy(PolicyEnum):
    ERROR = "error"
    PRESERVE_EXISTING = "preserve_existing"
    REPLACE_EXISTING = "replace_existing"


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
    "MissingDataPolicy",
    "NormalisationPolicy",
    "SiteMatrixDuplicateSitePolicy",
    "SiteMatrixMissingDataPolicy",
    "SiteMatrixPolicy",
    "SiteSequenceConflictPolicy",
    "SiteSequenceResolutionMode",
    "TotalProteinCorrectionIdentityMatchingPolicy",
    "TotalProteinCorrectionPolicy",
]
