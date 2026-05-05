"""Internal enum-backed policy models for high-impact scientific behaviors."""

from __future__ import annotations

from enum import Enum
from typing import TypeVar

from phospy.errors.input import PhosPyInputError
from phospy.validation.common.config_values import coerce_policy_enum

_PolicyEnumT = TypeVar("_PolicyEnumT", bound="_PolicyEnum")


class _PolicyEnum(str, Enum):
    """Base class for stable policy enums with strict parsing helpers."""

    def __str__(self) -> str:
        return self.value

    @classmethod
    def parse(
        cls: type[_PolicyEnumT],
        value: object,
        *,
        field_name: str,
    ) -> _PolicyEnumT:
        return coerce_policy_enum(
            cls,
            value,
            field_name=field_name,
            error_type=PhosPyInputError,
        )


class MissingDataPolicy(_PolicyEnum):
    FORBID = "forbid"
    IMPUTE_ROW_MEDIAN = "impute_row_median"
    IMPUTE_MINPROB = "impute_minprob"
    IMPUTE_KNN = "impute_knn"


class TotalProteinCorrectionPolicy(_PolicyEnum):
    NONE = "none"
    SUBTRACT_LOG_TOTAL = "subtract_log_total"


class TotalProteinCorrectionIdentityMatchingPolicy(_PolicyEnum):
    STRICT = "strict"
    GENE_SYMBOL_NORMALISED = "gene_symbol_normalised"


class IntensityTransformPolicy(_PolicyEnum):
    IDENTITY = "identity"
    LOG2 = "log2"


class NormalisationPolicy(_PolicyEnum):
    NONE = "none"
    MEDIAN_CENTER = "median_center"
    QUANTILE = "quantile"


class SiteMatrixPolicy(_PolicyEnum):
    AS_INPUT = "as_input"
    BUILD_FROM_METADATA = "build_from_metadata"


class SiteMatrixDuplicateSitePolicy(_PolicyEnum):
    MAX_MEAN_SIGNAL = "max_mean_signal"
    FIRST = "first"
    AGGREGATE_MEAN = "aggregate_mean"
    AGGREGATE_MEDIAN = "aggregate_median"
    ERROR = "error"


class SiteMatrixMissingDataPolicy(_PolicyEnum):
    DROP_ANY_MISSING = "drop_any_missing"
    RETAIN_MISSING = "retain_missing"
    REQUIRE_MIN_OBSERVED_VALUES = "require_min_observed_values"


class SiteSequenceResolutionMode(_PolicyEnum):
    VALIDATE_EXISTING_AND_FILL_MISSING = "validate_existing_and_fill_missing"
    FILL_MISSING_ONLY = "fill_missing_only"
    VALIDATE_EXISTING_ONLY = "validate_existing_only"
    REPLACE_EXISTING = "replace_existing"


class SiteSequenceConflictPolicy(_PolicyEnum):
    ERROR = "error"
    PRESERVE_EXISTING = "preserve_existing"
    REPLACE_EXISTING = "replace_existing"


class ComparisonBuildingPolicy(_PolicyEnum):
    NONE = "none"
    SAMPLE_METADATA_PAIRS = "sample_metadata_pairs"


class DownstreamScoreSource(_PolicyEnum):
    PROFILE_SCORES = "profile_scores"
    RANK_WEIGHTED_FUSION_SCORES = "rank_weighted_fusion_scores"


class ThresholdMode(_PolicyEnum):
    GREATER_THAN = "score > threshold"
    GREATER_THAN_OR_EQUAL = "score >= threshold"


__all__ = [
    "ComparisonBuildingPolicy",
    "DownstreamScoreSource",
    "IntensityTransformPolicy",
    "MissingDataPolicy",
    "NormalisationPolicy",
    "SiteMatrixDuplicateSitePolicy",
    "SiteMatrixMissingDataPolicy",
    "SiteMatrixPolicy",
    "SiteSequenceConflictPolicy",
    "SiteSequenceResolutionMode",
    "ThresholdMode",
    "TotalProteinCorrectionIdentityMatchingPolicy",
    "TotalProteinCorrectionPolicy",
]
