"""Policy values persisted by Dataset preparation records."""

from __future__ import annotations

from phospy.policies import PolicyEnum


class MissingDataPolicy(PolicyEnum):
    FORBID = "forbid"
    IMPUTE_ROW_MEDIAN = "impute_row_median"
    IMPUTE_MINPROB = "impute_minprob"
    IMPUTE_KNN = "impute_knn"
    IMPUTE_GROUP_AWARE = "impute_group_aware"


class TotalProteinCorrectionPolicy(PolicyEnum):
    NONE = "none"
    SUBTRACT_LOG_TOTAL = "subtract_log_total"


class TotalProteinCorrectionIdentityMatchingPolicy(PolicyEnum):
    STRICT = "strict"
    GENE_SYMBOL_NORMALISED = "gene_symbol_normalised"


class NormalisationPolicy(PolicyEnum):
    NONE = "none"
    MEDIAN_CENTER = "median_center"
    QUANTILE = "quantile"


class SiteSequenceResolutionMode(PolicyEnum):
    VALIDATE_EXISTING_AND_FILL_MISSING = "validate_existing_and_fill_missing"
    FILL_MISSING_ONLY = "fill_missing_only"
    VALIDATE_EXISTING_ONLY = "validate_existing_only"
    REPLACE_EXISTING = "replace_existing"


class SiteSequenceConflictPolicy(PolicyEnum):
    ERROR = "error"
    PRESERVE_EXISTING = "preserve_existing"
    REPLACE_EXISTING = "replace_existing"


__all__ = [
    "MissingDataPolicy",
    "NormalisationPolicy",
    "SiteSequenceConflictPolicy",
    "SiteSequenceResolutionMode",
    "TotalProteinCorrectionIdentityMatchingPolicy",
    "TotalProteinCorrectionPolicy",
]
