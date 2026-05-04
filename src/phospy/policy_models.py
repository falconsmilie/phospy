"""Internal enum-backed policy models for high-impact scientific behaviors."""

from __future__ import annotations

from enum import Enum
from typing import TypeVar

from phospy.errors.input import PhosPyInputError

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
        if isinstance(value, cls):
            return value
        if isinstance(value, str):
            normalized = value.strip()
            try:
                return cls(normalized)
            except ValueError:
                pass
        supported = ", ".join(member.value for member in cls)
        raise PhosPyInputError(
            f"{field_name} must be one of: {supported}; got {value!r}"
        )


class MissingDataPolicy(_PolicyEnum):
    FORBID = "forbid"
    IMPUTE_ROW_MEDIAN = "impute_row_median"
    IMPUTE_MINPROB = "impute_minprob"
    IMPUTE_KNN = "impute_knn"


class TotalProteinCorrectionPolicy(_PolicyEnum):
    NONE = "none"
    SUBTRACT_LOG_TOTAL = "subtract_log_total"


class DownstreamScoreSource(_PolicyEnum):
    PROFILE_SCORES = "profile_scores"
    RANK_WEIGHTED_FUSION_SCORES = "rank_weighted_fusion_scores"


class ThresholdMode(_PolicyEnum):
    GREATER_THAN = "score > threshold"
    GREATER_THAN_OR_EQUAL = "score >= threshold"


__all__ = [
    "DownstreamScoreSource",
    "MissingDataPolicy",
    "ThresholdMode",
    "TotalProteinCorrectionPolicy",
]
