from __future__ import annotations

from phospy.policies import PolicyEnum


class DownstreamScoreSource(PolicyEnum):
    PROFILE_SCORES = "profile_scores"
    RANK_WEIGHTED_FUSION_SCORES = "rank_weighted_fusion_scores"


class ThresholdMode(PolicyEnum):
    GREATER_THAN = "score > threshold"
    GREATER_THAN_OR_EQUAL = "score >= threshold"


__all__ = ["DownstreamScoreSource", "ThresholdMode"]
