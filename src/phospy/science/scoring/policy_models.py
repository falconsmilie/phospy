from __future__ import annotations

from phospy.policies import PolicyEnum


class DownstreamScoreSource(PolicyEnum):
    PROFILE_SCORES = "profile_scores"
    RANK_WEIGHTED_FUSION_SCORES = "rank_weighted_fusion_scores"
    KINASE_LIBRARY_MOTIF_SCORES = "kinase_library_motif_scores"
    COMBINED_PROFILE_MOTIF_SCORES = "combined_profile_motif_scores"


class ProfileSelfInclusionPolicy(PolicyEnum):
    ALLOW = "allow"
    LEAVE_ONE_OUT = "leave_one_out"


class ThresholdMode(PolicyEnum):
    GREATER_THAN = "score > threshold"
    GREATER_THAN_OR_EQUAL = "score >= threshold"


__all__ = [
    "DownstreamScoreSource",
    "ProfileSelfInclusionPolicy",
    "ThresholdMode",
]
