from __future__ import annotations

from phospy.policies import PolicyEnum


class TechnicalReplicatePolicy(PolicyEnum):
    """Policy for handling repeated biological replicate IDs."""

    REJECT = "reject"
    MEAN = "mean"
    MEDIAN = "median"


__all__ = ["TechnicalReplicatePolicy"]
