from __future__ import annotations

__phospy_contracts_facade_role__ = "science_owned_public_enum"

from phospy.policies import PolicyEnum


class TechnicalReplicatePolicy(PolicyEnum):
    """Policy for handling repeated biological replicate IDs."""

    REJECT = "reject"
    MEAN = "mean"
    MEDIAN = "median"


__all__ = ["TechnicalReplicatePolicy"]
