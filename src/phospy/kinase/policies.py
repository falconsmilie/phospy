"""Capability-neutral Kinase policy definitions."""

from __future__ import annotations

from phospy.core.policies import PolicyEnum


class ThresholdMode(PolicyEnum):
    GREATER_THAN = "score > threshold"
    GREATER_THAN_OR_EQUAL = "score >= threshold"


__all__ = ["ThresholdMode"]
