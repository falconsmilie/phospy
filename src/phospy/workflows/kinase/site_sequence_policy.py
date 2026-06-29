"""Kinase site-sequence policy parsing helpers."""

from __future__ import annotations

from phospy.contracts.configs import KinaseSiteSequenceConflictPolicy
from phospy.validation.common.config_values import coerce_policy_enum


def resolve_site_sequence_conflict_policy(
    value: object,
    *,
    field_name: str,
    error_type: type[Exception],
) -> KinaseSiteSequenceConflictPolicy:
    return coerce_policy_enum(
        KinaseSiteSequenceConflictPolicy,
        value,
        field_name=field_name,
        error_type=error_type,
    )


__all__ = ["resolve_site_sequence_conflict_policy"]
