"""Reliability-profile resolution for differential workflow validation."""

from __future__ import annotations

from typing import cast

from phospy.contracts.configs import DifferentialAnalysisConfig
from phospy.contracts.configs.differential import (
    DIFFERENTIAL_EXPLORATORY_MINIMUM_CONDITION_REPLICATES,
    DIFFERENTIAL_PRODUCTION_MINIMUM_CONDITION_REPLICATES,
    DIFFERENTIAL_RELIABILITY_PROFILE_EXPLORATORY_SINGLE_REPLICATE,
    DIFFERENTIAL_RELIABILITY_PROFILE_PRODUCTION,
    SUPPORTED_DIFFERENTIAL_RELIABILITY_PROFILES,
    DifferentialReliabilityProfile,
)
from phospy.errors.validation import WorkflowValidationError


def validate_differential_reliability_config(
    config: DifferentialAnalysisConfig,
) -> None:
    """Validate public reliability-profile semantics before design validation."""

    profile = config.reliability_profile
    if profile not in SUPPORTED_DIFFERENTIAL_RELIABILITY_PROFILES:
        supported = ", ".join(
            repr(value) for value in SUPPORTED_DIFFERENTIAL_RELIABILITY_PROFILES
        )
        raise WorkflowValidationError(
            "differential workflow request reliability_profile must be one of: "
            f"{supported}"
        )
    if not isinstance(
        cast(object, config.minimum_condition_replicates), int
    ) or isinstance(cast(object, config.minimum_condition_replicates), bool):
        raise WorkflowValidationError(
            "differential workflow request minimum_condition_replicates must be an int"
        )
    if config.minimum_condition_replicates < 1:
        raise WorkflowValidationError(
            "differential workflow request minimum_condition_replicates must be >= 1"
        )
    if (
        profile == DIFFERENTIAL_RELIABILITY_PROFILE_PRODUCTION
        and config.minimum_condition_replicates
        < DIFFERENTIAL_PRODUCTION_MINIMUM_CONDITION_REPLICATES
    ):
        raise WorkflowValidationError(
            "differential reliability_profile='production' requires at least two "
            "biological replicates per contrasted condition. "
            "minimum_condition_replicates=1 is allowed only through the explicit "
            "reliability_profile='exploratory_single_replicate' opt-in."
        )


def resolved_minimum_condition_replicates(
    config: DifferentialAnalysisConfig,
) -> int:
    """Return the biological-replicate threshold implied by the profile."""

    if (
        config.reliability_profile
        == DIFFERENTIAL_RELIABILITY_PROFILE_EXPLORATORY_SINGLE_REPLICATE
    ):
        return DIFFERENTIAL_EXPLORATORY_MINIMUM_CONDITION_REPLICATES
    return max(
        int(config.minimum_condition_replicates),
        DIFFERENTIAL_PRODUCTION_MINIMUM_CONDITION_REPLICATES,
    )


def resolved_reliability_profile(
    config: DifferentialAnalysisConfig,
) -> DifferentialReliabilityProfile:
    """Return the normalized reliability profile after contract validation."""

    return cast(DifferentialReliabilityProfile, str(config.reliability_profile))


__all__ = [
    "resolved_minimum_condition_replicates",
    "resolved_reliability_profile",
    "validate_differential_reliability_config",
]
