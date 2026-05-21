"""Intensity-transform preprocessing policy configuration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from phospy.validation.configs.preprocessing import (
    validate_intensity_transform_config,
)

DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY = "identity"
DATASET_INTENSITY_TRANSFORM_POLICY_LOG2 = "log2"
DatasetIntensityTransformPolicy = Literal["identity", "log2"]
DATASET_INTENSITY_TRANSFORM_POLICIES = frozenset(
    {
        DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY,
        DATASET_INTENSITY_TRANSFORM_POLICY_LOG2,
    }
)


@dataclass(frozen=True, slots=True)
class DatasetIntensityTransformConfig:
    """Public intensity transform policy options for dataset building.

    - `"identity"`: no transform (strict default).
    - `"log2"`: apply `log2(value + pseudocount)` to quantitative matrix values.

    `pseudocount` must be non-negative.
    """

    policy: DatasetIntensityTransformPolicy = (
        DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY
    )
    pseudocount: float = 1.0

    def __post_init__(self) -> None:
        validate_intensity_transform_config(
            policy=self.policy,
            pseudocount=self.pseudocount,
            supported_policies=DATASET_INTENSITY_TRANSFORM_POLICIES,
        )


__all__ = [
    "DATASET_INTENSITY_TRANSFORM_POLICIES",
    "DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY",
    "DATASET_INTENSITY_TRANSFORM_POLICY_LOG2",
    "DatasetIntensityTransformConfig",
    "DatasetIntensityTransformPolicy",
]
