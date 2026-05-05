"""Intensity-transform preprocessing policy configuration."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

from phospy.errors.input import PhosPyInputError

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
        policy = self.policy
        if policy not in DATASET_INTENSITY_TRANSFORM_POLICIES:
            supported = ", ".join(sorted(DATASET_INTENSITY_TRANSFORM_POLICIES))
            raise PhosPyInputError(
                "dataset build request preprocessing_config.intensity_transform."
                f"policy must be one of: {supported}"
            )

        pseudocount = self.pseudocount
        if isinstance(pseudocount, bool) or not isinstance(pseudocount, (int, float)):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.intensity_transform."
                "pseudocount must be a float or int"
            )
        if not math.isfinite(float(pseudocount)):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.intensity_transform."
                "pseudocount must be finite"
            )
        if pseudocount < 0:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.intensity_transform."
                "pseudocount must be greater than or equal to 0"
            )


__all__ = [
    "DATASET_INTENSITY_TRANSFORM_POLICIES",
    "DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY",
    "DATASET_INTENSITY_TRANSFORM_POLICY_LOG2",
    "DatasetIntensityTransformConfig",
    "DatasetIntensityTransformPolicy",
]
