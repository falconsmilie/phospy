"""Localisation-confidence preprocessing policy configuration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from phospy.contracts.configs.preprocessing._validation import (
    validate_localisation_config,
)

DATASET_LOCALISATION_MODE_REQUIRE_THRESHOLD = "require_threshold"
DATASET_LOCALISATION_MODE_ALLOW_MISSING_WITH_WAIVER = "allow_missing_with_waiver"
DATASET_LOCALISATION_MODE_IGNORE = "ignore"
DatasetLocalisationMode = Literal[
    "require_threshold",
    "allow_missing_with_waiver",
    "ignore",
]
DATASET_LOCALISATION_MODES = frozenset(
    {
        DATASET_LOCALISATION_MODE_REQUIRE_THRESHOLD,
        DATASET_LOCALISATION_MODE_ALLOW_MISSING_WITH_WAIVER,
        DATASET_LOCALISATION_MODE_IGNORE,
    }
)


@dataclass(frozen=True, slots=True)
class DatasetLocalisationConfig:
    """Public localisation-confidence eligibility policy for dataset building."""

    mode: DatasetLocalisationMode = DATASET_LOCALISATION_MODE_REQUIRE_THRESHOLD
    min_confidence: float = 0.75
    confidence_column: str = "localisation_confidence"
    waiver_reason: str | None = None

    def __post_init__(self) -> None:
        validate_localisation_config(
            mode=self.mode,
            min_confidence=self.min_confidence,
            confidence_column=self.confidence_column,
            waiver_reason=self.waiver_reason,
            supported_modes=DATASET_LOCALISATION_MODES,
            mode_allow_missing_with_waiver=DATASET_LOCALISATION_MODE_ALLOW_MISSING_WITH_WAIVER,
        )


__all__ = [
    "DATASET_LOCALISATION_MODES",
    "DATASET_LOCALISATION_MODE_ALLOW_MISSING_WITH_WAIVER",
    "DATASET_LOCALISATION_MODE_IGNORE",
    "DATASET_LOCALISATION_MODE_REQUIRE_THRESHOLD",
    "DatasetLocalisationConfig",
    "DatasetLocalisationMode",
]
