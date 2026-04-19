"""Public workflow and stage configuration models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

KINASE_SCORING_MIN_SUBSTRATES_FLOOR = 2
KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT = "strict"
KINASE_PROFILE_MISSING_VALUE_STRATEGY_MEDIAN_SKIPNA = "median_skipna"
KinaseProfileMissingValueStrategy = Literal[
    "strict",
    "median_skipna",
]
KINASE_PROFILE_MISSING_VALUE_STRATEGIES = frozenset(
    {
        KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT,
        KINASE_PROFILE_MISSING_VALUE_STRATEGY_MEDIAN_SKIPNA,
    }
)
KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR = 1
KINASE_ACTIVITY_TOP_N_SUBSTRATES_FLOOR = 1
KINASE_ACTIVITY_DEFAULT_THRESHOLD = 0.6
KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES = 3
KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES = 20


@dataclass(frozen=True, slots=True)
class KinaseScoringConfig:
    """Public scoring-stage configuration.

    `min_substrates` is constrained to the scientific support floor so one-site
    kinase profiles are not part of the default public lane.

    `include_diagnostic_scoring_tables` controls publication of non-authoritative
    diagnostic scoring outputs (`motif_scores`, `weights`). The authoritative
    downstream lane (`combined_scores` with profile fallback) is always computed.

    `profile_missing_value_strategy` controls column-wise median behavior when a
    kinase profile is built from multiple quantified substrates:

    - `"strict"` propagates missing values (`median(..., skipna=False)`)
    - `"median_skipna"` ignores missing values (`median(..., skipna=True)`)
    """

    min_substrates: int = KINASE_SCORING_MIN_SUBSTRATES_FLOOR
    include_diagnostic_scoring_tables: bool = False
    profile_missing_value_strategy: KinaseProfileMissingValueStrategy = (
        KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT
    )


@dataclass(frozen=True, slots=True)
class KinasePredictionConfig:
    """Public prediction-stage configuration."""

    top_k: int = 30
    ensemble_size: int = 10


@dataclass(frozen=True, slots=True)
class KinaseActivityConfig:
    """Configuration for the supported kinase activity stage.

    Activity runs inside `KinaseWorkflow` and can be disabled by setting either:

    - `activity_config=None` on `KinaseWorkflowRequest`, or
    - `activity_config.enabled=False`.
    """

    enabled: bool = True
    threshold: float = KINASE_ACTIVITY_DEFAULT_THRESHOLD
    min_substrates: int = KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES
    top_n_substrates: int = KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES


@dataclass(frozen=True, slots=True)
class SignalomeConfig:
    """Public signalome workflow configuration."""

    substrate_support_cutoff: float = 0.5
    network_correlation_threshold: float = 0.5


__all__ = [
    "KINASE_PROFILE_MISSING_VALUE_STRATEGIES",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGY_MEDIAN_SKIPNA",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT",
    "KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES",
    "KINASE_ACTIVITY_DEFAULT_THRESHOLD",
    "KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES",
    "KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR",
    "KINASE_ACTIVITY_TOP_N_SUBSTRATES_FLOOR",
    "KINASE_SCORING_MIN_SUBSTRATES_FLOOR",
    "KinaseProfileMissingValueStrategy",
    "KinaseActivityConfig",
    "KinasePredictionConfig",
    "KinaseScoringConfig",
    "SignalomeConfig",
]
