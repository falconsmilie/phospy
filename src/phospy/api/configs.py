"""Public workflow and stage configuration models."""

from __future__ import annotations

from dataclasses import dataclass

KINASE_SCORING_MIN_SUBSTRATES_FLOOR = 2
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
    """

    min_substrates: int = KINASE_SCORING_MIN_SUBSTRATES_FLOOR


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
    "KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES",
    "KINASE_ACTIVITY_DEFAULT_THRESHOLD",
    "KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES",
    "KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR",
    "KINASE_ACTIVITY_TOP_N_SUBSTRATES_FLOOR",
    "KINASE_SCORING_MIN_SUBSTRATES_FLOOR",
    "KinaseActivityConfig",
    "KinasePredictionConfig",
    "KinaseScoringConfig",
    "SignalomeConfig",
]
