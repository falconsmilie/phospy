"""Public workflow and stage configuration models."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class KinaseScoringConfig:
    """Public scoring-stage configuration."""

    min_substrates: int = 1


@dataclass(frozen=True, slots=True)
class KinasePredictionConfig:
    """Public prediction-stage configuration."""

    top_k: int = 30
    ensemble_size: int = 10


@dataclass(frozen=True, slots=True)
class KinaseActivityConfig:
    """Public kinase activity-stage configuration."""

    enabled: bool = True
    threshold: float = 0.6


@dataclass(frozen=True, slots=True)
class SignalomeConfig:
    """Public signalome workflow configuration."""

    substrate_support_cutoff: float = 0.5
    network_correlation_threshold: float = 0.5


__all__ = [
    "KinaseActivityConfig",
    "KinasePredictionConfig",
    "KinaseScoringConfig",
    "SignalomeConfig",
]
