"""Public workflow and stage configuration models."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.errors.validation import WorkflowValidationError


@dataclass(frozen=True, slots=True)
class KinaseScoringConfig:
    """Public scoring-stage configuration."""

    min_substrates: int = 1

    def __post_init__(self) -> None:
        if isinstance(self.min_substrates, bool) or not isinstance(
            self.min_substrates, int
        ):
            raise WorkflowValidationError(
                "scoring_config.min_substrates must be an int"
            )
        if self.min_substrates < 1:
            raise WorkflowValidationError(
                "scoring_config.min_substrates must be greater than or equal to 1"
            )


@dataclass(frozen=True, slots=True)
class KinasePredictionConfig:
    """Public prediction-stage configuration."""

    top_k: int = 30
    ensemble_size: int = 10

    def __post_init__(self) -> None:
        if isinstance(self.top_k, bool) or not isinstance(self.top_k, int):
            raise WorkflowValidationError("prediction_config.top_k must be an int")
        if self.top_k < 1:
            raise WorkflowValidationError(
                "prediction_config.top_k must be greater than or equal to 1"
            )
        if isinstance(self.ensemble_size, bool) or not isinstance(
            self.ensemble_size, int
        ):
            raise WorkflowValidationError(
                "prediction_config.ensemble_size must be an int"
            )
        if self.ensemble_size < 1:
            raise WorkflowValidationError(
                "prediction_config.ensemble_size must be greater than or equal to 1"
            )


@dataclass(frozen=True, slots=True)
class KinaseActivityConfig:
    """Public kinase activity-stage configuration."""

    enabled: bool = True
    threshold: float = 0.6

    def __post_init__(self) -> None:
        if not isinstance(self.enabled, bool):
            raise WorkflowValidationError("activity_config.enabled must be a bool")
        if isinstance(self.threshold, bool) or not isinstance(
            self.threshold, (int, float)
        ):
            raise WorkflowValidationError(
                "activity_config.threshold must be a float between 0.0 and 1.0"
            )
        if not 0.0 <= float(self.threshold) <= 1.0:
            raise WorkflowValidationError(
                "activity_config.threshold must be between 0.0 and 1.0"
            )


@dataclass(frozen=True, slots=True)
class SignalomeConfig:
    """Public signalome workflow configuration."""

    signalome_cutoff: float = 0.5


__all__ = [
    "KinaseActivityConfig",
    "KinasePredictionConfig",
    "KinaseScoringConfig",
    "SignalomeConfig",
]
