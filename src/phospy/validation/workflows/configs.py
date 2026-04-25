"""Shared workflow-config validation."""

from __future__ import annotations

from phospy.api.configs import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    SignalomeConfig,
)
from phospy.errors.validation import WorkflowValidationError


class KinaseWorkflowConfigValidator:
    """Validate kinase workflow config objects and local invariants."""

    def run(
        self,
        *,
        scoring_config: object,
        prediction_config: object,
        activity_config: object | None,
    ) -> tuple[
        KinaseScoringConfig,
        KinasePredictionConfig,
        KinaseActivityConfig | None,
    ]:
        validated_scoring = self._validate_scoring(scoring_config)
        validated_prediction = self._validate_prediction(prediction_config)
        validated_activity = self._validate_activity(activity_config)
        return validated_scoring, validated_prediction, validated_activity

    @staticmethod
    def _validate_scoring(config: object) -> KinaseScoringConfig:
        if not isinstance(config, KinaseScoringConfig):
            raise WorkflowValidationError(
                "kinase workflow request scoring_config must be KinaseScoringConfig"
            )
        return config

    @staticmethod
    def _validate_prediction(config: object) -> KinasePredictionConfig:
        if not isinstance(config, KinasePredictionConfig):
            raise WorkflowValidationError(
                "kinase workflow request prediction_config must be KinasePredictionConfig"
            )
        return config

    @staticmethod
    def _validate_activity(config: object | None) -> KinaseActivityConfig | None:
        if config is None:
            return None
        if not isinstance(config, KinaseActivityConfig):
            raise WorkflowValidationError(
                "kinase workflow request activity_config must be KinaseActivityConfig or None"
            )
        return config


class SignalomeConfigValidator:
    """Validate signalome workflow config objects and local invariants."""

    def run(self, config: object) -> SignalomeConfig:
        if not isinstance(config, SignalomeConfig):
            raise WorkflowValidationError(
                "signalome workflow request config must be SignalomeConfig"
            )
        return config


__all__ = ["KinaseWorkflowConfigValidator", "SignalomeConfigValidator"]
