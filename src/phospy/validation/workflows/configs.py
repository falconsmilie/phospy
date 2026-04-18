"""Shared workflow-config validation."""

from __future__ import annotations

from phospy.api.configs import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    SignalomeConfig,
)
from phospy.errors.validation import WorkflowValidationError
from phospy.validation.common.numbers import require_int_at_least, require_real_between


class WorkflowConfigValidator:
    """Validate workflow config objects and local invariants."""

    def run_kinase_scoring(self, config: object) -> KinaseScoringConfig:
        if not isinstance(config, KinaseScoringConfig):
            raise WorkflowValidationError(
                "simple kinase workflow request scoring_config must be KinaseScoringConfig"
            )
        require_int_at_least(
            config.min_substrates,
            field_name="scoring_config.min_substrates",
            minimum=1,
            error_type=WorkflowValidationError,
        )
        return config

    def run_kinase_prediction(self, config: object) -> KinasePredictionConfig:
        if not isinstance(config, KinasePredictionConfig):
            raise WorkflowValidationError(
                "simple kinase workflow request prediction_config must be KinasePredictionConfig"
            )
        require_int_at_least(
            config.top_k,
            field_name="prediction_config.top_k",
            minimum=1,
            error_type=WorkflowValidationError,
        )
        require_int_at_least(
            config.ensemble_size,
            field_name="prediction_config.ensemble_size",
            minimum=1,
            error_type=WorkflowValidationError,
        )
        return config

    def run_kinase_activity(self, config: object | None) -> KinaseActivityConfig | None:
        if config is None:
            return None
        if not isinstance(config, KinaseActivityConfig):
            raise WorkflowValidationError(
                "simple kinase workflow request activity_config must be KinaseActivityConfig or None"
            )
        if not isinstance(config.enabled, bool):
            raise WorkflowValidationError("activity_config.enabled must be a bool")
        require_real_between(
            config.threshold,
            field_name="activity_config.threshold",
            minimum=0.0,
            maximum=1.0,
            error_type=WorkflowValidationError,
        )
        return config

    def run_signalome(self, config: object) -> SignalomeConfig:
        if not isinstance(config, SignalomeConfig):
            raise WorkflowValidationError(
                "signalome workflow request config must be SignalomeConfig"
            )
        require_real_between(
            config.signalome_cutoff,
            field_name="signalome workflow request config.signalome_cutoff",
            minimum=0.0,
            maximum=1.0,
            error_type=WorkflowValidationError,
        )
        return config
