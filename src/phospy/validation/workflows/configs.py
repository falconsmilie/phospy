"""Shared workflow-config validation."""

from __future__ import annotations

from phospy.api.configs import (
    KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR,
    KINASE_ACTIVITY_TOP_N_SUBSTRATES_FLOOR,
    KINASE_ADAPTIVE_POLICIES,
    KINASE_PREDICTION_MODES,
    KINASE_PROFILE_MISSING_VALUE_STRATEGIES,
    KINASE_SCORING_MIN_SUBSTRATES_FLOOR,
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
                "kinase workflow request scoring_config must be KinaseScoringConfig"
            )
        if not isinstance(config.include_diagnostic_scoring_tables, bool):
            raise WorkflowValidationError(
                "scoring_config.include_diagnostic_scoring_tables must be a bool"
            )
        if (
            config.profile_missing_value_strategy
            not in KINASE_PROFILE_MISSING_VALUE_STRATEGIES
        ):
            allowed = ", ".join(sorted(KINASE_PROFILE_MISSING_VALUE_STRATEGIES))
            raise WorkflowValidationError(
                "scoring_config.profile_missing_value_strategy must be one of: "
                f"{allowed}"
            )
        require_int_at_least(
            config.min_substrates,
            field_name="scoring_config.min_substrates",
            minimum=KINASE_SCORING_MIN_SUBSTRATES_FLOOR,
            error_type=WorkflowValidationError,
        )
        return config

    def run_kinase_prediction(self, config: object) -> KinasePredictionConfig:
        if not isinstance(config, KinasePredictionConfig):
            raise WorkflowValidationError(
                "kinase workflow request prediction_config must be KinasePredictionConfig"
            )
        if config.mode not in KINASE_PREDICTION_MODES:
            allowed_modes = ", ".join(sorted(KINASE_PREDICTION_MODES))
            raise WorkflowValidationError(
                f"prediction_config.mode must be one of: {allowed_modes}"
            )
        if config.adaptive_policy not in KINASE_ADAPTIVE_POLICIES:
            allowed_policies = ", ".join(sorted(KINASE_ADAPTIVE_POLICIES))
            raise WorkflowValidationError(
                f"prediction_config.adaptive_policy must be one of: {allowed_policies}"
            )
        if config.random_state is not None:
            require_int_at_least(
                config.random_state,
                field_name="prediction_config.random_state",
                minimum=0,
                error_type=WorkflowValidationError,
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
        require_int_at_least(
            config.n_iterations,
            field_name="prediction_config.n_iterations",
            minimum=1,
            error_type=WorkflowValidationError,
        )
        return config

    def run_kinase_activity(self, config: object | None) -> KinaseActivityConfig | None:
        if config is None:
            return None
        if not isinstance(config, KinaseActivityConfig):
            raise WorkflowValidationError(
                "kinase workflow request activity_config must be KinaseActivityConfig or None"
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
        require_int_at_least(
            config.min_substrates,
            field_name="activity_config.min_substrates",
            minimum=KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR,
            error_type=WorkflowValidationError,
        )
        require_int_at_least(
            config.top_n_substrates,
            field_name="activity_config.top_n_substrates",
            minimum=KINASE_ACTIVITY_TOP_N_SUBSTRATES_FLOOR,
            error_type=WorkflowValidationError,
        )
        return config

    def run_signalome(self, config: object) -> SignalomeConfig:
        if not isinstance(config, SignalomeConfig):
            raise WorkflowValidationError(
                "signalome workflow request config must be SignalomeConfig"
            )
        require_real_between(
            config.substrate_support_cutoff,
            field_name="signalome workflow request config.substrate_support_cutoff",
            minimum=0.0,
            maximum=1.0,
            error_type=WorkflowValidationError,
        )
        require_real_between(
            config.network_correlation_threshold,
            field_name=(
                "signalome workflow request config.network_correlation_threshold"
            ),
            minimum=0.0,
            maximum=1.0,
            error_type=WorkflowValidationError,
        )
        return config
