"""Public kinase prediction configuration models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from phospy.api.configs.common import _require_int_at_least
from phospy.errors.validation import WorkflowValidationError

KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING = "deterministic_ranking"
KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE = "adaptive_ensemble"
KinasePredictionMode = Literal[
    "deterministic_ranking",
    "adaptive_ensemble",
]
KINASE_PREDICTION_MODES = frozenset(
    {
        KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
        KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE,
    }
)
KINASE_ADAPTIVE_POLICY_STABLE = "stable"
KINASE_ADAPTIVE_POLICY_R_PARITY = "r_parity"
KinaseAdaptivePolicy = Literal["stable", "r_parity"]
KINASE_ADAPTIVE_POLICIES = frozenset(
    {
        KINASE_ADAPTIVE_POLICY_STABLE,
        KINASE_ADAPTIVE_POLICY_R_PARITY,
    }
)
KINASE_PREDICTION_DEFAULT_ITERATIONS = 5
KINASE_PREDICTION_DEFAULT_DETERMINISTIC_MAX_SELECTED_KINASES = 10
KINASE_PREDICTION_DEFAULT_ADAPTIVE_ENSEMBLE_RUNS = 10


@dataclass(frozen=True, slots=True)
class KinasePredictionConfig:
    """Public prediction-stage configuration.

    `mode` selects the prediction lane:

    - `"deterministic_ranking"`: deterministic top-kinase selection from
      downstream scores.
    - `"adaptive_ensemble"`: real adaptive ensemble execution ported from donor
      science.

    `deterministic_max_selected_kinases` controls deterministic lane breadth.
    `adaptive_ensemble_runs` controls adaptive lane ensemble executions.
    """

    top_k: int = 30
    deterministic_max_selected_kinases: int = (
        KINASE_PREDICTION_DEFAULT_DETERMINISTIC_MAX_SELECTED_KINASES
    )
    adaptive_ensemble_runs: int = KINASE_PREDICTION_DEFAULT_ADAPTIVE_ENSEMBLE_RUNS
    mode: KinasePredictionMode = KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING
    adaptive_policy: KinaseAdaptivePolicy = KINASE_ADAPTIVE_POLICY_STABLE
    n_iterations: int = KINASE_PREDICTION_DEFAULT_ITERATIONS
    random_state: int | None = None

    def __post_init__(self) -> None:
        if self.mode not in KINASE_PREDICTION_MODES:
            allowed_modes = ", ".join(sorted(KINASE_PREDICTION_MODES))
            raise WorkflowValidationError(
                f"prediction_config.mode must be one of: {allowed_modes}"
            )
        if (
            self.mode == KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE
            and self.random_state is None
        ):
            raise WorkflowValidationError(
                "prediction_config.random_state must be provided when "
                "prediction_config.mode='adaptive_ensemble' so adaptive prediction "
                "runs are reproducible"
            )
        if self.adaptive_policy not in KINASE_ADAPTIVE_POLICIES:
            allowed_policies = ", ".join(sorted(KINASE_ADAPTIVE_POLICIES))
            raise WorkflowValidationError(
                f"prediction_config.adaptive_policy must be one of: {allowed_policies}"
            )
        if self.random_state is not None:
            _require_int_at_least(
                self.random_state,
                field_name="prediction_config.random_state",
                minimum=0,
                error_type=WorkflowValidationError,
            )
        _require_int_at_least(
            self.top_k,
            field_name="prediction_config.top_k",
            minimum=1,
            error_type=WorkflowValidationError,
        )
        _require_int_at_least(
            self.deterministic_max_selected_kinases,
            field_name="prediction_config.deterministic_max_selected_kinases",
            minimum=1,
            error_type=WorkflowValidationError,
        )
        _require_int_at_least(
            self.adaptive_ensemble_runs,
            field_name="prediction_config.adaptive_ensemble_runs",
            minimum=1,
            error_type=WorkflowValidationError,
        )
        _require_int_at_least(
            self.n_iterations,
            field_name="prediction_config.n_iterations",
            minimum=1,
            error_type=WorkflowValidationError,
        )


__all__ = [
    "KINASE_ADAPTIVE_POLICIES",
    "KINASE_ADAPTIVE_POLICY_R_PARITY",
    "KINASE_ADAPTIVE_POLICY_STABLE",
    "KINASE_PREDICTION_DEFAULT_ITERATIONS",
    "KINASE_PREDICTION_DEFAULT_ADAPTIVE_ENSEMBLE_RUNS",
    "KINASE_PREDICTION_DEFAULT_DETERMINISTIC_MAX_SELECTED_KINASES",
    "KINASE_PREDICTION_MODES",
    "KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE",
    "KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING",
    "KinaseAdaptivePolicy",
    "KinasePredictionConfig",
    "KinasePredictionMode",
]
