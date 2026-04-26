"""Config snapshot contract for kinase bundles."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING

from phospy.api.configs import (
    KINASE_ADAPTIVE_POLICY_STABLE,
    KINASE_PREDICTION_DEFAULT_ADAPTIVE_ENSEMBLE_RUNS,
    KINASE_PREDICTION_DEFAULT_DETERMINISTIC_MAX_SELECTED_KINASES,
    KINASE_PREDICTION_DEFAULT_ITERATIONS,
    KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
    KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
)
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.primitives import (
    require_bool,
    require_float,
    require_int,
    require_mapping,
    require_str,
)

if TYPE_CHECKING:
    from phospy.api.requests import KinaseWorkflowRequest


@dataclass(frozen=True, slots=True)
class KinaseWorkflowConfigSnapshot:
    """Serializable snapshot of the kinase workflow configuration."""

    scoring_config: KinaseScoringConfig
    prediction_config: KinasePredictionConfig
    activity_config: KinaseActivityConfig | None

    @classmethod
    def from_request(
        cls, request: KinaseWorkflowRequest
    ) -> KinaseWorkflowConfigSnapshot:
        """Create a config snapshot from a workflow request."""

        from phospy.api.requests import KinaseWorkflowRequest

        if not isinstance(request, KinaseWorkflowRequest):
            raise PhosPyInputError(
                "config snapshot request must be a KinaseWorkflowRequest"
            )
        return cls(
            scoring_config=request.scoring_config,
            prediction_config=request.prediction_config,
            activity_config=request.activity_config,
        )

    def to_payload(self) -> dict[str, object]:
        """Return a manifest-safe JSON payload for this config snapshot."""

        activity_payload: dict[str, object] | None
        if self.activity_config is None:
            activity_payload = None
        else:
            activity_payload = {
                "enabled": self.activity_config.enabled,
                "threshold": float(self.activity_config.threshold),
                "min_substrates": self.activity_config.min_substrates,
                "top_n_substrates": self.activity_config.top_n_substrates,
            }
        return {
            "scoring_config": {
                "min_substrates": self.scoring_config.min_substrates,
                "include_diagnostic_scoring_tables": (
                    self.scoring_config.include_diagnostic_scoring_tables
                ),
                "profile_missing_value_strategy": (
                    self.scoring_config.profile_missing_value_strategy
                ),
            },
            "prediction_config": {
                "top_k": self.prediction_config.top_k,
                "deterministic_max_selected_kinases": (
                    self.prediction_config.deterministic_max_selected_kinases
                ),
                "adaptive_ensemble_runs": self.prediction_config.adaptive_ensemble_runs,
                "mode": self.prediction_config.mode,
                "adaptive_policy": self.prediction_config.adaptive_policy,
                "n_iterations": self.prediction_config.n_iterations,
                "random_state": self.prediction_config.random_state,
            },
            "activity_config": activity_payload,
        }

    @classmethod
    def from_payload(
        cls, payload: Mapping[str, object]
    ) -> KinaseWorkflowConfigSnapshot:
        """Create a config snapshot from a decoded JSON payload."""

        scope = "config snapshot"
        scoring_payload = require_mapping(
            payload.get("scoring_config"),
            field_name=f"{scope}.scoring_config",
        )
        prediction_payload = require_mapping(
            payload.get("prediction_config"),
            field_name=f"{scope}.prediction_config",
        )
        activity_raw = payload.get("activity_config")
        if activity_raw is None:
            activity_config = None
        else:
            activity_payload = require_mapping(
                activity_raw,
                field_name=f"{scope}.activity_config",
            )
            activity_config = KinaseActivityConfig(
                enabled=require_bool(
                    activity_payload.get("enabled"),
                    field_name=f"{scope}.activity_config.enabled",
                ),
                threshold=require_float(
                    activity_payload.get("threshold"),
                    field_name=f"{scope}.activity_config.threshold",
                ),
                min_substrates=require_int(
                    activity_payload.get(
                        "min_substrates",
                        KinaseActivityConfig().min_substrates,
                    ),
                    field_name=f"{scope}.activity_config.min_substrates",
                ),
                top_n_substrates=require_int(
                    activity_payload.get(
                        "top_n_substrates",
                        KinaseActivityConfig().top_n_substrates,
                    ),
                    field_name=f"{scope}.activity_config.top_n_substrates",
                ),
            )
        legacy_ensemble_size = prediction_payload.get("ensemble_size")
        deterministic_max_selected_kinases: int
        adaptive_ensemble_runs: int
        if legacy_ensemble_size is not None:
            if (
                prediction_payload.get("deterministic_max_selected_kinases") is not None
                or prediction_payload.get("adaptive_ensemble_runs") is not None
            ):
                raise PhosPyInputError(
                    f"{scope}.prediction_config.ensemble_size cannot be combined with "
                    "deterministic_max_selected_kinases or adaptive_ensemble_runs"
                )
            deterministic_max_selected_kinases = require_int(
                legacy_ensemble_size,
                field_name=f"{scope}.prediction_config.ensemble_size",
            )
            adaptive_ensemble_runs = deterministic_max_selected_kinases
        else:
            deterministic_max_selected_kinases = require_int(
                prediction_payload.get(
                    "deterministic_max_selected_kinases",
                    KINASE_PREDICTION_DEFAULT_DETERMINISTIC_MAX_SELECTED_KINASES,
                ),
                field_name=(
                    f"{scope}.prediction_config.deterministic_max_selected_kinases"
                ),
            )
            adaptive_ensemble_runs = require_int(
                prediction_payload.get(
                    "adaptive_ensemble_runs",
                    KINASE_PREDICTION_DEFAULT_ADAPTIVE_ENSEMBLE_RUNS,
                ),
                field_name=f"{scope}.prediction_config.adaptive_ensemble_runs",
            )
        return cls(
            scoring_config=KinaseScoringConfig(
                min_substrates=require_int(
                    scoring_payload.get("min_substrates"),
                    field_name=f"{scope}.scoring_config.min_substrates",
                ),
                include_diagnostic_scoring_tables=require_bool(
                    scoring_payload.get("include_diagnostic_scoring_tables", True),
                    field_name=(
                        f"{scope}.scoring_config.include_diagnostic_scoring_tables"
                    ),
                ),
                profile_missing_value_strategy=require_str(
                    scoring_payload.get(
                        "profile_missing_value_strategy",
                        KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT,
                    ),
                    field_name=(
                        f"{scope}.scoring_config.profile_missing_value_strategy"
                    ),
                ),
            ),
            prediction_config=KinasePredictionConfig(
                top_k=require_int(
                    prediction_payload.get("top_k"),
                    field_name=f"{scope}.prediction_config.top_k",
                ),
                deterministic_max_selected_kinases=deterministic_max_selected_kinases,
                adaptive_ensemble_runs=adaptive_ensemble_runs,
                mode=require_str(
                    prediction_payload.get(
                        "mode",
                        KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
                    ),
                    field_name=f"{scope}.prediction_config.mode",
                ),
                adaptive_policy=require_str(
                    prediction_payload.get(
                        "adaptive_policy",
                        KINASE_ADAPTIVE_POLICY_STABLE,
                    ),
                    field_name=f"{scope}.prediction_config.adaptive_policy",
                ),
                n_iterations=require_int(
                    prediction_payload.get(
                        "n_iterations",
                        KINASE_PREDICTION_DEFAULT_ITERATIONS,
                    ),
                    field_name=f"{scope}.prediction_config.n_iterations",
                ),
                random_state=(
                    None
                    if prediction_payload.get("random_state") is None
                    else require_int(
                        prediction_payload.get("random_state"),
                        field_name=f"{scope}.prediction_config.random_state",
                    )
                ),
            ),
            activity_config=activity_config,
        )
