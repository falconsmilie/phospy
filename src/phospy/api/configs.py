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
SIGNALOME_MODULE_COUNT_FLOOR = 1
SIGNALOME_MODULE_SELECTION_PRIMARY_THRESHOLD_DEFAULT = 0.5
SIGNALOME_MODULE_SELECTION_FALLBACK_THRESHOLD_DEFAULT = 0.1
SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_FLOOR = 1
SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_DEFAULT = 10
SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY = "cutoff_binary"
SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP = "weighted_top"
SignalomeAssignmentPolicy = Literal["cutoff_binary", "weighted_top"]
SIGNALOME_ASSIGNMENT_POLICIES = frozenset(
    {
        SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
        SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
    }
)
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
    """Public prediction-stage configuration.

    `mode` selects the prediction lane:

    - `"deterministic_ranking"`: deterministic top-kinase selection from
      downstream scores (legacy rewrite shortcut lane).
    - `"adaptive_ensemble"`: real adaptive ensemble execution ported from donor
      science.

    `ensemble_size` is mode-dependent by design and should be interpreted with
    `mode`:

    - deterministic lane: maximum number of selected kinase columns.
    - adaptive lane: number of ensemble executions per kinase.
    """

    top_k: int = 30
    ensemble_size: int = 10
    mode: KinasePredictionMode = KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING
    adaptive_policy: KinaseAdaptivePolicy = KINASE_ADAPTIVE_POLICY_STABLE
    n_iterations: int = KINASE_PREDICTION_DEFAULT_ITERATIONS
    random_state: int | None = None


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
    """Public signalome workflow configuration.

    `assignment_policy` controls module-support attribution:

    - `"cutoff_binary"`: binary support per kinase from
      `substrate_support_cutoff`.
    - `"weighted_top"`: fractional support propagated from per-site
      `top_kinase_weights` ties.
    """

    substrate_support_cutoff: float = 0.5
    network_correlation_threshold: float = 0.5
    assignment_policy: SignalomeAssignmentPolicy = (
        SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY
    )
    module_count: int | None = None
    module_selection_primary_correlation_threshold: float = (
        SIGNALOME_MODULE_SELECTION_PRIMARY_THRESHOLD_DEFAULT
    )
    module_selection_fallback_correlation_threshold: float = (
        SIGNALOME_MODULE_SELECTION_FALLBACK_THRESHOLD_DEFAULT
    )
    module_selection_max_clusters: int = SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_DEFAULT


__all__ = [
    "KINASE_ADAPTIVE_POLICIES",
    "KINASE_ADAPTIVE_POLICY_R_PARITY",
    "KINASE_ADAPTIVE_POLICY_STABLE",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGIES",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGY_MEDIAN_SKIPNA",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT",
    "KINASE_PREDICTION_DEFAULT_ITERATIONS",
    "KINASE_PREDICTION_MODES",
    "KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE",
    "KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING",
    "KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES",
    "KINASE_ACTIVITY_DEFAULT_THRESHOLD",
    "KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES",
    "KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR",
    "KINASE_ACTIVITY_TOP_N_SUBSTRATES_FLOOR",
    "KINASE_SCORING_MIN_SUBSTRATES_FLOOR",
    "SIGNALOME_ASSIGNMENT_POLICIES",
    "SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY",
    "SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP",
    "SIGNALOME_MODULE_COUNT_FLOOR",
    "SIGNALOME_MODULE_SELECTION_FALLBACK_THRESHOLD_DEFAULT",
    "SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_DEFAULT",
    "SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_FLOOR",
    "SIGNALOME_MODULE_SELECTION_PRIMARY_THRESHOLD_DEFAULT",
    "KinaseAdaptivePolicy",
    "KinaseProfileMissingValueStrategy",
    "KinasePredictionMode",
    "KinaseActivityConfig",
    "KinasePredictionConfig",
    "KinaseScoringConfig",
    "SignalomeAssignmentPolicy",
    "SignalomeConfig",
]
