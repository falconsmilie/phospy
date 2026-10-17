"""Kinase activity configuration and execution contracts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from phospy.activities.contracts import ResolvedMethodQuantitativeInputContract
from phospy.errors.validation import ContractValidationError
from phospy.kinase.scoring.config import _require_bool
from phospy.policies.config_values import (
    require_int_at_least as _require_int_at_least,
)
from phospy.policies.config_values import (
    require_real_between as _require_real_between,
)

KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR = 1
KINASE_ACTIVITY_TOP_N_SUBSTRATES_FLOOR = 1
KINASE_ACTIVITY_DEFAULT_THRESHOLD = 0.6
KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES = 3
KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES = 20
KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY = (
    "simplified_weighted_substrate_activity"
)
KINASE_ACTIVITY_METHOD_KSEA_ZSCORE = "ksea_zscore"
KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT = "ssgsea_substrate_enrichment"
KinaseActivityMethod = Literal[
    "simplified_weighted_substrate_activity",
    "ksea_zscore",
    "ssgsea_substrate_enrichment",
]
KINASE_ACTIVITY_METHODS = frozenset(
    {
        KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
        KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
        KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
    }
)
KINASE_ACTIVITY_KSEA_DEFAULT_MIN_SUBSTRATES = 5
KINASE_ACTIVITY_KSEA_P_VALUE_METHOD_NORMAL_APPROXIMATION = "normal_approximation"
KinaseActivityPValueMethod = Literal["normal_approximation"]
KINASE_ACTIVITY_KSEA_P_VALUE_METHODS = frozenset(
    {KINASE_ACTIVITY_KSEA_P_VALUE_METHOD_NORMAL_APPROXIMATION}
)
KINASE_ACTIVITY_KSEA_DEFAULT_ADJUST_P_VALUES = True
KINASE_ACTIVITY_SSGSEA_DEFAULT_MIN_SUBSTRATES = 5
KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTION_DESCENDING = "descending"
KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTION_ASCENDING = "ascending"
KinaseActivitySsgseaRankingDirection = Literal["descending", "ascending"]
KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTIONS = frozenset(
    {
        KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTION_DESCENDING,
        KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTION_ASCENDING,
    }
)
KINASE_ACTIVITY_SSGSEA_DEFAULT_PERMUTATIONS = 0
KINASE_ACTIVITY_SSGSEA_DEFAULT_RANDOM_SEED = 0
KINASE_ACTIVITY_SSGSEA_DEFAULT_ADJUST_P_VALUES = True


@dataclass(frozen=True, slots=True)
class KinaseActivityConfig:
    """Configuration for the supported kinase activity score stage.

    Activity-like scoring runs inside `KinaseWorkflow` and can be disabled by
    setting either:

    - `activity_config=None` on `KinaseWorkflowRequest`, or
    - `activity_config.enabled=False`.

    These outputs are exploratory kinase activity scores or activity-like
    substrate summaries, not direct proof that a kinase is causally active.
    Support is method-specific and separate from scoring support. The default
    simplified weighted score floor is `3` predicted substrates; the KSEA-style
    and ssGSEA-style floors default to `5`. Weighted and KSEA-style membership
    uses finite prediction support at or above the configured threshold. Missing
    or sparse substrate support weakens interpretation. When a kinase or
    kinase-profile pair has too few usable substrates for the selected method,
    score values are omitted or reported as not computable through diagnostics
    rather than silently filled.
    """

    enabled: bool = True
    method: KinaseActivityMethod = (
        KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY
    )
    threshold: float = KINASE_ACTIVITY_DEFAULT_THRESHOLD
    min_substrates: int = KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES
    top_n_substrates: int = KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES
    ksea_min_substrates: int = KINASE_ACTIVITY_KSEA_DEFAULT_MIN_SUBSTRATES
    ksea_evidence_threshold: float | None = None
    ksea_p_value_method: KinaseActivityPValueMethod = (
        KINASE_ACTIVITY_KSEA_P_VALUE_METHOD_NORMAL_APPROXIMATION
    )
    ksea_adjust_p_values: bool = KINASE_ACTIVITY_KSEA_DEFAULT_ADJUST_P_VALUES
    ssgsea_min_substrates: int = KINASE_ACTIVITY_SSGSEA_DEFAULT_MIN_SUBSTRATES
    ssgsea_ranking_direction: KinaseActivitySsgseaRankingDirection = (
        KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTION_DESCENDING
    )
    ssgsea_permutations: int = KINASE_ACTIVITY_SSGSEA_DEFAULT_PERMUTATIONS
    ssgsea_random_seed: int | None = KINASE_ACTIVITY_SSGSEA_DEFAULT_RANDOM_SEED
    ssgsea_adjust_p_values: bool = KINASE_ACTIVITY_SSGSEA_DEFAULT_ADJUST_P_VALUES

    def __post_init__(self) -> None:
        enabled = _require_bool(
            self.enabled,
            field_name="activity_config.enabled",
        )
        if self.method not in KINASE_ACTIVITY_METHODS:
            allowed_methods = ", ".join(sorted(KINASE_ACTIVITY_METHODS))
            raise ContractValidationError(
                f"activity_config.method must be one of: {allowed_methods}"
            )
        _require_real_between(
            self.threshold,
            field_name="activity_config.threshold",
            minimum=0.0,
            maximum=1.0,
            error_type=ContractValidationError,
        )
        _require_int_at_least(
            self.min_substrates,
            field_name="activity_config.min_substrates",
            minimum=KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR,
            error_type=ContractValidationError,
        )
        _require_int_at_least(
            self.top_n_substrates,
            field_name="activity_config.top_n_substrates",
            minimum=KINASE_ACTIVITY_TOP_N_SUBSTRATES_FLOOR,
            error_type=ContractValidationError,
        )
        _require_int_at_least(
            self.ksea_min_substrates,
            field_name="activity_config.ksea_min_substrates",
            minimum=KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR,
            error_type=ContractValidationError,
        )
        if self.ksea_evidence_threshold is not None:
            _require_real_between(
                self.ksea_evidence_threshold,
                field_name="activity_config.ksea_evidence_threshold",
                minimum=0.0,
                maximum=1.0,
                error_type=ContractValidationError,
            )
        if self.ksea_p_value_method not in KINASE_ACTIVITY_KSEA_P_VALUE_METHODS:
            allowed_methods = ", ".join(sorted(KINASE_ACTIVITY_KSEA_P_VALUE_METHODS))
            raise ContractValidationError(
                f"activity_config.ksea_p_value_method must be one of: {allowed_methods}"
            )
        ksea_adjust_p_values = _require_bool(
            self.ksea_adjust_p_values,
            field_name="activity_config.ksea_adjust_p_values",
        )
        _require_int_at_least(
            self.ssgsea_min_substrates,
            field_name="activity_config.ssgsea_min_substrates",
            minimum=KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR,
            error_type=ContractValidationError,
        )
        if (
            self.ssgsea_ranking_direction
            not in KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTIONS
        ):
            allowed = ", ".join(sorted(KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTIONS))
            raise ContractValidationError(
                f"activity_config.ssgsea_ranking_direction must be one of: {allowed}"
            )
        _require_int_at_least(
            self.ssgsea_permutations,
            field_name="activity_config.ssgsea_permutations",
            minimum=0,
            error_type=ContractValidationError,
        )
        if self.ssgsea_random_seed is not None:
            _require_int_at_least(
                self.ssgsea_random_seed,
                field_name="activity_config.ssgsea_random_seed",
                minimum=0,
                error_type=ContractValidationError,
            )
        if self.ssgsea_permutations > 0 and self.ssgsea_random_seed is None:
            raise ContractValidationError(
                "activity_config.ssgsea_random_seed must be set when "
                "activity_config.ssgsea_permutations is greater than 0"
            )
        ssgsea_adjust_p_values = _require_bool(
            self.ssgsea_adjust_p_values,
            field_name="activity_config.ssgsea_adjust_p_values",
        )
        object.__setattr__(self, "enabled", enabled)
        object.__setattr__(self, "ksea_adjust_p_values", ksea_adjust_p_values)
        object.__setattr__(self, "ssgsea_adjust_p_values", ssgsea_adjust_p_values)

    @classmethod
    def ssgsea_with_permutation_significance(
        cls,
        *,
        permutations: int,
        random_seed: int,
        ssgsea_min_substrates: int = KINASE_ACTIVITY_SSGSEA_DEFAULT_MIN_SUBSTRATES,
        ssgsea_ranking_direction: KinaseActivitySsgseaRankingDirection = (
            KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTION_DESCENDING
        ),
        ssgsea_adjust_p_values: bool = KINASE_ACTIVITY_SSGSEA_DEFAULT_ADJUST_P_VALUES,
        enabled: bool = True,
    ) -> KinaseActivityConfig:
        """Return ssGSEA config with explicit seeded permutation p-values."""

        config = cls(
            enabled=enabled,
            method=KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
            ssgsea_min_substrates=ssgsea_min_substrates,
            ssgsea_ranking_direction=ssgsea_ranking_direction,
            ssgsea_permutations=permutations,
            ssgsea_random_seed=random_seed,
            ssgsea_adjust_p_values=ssgsea_adjust_p_values,
        )
        if int(config.ssgsea_permutations) <= 0:
            raise ContractValidationError(
                "activity_config.ssgsea_permutations must be greater than 0 "
                "when requesting ssGSEA permutation significance"
            )
        return config


@dataclass(frozen=True, slots=True)
class ResolvedKinaseActivityExecutionConfig:
    """Execution-ready kinase activity-stage config."""

    method: KinaseActivityMethod
    threshold: float
    min_substrates: int
    top_n_substrates: int
    ksea_min_substrates: int
    ksea_evidence_threshold: float
    ksea_p_value_method: KinaseActivityPValueMethod
    ksea_adjust_p_values: bool
    ssgsea_min_substrates: int = 5
    ssgsea_ranking_direction: KinaseActivitySsgseaRankingDirection = "descending"
    ssgsea_permutations: int = 0
    ssgsea_random_seed: int | None = 0
    ssgsea_adjust_p_values: bool = True
    method_input_contract: ResolvedMethodQuantitativeInputContract | None = None


def resolve_kinase_activity_execution_config(
    config: KinaseActivityConfig | None,
    *,
    method_input_contract: ResolvedMethodQuantitativeInputContract | None,
) -> ResolvedKinaseActivityExecutionConfig | None:
    """Resolve an optional public activity config for stage execution."""

    if config is None or not config.enabled:
        return None
    return ResolvedKinaseActivityExecutionConfig(
        method=config.method,
        threshold=float(config.threshold),
        min_substrates=int(config.min_substrates),
        top_n_substrates=int(config.top_n_substrates),
        ksea_min_substrates=int(config.ksea_min_substrates),
        ksea_evidence_threshold=float(
            config.threshold
            if config.ksea_evidence_threshold is None
            else config.ksea_evidence_threshold
        ),
        ksea_p_value_method=config.ksea_p_value_method,
        ksea_adjust_p_values=bool(config.ksea_adjust_p_values),
        ssgsea_min_substrates=int(config.ssgsea_min_substrates),
        ssgsea_ranking_direction=config.ssgsea_ranking_direction,
        ssgsea_permutations=int(config.ssgsea_permutations),
        ssgsea_random_seed=(
            None
            if config.ssgsea_random_seed is None
            else int(config.ssgsea_random_seed)
        ),
        ssgsea_adjust_p_values=bool(config.ssgsea_adjust_p_values),
        method_input_contract=method_input_contract,
    )


__all__ = [
    "KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES",
    "KINASE_ACTIVITY_DEFAULT_THRESHOLD",
    "KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES",
    "KINASE_ACTIVITY_KSEA_DEFAULT_ADJUST_P_VALUES",
    "KINASE_ACTIVITY_KSEA_DEFAULT_MIN_SUBSTRATES",
    "KINASE_ACTIVITY_KSEA_P_VALUE_METHOD_NORMAL_APPROXIMATION",
    "KINASE_ACTIVITY_KSEA_P_VALUE_METHODS",
    "KINASE_ACTIVITY_METHODS",
    "KINASE_ACTIVITY_METHOD_KSEA_ZSCORE",
    "KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY",
    "KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT",
    "KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR",
    "KINASE_ACTIVITY_SSGSEA_DEFAULT_ADJUST_P_VALUES",
    "KINASE_ACTIVITY_SSGSEA_DEFAULT_MIN_SUBSTRATES",
    "KINASE_ACTIVITY_SSGSEA_DEFAULT_PERMUTATIONS",
    "KINASE_ACTIVITY_SSGSEA_DEFAULT_RANDOM_SEED",
    "KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTION_ASCENDING",
    "KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTION_DESCENDING",
    "KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTIONS",
    "KINASE_ACTIVITY_TOP_N_SUBSTRATES_FLOOR",
    "KinaseActivityConfig",
    "KinaseActivityMethod",
    "KinaseActivityPValueMethod",
    "KinaseActivitySsgseaRankingDirection",
    "ResolvedKinaseActivityExecutionConfig",
    "resolve_kinase_activity_execution_config",
]
