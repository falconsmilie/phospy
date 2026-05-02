"""Public kinase workflow configuration models."""
# pyright: reportUnnecessaryIsInstance=false
# Runtime boundary guards are intentionally retained for untyped external callers.

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from phospy.api.configs.common import _require_int_at_least, _require_real_between
from phospy.errors.validation import WorkflowValidationError

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
KINASE_ALLOW_MIXED_TOTAL_PROTEIN_QUANTITATIVE_MEANING_DEFAULT = False


@dataclass(frozen=True, slots=True)
class KinaseScoringConfig:
    """Public scoring-stage configuration.

    `min_substrates` is constrained to the public scoring support floor used by
    the supported rewrite contract.

    Supported scoring semantics are stage-pure: score generation is determined
    only by analysis-ready dataset values, resolved reference content, and this
    explicit scoring configuration. Prediction mode and reference input
    provenance (preset vs explicit bundle) do not redefine scoring behavior.

    `include_diagnostic_scoring_tables` controls publication of non-authoritative
    diagnostic scoring outputs (`motif_scores`, `score_fusion_weights`). The
    authoritative downstream lane (`rank_weighted_fusion_scores` with profile
    fallback) is always computed.

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
    allow_mixed_total_protein_quantitative_meaning: bool = (
        KINASE_ALLOW_MIXED_TOTAL_PROTEIN_QUANTITATIVE_MEANING_DEFAULT
    )

    def __post_init__(self) -> None:
        if not isinstance(self.include_diagnostic_scoring_tables, bool):
            raise WorkflowValidationError(
                "scoring_config.include_diagnostic_scoring_tables must be a bool"
            )
        if not isinstance(self.allow_mixed_total_protein_quantitative_meaning, bool):
            raise WorkflowValidationError(
                "scoring_config.allow_mixed_total_protein_quantitative_meaning "
                "must be a bool"
            )
        if (
            self.profile_missing_value_strategy
            not in KINASE_PROFILE_MISSING_VALUE_STRATEGIES
        ):
            allowed = ", ".join(sorted(KINASE_PROFILE_MISSING_VALUE_STRATEGIES))
            raise WorkflowValidationError(
                "scoring_config.profile_missing_value_strategy must be one of: "
                f"{allowed}"
            )
        _require_int_at_least(
            self.min_substrates,
            field_name="scoring_config.min_substrates",
            minimum=KINASE_SCORING_MIN_SUBSTRATES_FLOOR,
            error_type=WorkflowValidationError,
        )

    @classmethod
    def default(cls) -> KinaseScoringConfig:
        """Return the package default kinase scoring profile."""
        return cls()

    @classmethod
    def strict_missing_values(cls) -> KinaseScoringConfig:
        """Return strict missing-value handling for profile aggregation."""
        return cls(
            profile_missing_value_strategy=(
                KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT
            )
        )


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

    def __post_init__(self) -> None:
        if not isinstance(self.enabled, bool):
            raise WorkflowValidationError("activity_config.enabled must be a bool")
        _require_real_between(
            self.threshold,
            field_name="activity_config.threshold",
            minimum=0.0,
            maximum=1.0,
            error_type=WorkflowValidationError,
        )
        _require_int_at_least(
            self.min_substrates,
            field_name="activity_config.min_substrates",
            minimum=KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR,
            error_type=WorkflowValidationError,
        )
        _require_int_at_least(
            self.top_n_substrates,
            field_name="activity_config.top_n_substrates",
            minimum=KINASE_ACTIVITY_TOP_N_SUBSTRATES_FLOOR,
            error_type=WorkflowValidationError,
        )


__all__ = [
    "KINASE_ALLOW_MIXED_TOTAL_PROTEIN_QUANTITATIVE_MEANING_DEFAULT",
    "KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES",
    "KINASE_ACTIVITY_DEFAULT_THRESHOLD",
    "KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES",
    "KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR",
    "KINASE_ACTIVITY_TOP_N_SUBSTRATES_FLOOR",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGIES",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGY_MEDIAN_SKIPNA",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT",
    "KINASE_SCORING_MIN_SUBSTRATES_FLOOR",
    "KinaseActivityConfig",
    "KinaseProfileMissingValueStrategy",
    "KinaseScoringConfig",
]
