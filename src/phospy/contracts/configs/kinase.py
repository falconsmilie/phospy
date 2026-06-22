"""Public kinase workflow configuration models."""
# pyright: reportUnnecessaryIsInstance=false
# Runtime boundary guards are intentionally retained for untyped external callers.

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

from phospy.contracts.configs.common import _require_int_at_least, _require_real_between
from phospy.contracts.configs.localisation import LocalisationRequirement
from phospy.errors.validation import WorkflowValidationError

KINASE_SCORING_MIN_SUBSTRATES_FLOOR = 2
KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED = "phosr_rank_weighted"
KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF = "kinase_library_motif"
KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF = "combined_profile_motif"
KinaseScoringMode = Literal[
    "phosr_rank_weighted",
    "kinase_library_motif",
    "combined_profile_motif",
]
KINASE_SCORING_MODES = frozenset(
    {
        KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
        KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF,
        KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
    }
)
KINASE_SCORING_MODES_REQUIRING_KINASE_LIBRARY = frozenset(
    {
        KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF,
        KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
    }
)
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
KINASE_ALLOW_MIXED_TOTAL_PROTEIN_QUANTITATIVE_MEANING_DEFAULT = False
KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR = "error"
KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE = "prefer_reference"
KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET = "prefer_dataset"
KinaseSiteSequenceConflictPolicy = Literal[
    "error",
    "prefer_reference",
    "prefer_dataset",
]
KINASE_SITE_SEQUENCE_CONFLICT_POLICIES = frozenset(
    {
        KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
        KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
        KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET,
    }
)
KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR = "error"
KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS = (
    "allow_with_diagnostics"
)
KinaseReferenceDisplayAmbiguityPolicy = Literal[
    "error",
    "allow_with_diagnostics",
]
KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICIES = frozenset(
    {
        KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR,
        KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS,
    }
)


@dataclass(frozen=True, slots=True)
class KinaseScoringConfig:
    """Public scoring-stage configuration.

    `min_substrates` is constrained to the public scoring support floor used by
    the supported rewrite contract. The default and minimum public value is `2`.
    A scoring substrate counts only once per kinase after reference projection:
    it must resolve to a dataset `site_key`, remain in the scoring phospho
    matrix, and have usable workflow site-sequence support. Reference-only,
    unmapped, missing, or duplicate substrate-map rows do not increase the
    support count.

    Kinases with fewer than `min_substrates` usable substrates are excluded from
    scoring profiles and downstream scoring columns. Kinases exactly at the
    floor are included, but two-substrate profiles are still minimal evidence
    and should be interpreted cautiously. Single-substrate profiles are not part
    of the public scoring contract because a one-row profile is especially
    sensitive to that substrate's values.

    Supported scoring semantics are stage-pure: score generation is determined
    only by analysis-ready dataset values, resolved reference content, and this
    explicit scoring configuration. Prediction mode and reference input
    provenance (preset vs explicit bundle) do not redefine scoring behavior.

    `scoring_mode` selects the authoritative downstream scoring source. The
    default (`"phosr_rank_weighted"`) preserves the existing PhosR-inspired
    rank-weighted scoring lane. PhosPy builds profiles from available
    substrate/reference evidence, uses motif-frequency support when sequence
    and reference evidence allow, and combines available profile and motif
    evidence with rank-derived weights under the configured substrate/support
    rules. This is a PhosPy-specific scoring mode, not an exact PhosR
    implementation and not a numerical compatibility mode.

    `"kinase_library_motif"` and `"combined_profile_motif"` are workflow-level
    opt-ins for caller-supplied Kinase Library-style resources. They still run
    inside normal kinase workflow interpretation: reference resolution,
    display-ID projection, site-sequence support, and eligible
    kinase-substrate-map context remain required. The local
    `KinaseLibraryResource` supplies motif matrices for workflow support
    scores; it is not a replacement for the workflow reference bundle and does
    not imply official Kinase Library predictor parity.

    `include_diagnostic_scoring_tables` controls publication of non-authoritative
    diagnostic scoring outputs (`motif_scores`, `score_fusion_weights`). The
    authoritative downstream lane for the default mode (`rank_weighted_fusion_scores`
    with profile fallback) is always computed in that mode.

    `profile_missing_value_strategy` controls column-wise median behavior when a
    kinase profile is built from multiple quantified substrates:

    - `"strict"` propagates missing values (`median(..., skipna=False)`)
    - `"median_skipna"` ignores missing values (`median(..., skipna=True)`)
    """

    min_substrates: int = KINASE_SCORING_MIN_SUBSTRATES_FLOOR
    scoring_mode: KinaseScoringMode = KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED
    include_diagnostic_scoring_tables: bool = False
    profile_missing_value_strategy: KinaseProfileMissingValueStrategy = (
        KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT
    )
    localisation_requirement: LocalisationRequirement = field(
        default_factory=LocalisationRequirement
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
        if self.scoring_mode not in KINASE_SCORING_MODES:
            allowed = ", ".join(sorted(KINASE_SCORING_MODES))
            raise WorkflowValidationError(
                f"scoring_config.scoring_mode must be one of: {allowed}"
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
        if not isinstance(self.localisation_requirement, LocalisationRequirement):
            raise WorkflowValidationError(
                "scoring_config.localisation_requirement must be "
                "LocalisationRequirement"
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

    Activity support is method-specific and separate from scoring support. The
    default simplified weighted activity floor is `3` predicted substrates; the
    KSEA-style and ssGSEA-style floors default to `5`. Weighted and KSEA-style
    membership uses finite prediction support at or above the configured
    threshold. When a kinase or kinase-condition pair has too few usable
    substrates for the selected method, activity values are omitted or reported
    as not computable through activity diagnostics rather than silently filled.
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
        if not isinstance(self.enabled, bool):
            raise WorkflowValidationError("activity_config.enabled must be a bool")
        if self.method not in KINASE_ACTIVITY_METHODS:
            allowed_methods = ", ".join(sorted(KINASE_ACTIVITY_METHODS))
            raise WorkflowValidationError(
                f"activity_config.method must be one of: {allowed_methods}"
            )
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
        _require_int_at_least(
            self.ksea_min_substrates,
            field_name="activity_config.ksea_min_substrates",
            minimum=KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR,
            error_type=WorkflowValidationError,
        )
        if self.ksea_evidence_threshold is not None:
            _require_real_between(
                self.ksea_evidence_threshold,
                field_name="activity_config.ksea_evidence_threshold",
                minimum=0.0,
                maximum=1.0,
                error_type=WorkflowValidationError,
            )
        if self.ksea_p_value_method not in KINASE_ACTIVITY_KSEA_P_VALUE_METHODS:
            allowed_methods = ", ".join(sorted(KINASE_ACTIVITY_KSEA_P_VALUE_METHODS))
            raise WorkflowValidationError(
                f"activity_config.ksea_p_value_method must be one of: {allowed_methods}"
            )
        if not isinstance(self.ksea_adjust_p_values, bool):
            raise WorkflowValidationError(
                "activity_config.ksea_adjust_p_values must be a bool"
            )
        _require_int_at_least(
            self.ssgsea_min_substrates,
            field_name="activity_config.ssgsea_min_substrates",
            minimum=KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR,
            error_type=WorkflowValidationError,
        )
        if (
            self.ssgsea_ranking_direction
            not in KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTIONS
        ):
            allowed = ", ".join(sorted(KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTIONS))
            raise WorkflowValidationError(
                f"activity_config.ssgsea_ranking_direction must be one of: {allowed}"
            )
        _require_int_at_least(
            self.ssgsea_permutations,
            field_name="activity_config.ssgsea_permutations",
            minimum=0,
            error_type=WorkflowValidationError,
        )
        if self.ssgsea_random_seed is not None:
            _require_int_at_least(
                self.ssgsea_random_seed,
                field_name="activity_config.ssgsea_random_seed",
                minimum=0,
                error_type=WorkflowValidationError,
            )
        if self.ssgsea_permutations > 0 and self.ssgsea_random_seed is None:
            raise WorkflowValidationError(
                "activity_config.ssgsea_random_seed must be set when "
                "activity_config.ssgsea_permutations is greater than 0"
            )
        if not isinstance(self.ssgsea_adjust_p_values, bool):
            raise WorkflowValidationError(
                "activity_config.ssgsea_adjust_p_values must be a bool"
            )


__all__ = [
    "KINASE_ALLOW_MIXED_TOTAL_PROTEIN_QUANTITATIVE_MEANING_DEFAULT",
    "KINASE_SITE_SEQUENCE_CONFLICT_POLICIES",
    "KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR",
    "KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET",
    "KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE",
    "KINASE_ACTIVITY_KSEA_DEFAULT_ADJUST_P_VALUES",
    "KINASE_ACTIVITY_KSEA_DEFAULT_MIN_SUBSTRATES",
    "KINASE_ACTIVITY_KSEA_P_VALUE_METHOD_NORMAL_APPROXIMATION",
    "KINASE_ACTIVITY_KSEA_P_VALUE_METHODS",
    "KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES",
    "KINASE_ACTIVITY_DEFAULT_THRESHOLD",
    "KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES",
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
    "KINASE_PROFILE_MISSING_VALUE_STRATEGIES",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGY_MEDIAN_SKIPNA",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT",
    "KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICIES",
    "KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS",
    "KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR",
    "KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF",
    "KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF",
    "KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED",
    "KINASE_SCORING_MODES",
    "KINASE_SCORING_MODES_REQUIRING_KINASE_LIBRARY",
    "KINASE_SCORING_MIN_SUBSTRATES_FLOOR",
    "KinaseActivityConfig",
    "KinaseActivityMethod",
    "KinaseActivityPValueMethod",
    "KinaseActivitySsgseaRankingDirection",
    "KinaseProfileMissingValueStrategy",
    "KinaseReferenceDisplayAmbiguityPolicy",
    "KinaseScoringMode",
    "KinaseSiteSequenceConflictPolicy",
    "LocalisationRequirement",
    "KinaseScoringConfig",
]
