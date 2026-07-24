"""Public kinase workflow configuration models."""
# pyright: reportUnnecessaryIsInstance=false
# Runtime boundary guards are intentionally retained for untyped external callers.

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import cast

from phospy.contracts.configs._validation import coerce_policy_enum
from phospy.contracts.configs.common import _require_int_at_least, _require_real_between
from phospy.contracts.configs.localisation import LocalisationRequirement
from phospy.contracts.configs.reference_context import (
    REFERENCE_CONTEXT_COMPATIBILITY_POLICY_REQUIRE_KNOWN_MATCH,
    ReferenceContextCompatibilityPolicy,
)
from phospy.errors.validation import ContractValidationError
from phospy.science.configs.kinase import (
    KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES,
    KINASE_ACTIVITY_DEFAULT_THRESHOLD,
    KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES,
    KINASE_ACTIVITY_KSEA_DEFAULT_ADJUST_P_VALUES,
    KINASE_ACTIVITY_KSEA_DEFAULT_MIN_SUBSTRATES,
    KINASE_ACTIVITY_KSEA_P_VALUE_METHOD_NORMAL_APPROXIMATION,
    KINASE_ACTIVITY_KSEA_P_VALUE_METHODS,
    KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
    KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
    KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
    KINASE_ACTIVITY_METHODS,
    KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR,
    KINASE_ACTIVITY_SSGSEA_DEFAULT_ADJUST_P_VALUES,
    KINASE_ACTIVITY_SSGSEA_DEFAULT_MIN_SUBSTRATES,
    KINASE_ACTIVITY_SSGSEA_DEFAULT_PERMUTATIONS,
    KINASE_ACTIVITY_SSGSEA_DEFAULT_RANDOM_SEED,
    KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTION_ASCENDING,
    KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTION_DESCENDING,
    KINASE_ACTIVITY_SSGSEA_RANKING_DIRECTIONS,
    KINASE_ACTIVITY_TOP_N_SUBSTRATES_FLOOR,
    KINASE_ALLOW_MIXED_TOTAL_PROTEIN_QUANTITATIVE_MEANING_DEFAULT,
    KINASE_ATTRITION_POLICY_ON_VIOLATION_ERROR,
    KINASE_ATTRITION_POLICY_ON_VIOLATION_MODES,
    KINASE_ATTRITION_POLICY_ON_VIOLATION_WARN,
    KINASE_LIBRARY_MOTIF_ALIAS_DEPRECATION_MESSAGE,
    KINASE_PROFILE_MISSING_VALUE_STRATEGIES,
    KINASE_PROFILE_MISSING_VALUE_STRATEGY_MEDIAN_SKIPNA,
    KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT,
    KINASE_PROFILE_SELF_INCLUSION_POLICIES,
    KINASE_PROFILE_SELF_INCLUSION_POLICY_ALLOW,
    KINASE_PROFILE_SELF_INCLUSION_POLICY_LEAVE_ONE_OUT,
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICIES,
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS,
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR,
    KINASE_RELIABILITY_PROFILE_CUSTOM,
    KINASE_RELIABILITY_PROFILE_EXPLORATORY,
    KINASE_RELIABILITY_PROFILE_PRODUCTION,
    KINASE_RELIABILITY_PROFILES,
    KINASE_SCORING_MIN_SUBSTRATES_FLOOR,
    KINASE_SCORING_MODE_ALIASES,
    KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
    KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF,
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
    KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
    KINASE_SCORING_MODES,
    KINASE_SCORING_MODES_REQUIRING_KINASE_LIBRARY,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICIES,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
    KinaseActivityMethod,
    KinaseActivityPValueMethod,
    KinaseActivitySsgseaRankingDirection,
    KinaseAttritionViolationMode,
    KinaseProfileMissingValueStrategy,
    KinaseReferenceDisplayAmbiguityPolicy,
    KinaseReliabilityProfile,
    KinaseScoringMode,
    KinaseSiteSequenceConflictPolicy,
    ProfileSelfInclusionPolicy,
    normalize_kinase_scoring_mode,
)

_UNSET = object()


@dataclass(frozen=True, slots=True)
class KinaseAttritionPolicy:
    """Policy thresholds for reporting unacceptable kinase site attrition.

    This is an explicit reliability policy only. It records caller thresholds
    for reference-overlap, sequence-support, and final scored-site retention,
    plus whether threshold violations fail the workflow or continue with
    structured caveats. The workflow never lowers these thresholds or drops
    additional sites to satisfy them.
    """

    minimum_reference_overlap_fraction: float = 0.0
    minimum_sequence_supported_fraction: float = 0.0
    minimum_scored_fraction: float = 0.0
    on_violation: KinaseAttritionViolationMode = (
        KINASE_ATTRITION_POLICY_ON_VIOLATION_WARN
    )

    def __post_init__(self) -> None:
        reference_overlap = _require_attrition_fraction(
            self.minimum_reference_overlap_fraction,
            field_name="attrition_policy.minimum_reference_overlap_fraction",
        )
        sequence_supported = _require_attrition_fraction(
            self.minimum_sequence_supported_fraction,
            field_name="attrition_policy.minimum_sequence_supported_fraction",
        )
        scored = _require_attrition_fraction(
            self.minimum_scored_fraction,
            field_name="attrition_policy.minimum_scored_fraction",
        )
        if self.on_violation not in KINASE_ATTRITION_POLICY_ON_VIOLATION_MODES:
            allowed = ", ".join(sorted(KINASE_ATTRITION_POLICY_ON_VIOLATION_MODES))
            raise ContractValidationError(
                f"attrition_policy.on_violation must be one of: {allowed}"
            )
        object.__setattr__(
            self,
            "minimum_reference_overlap_fraction",
            reference_overlap,
        )
        object.__setattr__(
            self,
            "minimum_sequence_supported_fraction",
            sequence_supported,
        )
        object.__setattr__(self, "minimum_scored_fraction", scored)


@dataclass(frozen=True, slots=True, init=False)
class KinaseScoringConfig:
    """Public scoring-stage configuration.

    ``reliability_profile`` exposes whether scoring is using the explicit
    exploratory preset, the stricter production preset, or custom caller values.
    Direct construction without ``reliability_profile`` keeps historical
    numerical defaults and is classified from the supplied values: exact old
    defaults resolve to ``EXPLORATORY`` and modified values resolve to
    ``CUSTOM``. ``PRODUCTION`` is never inferred from values alone.

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

    `"kinase_library_motif_only"` is a motif-only scoring mode for
    caller-supplied Kinase Library-style resources. It uses workflow-validated
    site identity and centered sequence context but does not construct
    substrate-derived kinase profiles or require quantified known-substrate
    profile overlap.

    `"kinase_library_contextual_motif"` and `"combined_profile_motif"` are
    workflow-level opt-ins for caller-supplied Kinase Library-style resources.
    They still run inside normal kinase workflow interpretation: reference
    resolution, display-ID projection, site-sequence support, and eligible
    kinase-substrate-map context remain required. The local
    `KinaseLibraryResource` supplies motif matrices for workflow support
    scores; it is not a replacement for the workflow reference bundle and does
    not imply official Kinase Library predictor parity. The legacy
    `"kinase_library_motif"` string is accepted as a deprecated alias for
    `"kinase_library_contextual_motif"` during migration.

    `include_diagnostic_scoring_tables` controls publication of non-authoritative
    diagnostic scoring outputs (`motif_scores`, `score_fusion_weights`). The
    authoritative downstream lane for the default mode (`rank_weighted_fusion_scores`
    with profile fallback) is always computed in that mode.

    `include_substrate_contributions` controls assembly and publication of
    optional substrate-level contribution rows on the workflow result. It
    defaults to `False` so routine runs do not build a large evidence table.

    `attrition_policy` records caller-defined minimum retained-site fractions
    for reference overlap, sequence support, and final scoring. The interpreter
    evaluates these thresholds after reference projection and sequence support
    are known, before scoring starts.

    `profile_missing_value_strategy` controls column-wise median behavior when a
    kinase profile is built from multiple quantified substrates:

    - `"strict"` propagates missing values (`median(..., skipna=False)`)
    - `"median_skipna"` ignores missing values (`median(..., skipna=True)`)

    `profile_self_inclusion_policy` declares whether a known substrate site is
    allowed to contribute to the kinase profile used to score that same site.
    The default (`"allow"`) preserves historical scoring behavior. The
    `"leave_one_out"` opt-in recomputes a site's profile score without that
    site when it is part of the kinase's quantified substrate profile.

    `reference_context_compatibility_policy` controls how the workflow handles
    missing dataset/reference biological reference context. The default requires
    known matching contexts. The explicit `"allow_unknown_with_caveat"` override
    permits unknown context only with a result caveat; mismatched known contexts
    are always rejected.
    """

    min_substrates: int = KINASE_SCORING_MIN_SUBSTRATES_FLOOR
    scoring_mode: KinaseScoringMode = KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED
    include_diagnostic_scoring_tables: bool = False
    include_substrate_contributions: bool = False
    profile_missing_value_strategy: KinaseProfileMissingValueStrategy = (
        KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT
    )
    profile_self_inclusion_policy: ProfileSelfInclusionPolicy = (
        ProfileSelfInclusionPolicy.ALLOW
    )
    attrition_policy: KinaseAttritionPolicy = field(
        default_factory=KinaseAttritionPolicy
    )
    localisation_requirement: LocalisationRequirement = field(
        default_factory=LocalisationRequirement
    )
    reference_context_compatibility_policy: ReferenceContextCompatibilityPolicy = (
        REFERENCE_CONTEXT_COMPATIBILITY_POLICY_REQUIRE_KNOWN_MATCH
    )
    allow_mixed_total_protein_quantitative_meaning: bool = (
        KINASE_ALLOW_MIXED_TOTAL_PROTEIN_QUANTITATIVE_MEANING_DEFAULT
    )
    reliability_profile: KinaseReliabilityProfile = field(init=False)
    requested_reliability_profile: KinaseReliabilityProfile | None = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __init__(
        self,
        min_substrates: int = KINASE_SCORING_MIN_SUBSTRATES_FLOOR,
        scoring_mode: KinaseScoringMode = KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
        include_diagnostic_scoring_tables: bool = False,
        include_substrate_contributions: bool = False,
        profile_missing_value_strategy: KinaseProfileMissingValueStrategy = (
            KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT
        ),
        profile_self_inclusion_policy: ProfileSelfInclusionPolicy = (
            ProfileSelfInclusionPolicy.ALLOW
        ),
        attrition_policy: KinaseAttritionPolicy | object = _UNSET,
        localisation_requirement: LocalisationRequirement | object = _UNSET,
        reference_context_compatibility_policy: ReferenceContextCompatibilityPolicy = (
            REFERENCE_CONTEXT_COMPATIBILITY_POLICY_REQUIRE_KNOWN_MATCH
        ),
        allow_mixed_total_protein_quantitative_meaning: bool = (
            KINASE_ALLOW_MIXED_TOTAL_PROTEIN_QUANTITATIVE_MEANING_DEFAULT
        ),
        *,
        reliability_profile: KinaseReliabilityProfile | str | None | object = _UNSET,
    ) -> None:
        object.__setattr__(self, "min_substrates", min_substrates)
        object.__setattr__(self, "scoring_mode", scoring_mode)
        object.__setattr__(
            self,
            "include_diagnostic_scoring_tables",
            include_diagnostic_scoring_tables,
        )
        object.__setattr__(
            self,
            "include_substrate_contributions",
            include_substrate_contributions,
        )
        object.__setattr__(
            self,
            "profile_missing_value_strategy",
            profile_missing_value_strategy,
        )
        object.__setattr__(
            self,
            "profile_self_inclusion_policy",
            profile_self_inclusion_policy,
        )
        object.__setattr__(
            self,
            "attrition_policy",
            (
                KinaseAttritionPolicy()
                if attrition_policy is _UNSET
                else attrition_policy
            ),
        )
        object.__setattr__(
            self,
            "localisation_requirement",
            (
                LocalisationRequirement()
                if localisation_requirement is _UNSET
                else localisation_requirement
            ),
        )
        object.__setattr__(
            self,
            "reference_context_compatibility_policy",
            reference_context_compatibility_policy,
        )
        object.__setattr__(
            self,
            "allow_mixed_total_protein_quantitative_meaning",
            allow_mixed_total_protein_quantitative_meaning,
        )
        self._validate_and_resolve_profile(reliability_profile=reliability_profile)

    def _validate_and_resolve_profile(
        self,
        *,
        reliability_profile: KinaseReliabilityProfile | str | None | object,
    ) -> None:
        if not isinstance(self.include_diagnostic_scoring_tables, bool):
            raise ContractValidationError(
                "scoring_config.include_diagnostic_scoring_tables must be a bool"
            )
        if not isinstance(self.include_substrate_contributions, bool):
            raise ContractValidationError(
                "scoring_config.include_substrate_contributions must be a bool"
            )
        if not isinstance(self.allow_mixed_total_protein_quantitative_meaning, bool):
            raise ContractValidationError(
                "scoring_config.allow_mixed_total_protein_quantitative_meaning "
                "must be a bool"
            )
        normalized_scoring_mode = normalize_kinase_scoring_mode(
            self.scoring_mode,
            warn_on_deprecated_alias=True,
        )
        if normalized_scoring_mode not in KINASE_SCORING_MODES:
            allowed = ", ".join(sorted(KINASE_SCORING_MODES))
            raise ContractValidationError(
                f"scoring_config.scoring_mode must be one of: {allowed}"
            )
        object.__setattr__(self, "scoring_mode", normalized_scoring_mode)
        if (
            self.profile_missing_value_strategy
            not in KINASE_PROFILE_MISSING_VALUE_STRATEGIES
        ):
            allowed = ", ".join(sorted(KINASE_PROFILE_MISSING_VALUE_STRATEGIES))
            raise ContractValidationError(
                "scoring_config.profile_missing_value_strategy must be one of: "
                f"{allowed}"
            )
        profile_self_inclusion_policy = coerce_policy_enum(
            ProfileSelfInclusionPolicy,
            self.profile_self_inclusion_policy,
            field_name="scoring_config.profile_self_inclusion_policy",
            error_type=ContractValidationError,
        )
        object.__setattr__(
            self,
            "profile_self_inclusion_policy",
            profile_self_inclusion_policy,
        )
        reference_context_compatibility_policy = coerce_policy_enum(
            ReferenceContextCompatibilityPolicy,
            self.reference_context_compatibility_policy,
            field_name="scoring_config.reference_context_compatibility_policy",
            error_type=ContractValidationError,
        )
        object.__setattr__(
            self,
            "reference_context_compatibility_policy",
            reference_context_compatibility_policy,
        )
        if not isinstance(self.attrition_policy, KinaseAttritionPolicy):
            raise ContractValidationError(
                "scoring_config.attrition_policy must be KinaseAttritionPolicy"
            )
        _require_int_at_least(
            self.min_substrates,
            field_name="scoring_config.min_substrates",
            minimum=KINASE_SCORING_MIN_SUBSTRATES_FLOOR,
            error_type=ContractValidationError,
        )
        if not isinstance(self.localisation_requirement, LocalisationRequirement):
            raise ContractValidationError(
                "scoring_config.localisation_requirement must be "
                "LocalisationRequirement"
            )
        requested_profile = _parse_requested_reliability_profile(
            reliability_profile,
            field_name="scoring_config.reliability_profile",
        )
        effective_profile = _resolve_reliability_profile(
            config=self,
            requested_profile=requested_profile,
        )
        object.__setattr__(
            self,
            "requested_reliability_profile",
            requested_profile,
        )
        object.__setattr__(self, "reliability_profile", effective_profile)

    @property
    def effective_reliability_profile(self) -> KinaseReliabilityProfile:
        """Return the resolved reliability profile for this scoring config."""

        return self.reliability_profile

    @classmethod
    def default(cls) -> KinaseScoringConfig:
        """Deprecated alias for :meth:`exploratory`."""
        warnings.warn(
            (
                "KinaseScoringConfig.default() is deprecated because the name is "
                "ambiguous; use KinaseScoringConfig.exploratory() for the "
                "historical exploratory behavior."
            ),
            DeprecationWarning,
            stacklevel=2,
        )
        return cls.exploratory()

    @classmethod
    def exploratory(cls) -> KinaseScoringConfig:
        """Return the explicit exploratory profile preserving old defaults."""
        return cls(reliability_profile=KINASE_RELIABILITY_PROFILE_EXPLORATORY)

    @classmethod
    def strict_missing_values(cls) -> KinaseScoringConfig:
        """Return strict missing-value handling for profile aggregation."""
        return cls(
            profile_missing_value_strategy=(
                KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT
            )
        )

    @classmethod
    def production(
        cls,
        *,
        minimum_reference_overlap_fraction: float,
        minimum_sequence_supported_fraction: float,
        minimum_scored_fraction: float,
    ) -> KinaseScoringConfig:
        """Return production scoring config with caller-chosen coverage thresholds."""
        return cls(
            min_substrates=5,
            profile_self_inclusion_policy=(
                KINASE_PROFILE_SELF_INCLUSION_POLICY_LEAVE_ONE_OUT
            ),
            localisation_requirement=LocalisationRequirement.production_site_level(),
            attrition_policy=KinaseAttritionPolicy(
                minimum_reference_overlap_fraction=(minimum_reference_overlap_fraction),
                minimum_sequence_supported_fraction=(
                    minimum_sequence_supported_fraction
                ),
                minimum_scored_fraction=minimum_scored_fraction,
                on_violation=KINASE_ATTRITION_POLICY_ON_VIOLATION_ERROR,
            ),
            reliability_profile=KINASE_RELIABILITY_PROFILE_PRODUCTION,
        )


def validate_kinase_production_reliability_invariants(
    *,
    min_substrates: int,
    profile_self_inclusion_policy: ProfileSelfInclusionPolicy,
    localisation_requirement: LocalisationRequirement,
    attrition_policy: KinaseAttritionPolicy,
    field_name: str,
    error_type: type[Exception],
) -> None:
    """Validate production reliability invariants for config and workflow seams."""

    if int(min_substrates) < 5:
        raise error_type(f"{field_name}.min_substrates must be at least 5")
    if profile_self_inclusion_policy is not ProfileSelfInclusionPolicy.LEAVE_ONE_OUT:
        raise error_type(
            f"{field_name}.profile_self_inclusion_policy must be leave_one_out "
            "for production reliability"
        )
    requirement = localisation_requirement
    production_requirement = LocalisationRequirement.production_site_level()
    if not bool(requirement.require_present):
        raise error_type(
            f"{field_name}.localisation_requirement must reject unknown "
            "localisation for production reliability"
        )
    if requirement.minimum_probability is None:
        raise error_type(
            f"{field_name}.localisation_requirement.minimum_probability is "
            "required for production reliability"
        )
    if float(requirement.minimum_probability) < float(
        cast(float, production_requirement.minimum_probability)
    ):
        raise error_type(
            f"{field_name}.localisation_requirement.minimum_probability must be "
            f"at least {production_requirement.minimum_probability} for "
            "production reliability"
        )
    if attrition_policy.on_violation != KINASE_ATTRITION_POLICY_ON_VIOLATION_ERROR:
        raise error_type(
            f"{field_name}.attrition_policy.on_violation must be error for "
            "production reliability"
        )
    for threshold_name, threshold_value in _attrition_threshold_items(attrition_policy):
        if float(threshold_value) <= 0.0:
            raise error_type(
                f"{field_name}.attrition_policy.{threshold_name} must be set "
                "above 0.0 for production reliability"
            )


def _parse_requested_reliability_profile(
    value: KinaseReliabilityProfile | str | None | object,
    *,
    field_name: str,
) -> KinaseReliabilityProfile | None:
    if value is _UNSET or value is None:
        return None
    return coerce_policy_enum(
        KinaseReliabilityProfile,
        value,
        field_name=field_name,
        error_type=ContractValidationError,
    )


def _resolve_reliability_profile(
    *,
    config: KinaseScoringConfig,
    requested_profile: KinaseReliabilityProfile | None,
) -> KinaseReliabilityProfile:
    if requested_profile is None:
        if _matches_exploratory_scoring_preset(config):
            return KINASE_RELIABILITY_PROFILE_EXPLORATORY
        return KINASE_RELIABILITY_PROFILE_CUSTOM
    if requested_profile is KINASE_RELIABILITY_PROFILE_CUSTOM:
        return KINASE_RELIABILITY_PROFILE_CUSTOM
    if requested_profile is KINASE_RELIABILITY_PROFILE_EXPLORATORY:
        if not _matches_exploratory_scoring_preset(config):
            raise ContractValidationError(
                "scoring_config.reliability_profile='exploratory' requires the "
                "exploratory preset values; use reliability_profile='custom' for "
                "modified exploratory settings"
            )
        return KINASE_RELIABILITY_PROFILE_EXPLORATORY
    if requested_profile is KINASE_RELIABILITY_PROFILE_PRODUCTION:
        validate_kinase_production_reliability_invariants(
            min_substrates=int(config.min_substrates),
            profile_self_inclusion_policy=config.profile_self_inclusion_policy,
            localisation_requirement=config.localisation_requirement,
            attrition_policy=config.attrition_policy,
            field_name="scoring_config",
            error_type=ContractValidationError,
        )
        return KINASE_RELIABILITY_PROFILE_PRODUCTION
    allowed = ", ".join(str(profile) for profile in KINASE_RELIABILITY_PROFILES)
    raise ContractValidationError(
        f"scoring_config.reliability_profile must be one of: {allowed}"
    )


def _matches_exploratory_scoring_preset(config: KinaseScoringConfig) -> bool:
    return (
        int(config.min_substrates) == KINASE_SCORING_MIN_SUBSTRATES_FLOOR
        and str(config.scoring_mode) == KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED
        and bool(config.include_diagnostic_scoring_tables) is False
        and bool(config.include_substrate_contributions) is False
        and config.profile_missing_value_strategy
        == KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT
        and config.profile_self_inclusion_policy is ProfileSelfInclusionPolicy.ALLOW
        and config.attrition_policy == KinaseAttritionPolicy()
        and config.localisation_requirement == LocalisationRequirement()
        and config.reference_context_compatibility_policy
        is ReferenceContextCompatibilityPolicy.REQUIRE_KNOWN_MATCH
        and bool(config.allow_mixed_total_protein_quantitative_meaning) is False
    )


def _attrition_threshold_items(
    policy: KinaseAttritionPolicy,
) -> tuple[tuple[str, float], ...]:
    return (
        (
            "minimum_reference_overlap_fraction",
            policy.minimum_reference_overlap_fraction,
        ),
        (
            "minimum_sequence_supported_fraction",
            policy.minimum_sequence_supported_fraction,
        ),
        ("minimum_scored_fraction", policy.minimum_scored_fraction),
    )


def _require_attrition_fraction(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractValidationError(f"{field_name} must be a float or int")
    fraction = float(value)
    if not math.isfinite(fraction):
        raise ContractValidationError(f"{field_name} must be finite")
    if not 0.0 <= fraction <= 1.0:
        raise ContractValidationError(f"{field_name} must be between 0.0 and 1.0")
    return fraction


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
    kinase-condition pair has too few usable substrates for the selected method,
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
        if not isinstance(self.enabled, bool):
            raise ContractValidationError("activity_config.enabled must be a bool")
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
        if not isinstance(self.ksea_adjust_p_values, bool):
            raise ContractValidationError(
                "activity_config.ksea_adjust_p_values must be a bool"
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
        if not isinstance(self.ssgsea_adjust_p_values, bool):
            raise ContractValidationError(
                "activity_config.ssgsea_adjust_p_values must be a bool"
            )

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
    "KINASE_ATTRITION_POLICY_ON_VIOLATION_ERROR",
    "KINASE_ATTRITION_POLICY_ON_VIOLATION_MODES",
    "KINASE_ATTRITION_POLICY_ON_VIOLATION_WARN",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGIES",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGY_MEDIAN_SKIPNA",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT",
    "KINASE_PROFILE_SELF_INCLUSION_POLICIES",
    "KINASE_PROFILE_SELF_INCLUSION_POLICY_ALLOW",
    "KINASE_PROFILE_SELF_INCLUSION_POLICY_LEAVE_ONE_OUT",
    "KINASE_LIBRARY_MOTIF_ALIAS_DEPRECATION_MESSAGE",
    "KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICIES",
    "KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS",
    "KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR",
    "KINASE_RELIABILITY_PROFILE_CUSTOM",
    "KINASE_RELIABILITY_PROFILE_EXPLORATORY",
    "KINASE_RELIABILITY_PROFILE_PRODUCTION",
    "KINASE_RELIABILITY_PROFILES",
    "KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF",
    "KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF",
    "KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF",
    "KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY",
    "KINASE_SCORING_MODE_ALIASES",
    "KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED",
    "KINASE_SCORING_MODES",
    "KINASE_SCORING_MODES_REQUIRING_KINASE_LIBRARY",
    "KINASE_SCORING_MIN_SUBSTRATES_FLOOR",
    "KinaseActivityConfig",
    "KinaseActivityMethod",
    "KinaseActivityPValueMethod",
    "KinaseActivitySsgseaRankingDirection",
    "KinaseAttritionPolicy",
    "KinaseAttritionViolationMode",
    "KinaseProfileMissingValueStrategy",
    "KinaseReferenceDisplayAmbiguityPolicy",
    "KinaseReliabilityProfile",
    "KinaseScoringMode",
    "KinaseSiteSequenceConflictPolicy",
    "LocalisationRequirement",
    "ProfileSelfInclusionPolicy",
    "ReferenceContextCompatibilityPolicy",
    "KinaseScoringConfig",
    "normalize_kinase_scoring_mode",
]
