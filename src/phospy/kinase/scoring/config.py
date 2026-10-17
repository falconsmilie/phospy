"""Public kinase scoring configuration models."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal, cast

from phospy.errors.validation import ContractValidationError
from phospy.evidence.config import LocalisationRequirement
from phospy.kinase.scoring.policies import ProfileSelfInclusionPolicy
from phospy.policies import PolicyEnum
from phospy.policies.config_values import (
    coerce_policy_enum,
)
from phospy.policies.config_values import (
    require_int_at_least as _require_int_at_least,
)
from phospy.references.config import (
    REFERENCE_CONTEXT_COMPATIBILITY_POLICY_REQUIRE_KNOWN_MATCH,
    ReferenceContextCompatibilityPolicy,
)

KINASE_SCORING_MIN_SUBSTRATES_FLOOR = 2
KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED = "phosr_rank_weighted"
KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY = "kinase_library_motif_only"
KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF = "kinase_library_contextual_motif"
KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF = "combined_profile_motif"
KinaseScoringMode = Literal[
    "phosr_rank_weighted",
    "kinase_library_contextual_motif",
    "kinase_library_motif_only",
    "combined_profile_motif",
]
KINASE_SCORING_MODES = frozenset(
    {
        KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
        KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
        KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
        KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
    }
)
KINASE_SCORING_MODES_REQUIRING_KINASE_LIBRARY = frozenset(
    {
        KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
        KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
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
KINASE_PROFILE_SELF_INCLUSION_POLICY_ALLOW = ProfileSelfInclusionPolicy.ALLOW
KINASE_PROFILE_SELF_INCLUSION_POLICY_LEAVE_ONE_OUT = (
    ProfileSelfInclusionPolicy.LEAVE_ONE_OUT
)
KINASE_PROFILE_SELF_INCLUSION_POLICIES = frozenset(
    {
        KINASE_PROFILE_SELF_INCLUSION_POLICY_ALLOW,
        KINASE_PROFILE_SELF_INCLUSION_POLICY_LEAVE_ONE_OUT,
    }
)
KINASE_ALLOW_MIXED_TOTAL_PROTEIN_QUANTITATIVE_MEANING_DEFAULT = False


class KinaseReliabilityProfile(PolicyEnum):
    EXPLORATORY = "exploratory"
    PRODUCTION = "production"
    CUSTOM = "custom"


KINASE_RELIABILITY_PROFILE_EXPLORATORY = KinaseReliabilityProfile.EXPLORATORY
KINASE_RELIABILITY_PROFILE_PRODUCTION = KinaseReliabilityProfile.PRODUCTION
KINASE_RELIABILITY_PROFILE_CUSTOM = KinaseReliabilityProfile.CUSTOM
KINASE_RELIABILITY_PROFILES = frozenset(
    {
        KINASE_RELIABILITY_PROFILE_EXPLORATORY,
        KINASE_RELIABILITY_PROFILE_PRODUCTION,
        KINASE_RELIABILITY_PROFILE_CUSTOM,
    }
)


class KinaseSiteSequenceConflictPolicy(PolicyEnum):
    ERROR = "error"
    PREFER_REFERENCE = "prefer_reference"
    PREFER_DATASET = "prefer_dataset"


KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR = KinaseSiteSequenceConflictPolicy.ERROR
KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE = (
    KinaseSiteSequenceConflictPolicy.PREFER_REFERENCE
)
KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET = (
    KinaseSiteSequenceConflictPolicy.PREFER_DATASET
)
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
KINASE_ATTRITION_POLICY_ON_VIOLATION_ERROR = "error"
KINASE_ATTRITION_POLICY_ON_VIOLATION_WARN = "warn"
KinaseAttritionViolationMode = Literal["error", "warn"]
KINASE_ATTRITION_POLICY_ON_VIOLATION_MODES = frozenset(
    {
        KINASE_ATTRITION_POLICY_ON_VIOLATION_ERROR,
        KINASE_ATTRITION_POLICY_ON_VIOLATION_WARN,
    }
)


def normalize_kinase_scoring_mode(value: object) -> KinaseScoringMode:
    """Return a typed scoring-mode value without accepting legacy aliases."""

    return cast(KinaseScoringMode, str(value))


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

    ``reliability_profile`` is required for direct construction and exposes
    whether scoring is using the explicit exploratory preset, the stricter
    production preset, or custom caller values. Use
    ``KinaseScoringConfig.exploratory()`` for the historical permissive preset,
    ``KinaseScoringConfig.production(...)`` for production reliability with
    caller-supplied coverage thresholds, or direct construction with
    ``reliability_profile="custom"`` for caller-defined settings.

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
    not imply official Kinase Library predictor parity.

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
        include_diagnostic_scoring_tables = _require_bool(
            self.include_diagnostic_scoring_tables,
            field_name="scoring_config.include_diagnostic_scoring_tables",
        )
        include_substrate_contributions = _require_bool(
            self.include_substrate_contributions,
            field_name="scoring_config.include_substrate_contributions",
        )
        allow_mixed_total_protein_quantitative_meaning = _require_bool(
            self.allow_mixed_total_protein_quantitative_meaning,
            field_name=(
                "scoring_config.allow_mixed_total_protein_quantitative_meaning"
            ),
        )
        normalized_scoring_mode = normalize_kinase_scoring_mode(self.scoring_mode)
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
        attrition_policy = _require_kinase_attrition_policy(
            self.attrition_policy,
            field_name="scoring_config.attrition_policy",
        )
        localisation_requirement = _require_localisation_requirement(
            self.localisation_requirement,
            field_name="scoring_config.localisation_requirement",
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
            "allow_mixed_total_protein_quantitative_meaning",
            allow_mixed_total_protein_quantitative_meaning,
        )
        object.__setattr__(self, "attrition_policy", attrition_policy)
        object.__setattr__(self, "localisation_requirement", localisation_requirement)
        _require_int_at_least(
            self.min_substrates,
            field_name="scoring_config.min_substrates",
            minimum=KINASE_SCORING_MIN_SUBSTRATES_FLOOR,
            error_type=ContractValidationError,
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
    def exploratory(cls) -> KinaseScoringConfig:
        """Return the explicit exploratory profile preserving old defaults."""
        return cls(reliability_profile=KINASE_RELIABILITY_PROFILE_EXPLORATORY)

    @classmethod
    def strict_missing_values(cls) -> KinaseScoringConfig:
        """Return strict missing-value handling for profile aggregation."""
        return cls(
            profile_missing_value_strategy=(
                KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT
            ),
            reliability_profile=KINASE_RELIABILITY_PROFILE_EXPLORATORY,
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

    @classmethod
    def _from_legacy_inferred_profile(
        cls,
        *,
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
    ) -> KinaseScoringConfig:
        """Restore historical payloads that predate explicit profile selection."""

        custom_config = cls(
            min_substrates=min_substrates,
            scoring_mode=scoring_mode,
            include_diagnostic_scoring_tables=include_diagnostic_scoring_tables,
            include_substrate_contributions=include_substrate_contributions,
            profile_missing_value_strategy=profile_missing_value_strategy,
            profile_self_inclusion_policy=profile_self_inclusion_policy,
            attrition_policy=attrition_policy,
            localisation_requirement=localisation_requirement,
            reference_context_compatibility_policy=(
                reference_context_compatibility_policy
            ),
            allow_mixed_total_protein_quantitative_meaning=(
                allow_mixed_total_protein_quantitative_meaning
            ),
            reliability_profile=KINASE_RELIABILITY_PROFILE_CUSTOM,
        )
        if not _matches_exploratory_scoring_preset(custom_config):
            return custom_config
        return cls(
            min_substrates=min_substrates,
            scoring_mode=scoring_mode,
            include_diagnostic_scoring_tables=include_diagnostic_scoring_tables,
            include_substrate_contributions=include_substrate_contributions,
            profile_missing_value_strategy=profile_missing_value_strategy,
            profile_self_inclusion_policy=profile_self_inclusion_policy,
            attrition_policy=attrition_policy,
            localisation_requirement=localisation_requirement,
            reference_context_compatibility_policy=(
                reference_context_compatibility_policy
            ),
            allow_mixed_total_protein_quantitative_meaning=(
                allow_mixed_total_protein_quantitative_meaning
            ),
            reliability_profile=KINASE_RELIABILITY_PROFILE_EXPLORATORY,
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
) -> KinaseReliabilityProfile:
    if value is _UNSET or value is None:
        raise ContractValidationError(
            "scoring_config.reliability_profile is required; choose "
            "KinaseScoringConfig.exploratory(), "
            "KinaseScoringConfig.production(...), or direct construction with "
            "reliability_profile='custom'"
        )
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


def _require_bool(value: object, *, field_name: str) -> bool:
    if not isinstance(value, bool):
        raise ContractValidationError(f"{field_name} must be a bool")
    return value


def _require_kinase_attrition_policy(
    value: object,
    *,
    field_name: str,
) -> KinaseAttritionPolicy:
    if not isinstance(value, KinaseAttritionPolicy):
        raise ContractValidationError(f"{field_name} must be KinaseAttritionPolicy")
    return value


def _require_localisation_requirement(
    value: object,
    *,
    field_name: str,
) -> LocalisationRequirement:
    if not isinstance(value, LocalisationRequirement):
        raise ContractValidationError(f"{field_name} must be LocalisationRequirement")
    return value


__all__ = [
    "KINASE_ALLOW_MIXED_TOTAL_PROTEIN_QUANTITATIVE_MEANING_DEFAULT",
    "KINASE_ATTRITION_POLICY_ON_VIOLATION_ERROR",
    "KINASE_ATTRITION_POLICY_ON_VIOLATION_MODES",
    "KINASE_ATTRITION_POLICY_ON_VIOLATION_WARN",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGIES",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGY_MEDIAN_SKIPNA",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT",
    "KINASE_PROFILE_SELF_INCLUSION_POLICIES",
    "KINASE_PROFILE_SELF_INCLUSION_POLICY_ALLOW",
    "KINASE_PROFILE_SELF_INCLUSION_POLICY_LEAVE_ONE_OUT",
    "KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICIES",
    "KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS",
    "KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR",
    "KINASE_RELIABILITY_PROFILE_CUSTOM",
    "KINASE_RELIABILITY_PROFILE_EXPLORATORY",
    "KINASE_RELIABILITY_PROFILE_PRODUCTION",
    "KINASE_RELIABILITY_PROFILES",
    "KINASE_SCORING_MIN_SUBSTRATES_FLOOR",
    "KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF",
    "KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF",
    "KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY",
    "KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED",
    "KINASE_SCORING_MODES",
    "KINASE_SCORING_MODES_REQUIRING_KINASE_LIBRARY",
    "KINASE_SITE_SEQUENCE_CONFLICT_POLICIES",
    "KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR",
    "KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET",
    "KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE",
    "KinaseAttritionPolicy",
    "KinaseAttritionViolationMode",
    "KinaseProfileMissingValueStrategy",
    "KinaseReferenceDisplayAmbiguityPolicy",
    "KinaseReliabilityProfile",
    "KinaseScoringConfig",
    "KinaseScoringMode",
    "KinaseSiteSequenceConflictPolicy",
    "LocalisationRequirement",
    "ProfileSelfInclusionPolicy",
    "ReferenceContextCompatibilityPolicy",
    "normalize_kinase_scoring_mode",
    "validate_kinase_production_reliability_invariants",
]
