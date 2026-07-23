"""Public signalome workflow configuration models."""
# pyright: reportUnnecessaryIsInstance=false
# Runtime boundary guards are intentionally retained for untyped external callers.

from __future__ import annotations

from dataclasses import dataclass, field

from phospy.contracts.configs._validation import coerce_policy_enum
from phospy.contracts.configs.common import _require_int_at_least, _require_real_between
from phospy.contracts.configs.localisation import LocalisationRequirement
from phospy.contracts.configs.reference_context import (
    REFERENCE_CONTEXT_COMPATIBILITY_POLICY_REQUIRE_KNOWN_MATCH,
    ReferenceContextCompatibilityPolicy,
)
from phospy.errors.validation import ContractValidationError
from phospy.science.configs.signalome import (
    SIGNALOME_ALLOW_MIXED_TOTAL_PROTEIN_QUANTITATIVE_MEANING_DEFAULT,
    SIGNALOME_ASSIGNMENT_POLICIES,
    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
    SIGNALOME_CANDIDATE_SCORING_POLICIES,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
    SIGNALOME_CLUSTERING_ENGINES,
    SIGNALOME_KINASE_NETWORK_POLICIES,
    SIGNALOME_KINASE_NETWORK_POLICY_ABSOLUTE_THRESHOLD,
    SIGNALOME_KINASE_NETWORK_POLICY_POSITIVE_ONLY,
    SIGNALOME_KINASE_NETWORK_POLICY_SIGNED,
    SIGNALOME_MAX_EXACT_TREE_SITES_DEFAULT,
    SIGNALOME_MAX_EXACT_TREE_SITES_FLOOR,
    SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_DEFAULT,
    SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_FLOOR,
    SIGNALOME_MODULE_COUNT_FLOOR,
    SIGNALOME_MODULE_SELECTION_FALLBACK_THRESHOLD_DEFAULT,
    SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_DEFAULT,
    SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_FLOOR,
    SIGNALOME_MODULE_SELECTION_PRIMARY_THRESHOLD_DEFAULT,
    SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT,
    SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_FLOOR,
    SIGNALOME_SCORE_PRECONDITIONING_POLICIES,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
    SignalomeAssignmentPolicy,
    SignalomeCandidateScoringPolicy,
    SignalomeClusteringEngine,
    SignalomeKinaseNetworkPolicy,
    SignalomeScorePreconditioningPolicy,
)


@dataclass(frozen=True, slots=True)
class SignalomeScientificConfig:
    """Scientific interpretation choices for score-derived signalome summaries."""

    substrate_support_cutoff: float = 0.5
    assignment_policy: SignalomeAssignmentPolicy = (
        SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY
    )

    def __post_init__(self) -> None:
        _require_real_between(
            self.substrate_support_cutoff,
            field_name=(
                "signalome workflow request config.scientific.substrate_support_cutoff"
            ),
            minimum=0.0,
            maximum=1.0,
            error_type=ContractValidationError,
        )
        if self.assignment_policy not in SIGNALOME_ASSIGNMENT_POLICIES:
            allowed_policies = ", ".join(sorted(SIGNALOME_ASSIGNMENT_POLICIES))
            raise ContractValidationError(
                "signalome workflow request config.scientific.assignment_policy "
                f"must be one of: {allowed_policies}"
            )


@dataclass(frozen=True, slots=True)
class SignalomeClusteringConfig:
    """Clustering and module-selection behaviour for the signalome workflow."""

    module_count: int | None = None
    module_selection_primary_correlation_threshold: float = (
        SIGNALOME_MODULE_SELECTION_PRIMARY_THRESHOLD_DEFAULT
    )
    module_selection_fallback_correlation_threshold: float = (
        SIGNALOME_MODULE_SELECTION_FALLBACK_THRESHOLD_DEFAULT
    )
    module_selection_max_clusters: int = SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_DEFAULT
    candidate_scoring_policy: SignalomeCandidateScoringPolicy = (
        SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    )
    clustering_engine: SignalomeClusteringEngine = (
        SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL
    )

    def __post_init__(self) -> None:
        if self.module_count is not None:
            _require_int_at_least(
                self.module_count,
                field_name=(
                    "signalome workflow request config.clustering.module_count"
                ),
                minimum=SIGNALOME_MODULE_COUNT_FLOOR,
                error_type=ContractValidationError,
            )
        _require_real_between(
            self.module_selection_primary_correlation_threshold,
            field_name=(
                "signalome workflow request config.clustering."
                "module_selection_primary_correlation_threshold"
            ),
            minimum=0.0,
            maximum=1.0,
            error_type=ContractValidationError,
        )
        _require_real_between(
            self.module_selection_fallback_correlation_threshold,
            field_name=(
                "signalome workflow request config.clustering."
                "module_selection_fallback_correlation_threshold"
            ),
            minimum=0.0,
            maximum=1.0,
            error_type=ContractValidationError,
        )
        _require_int_at_least(
            self.module_selection_max_clusters,
            field_name=(
                "signalome workflow request config.clustering."
                "module_selection_max_clusters"
            ),
            minimum=SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_FLOOR,
            error_type=ContractValidationError,
        )
        if self.candidate_scoring_policy not in SIGNALOME_CANDIDATE_SCORING_POLICIES:
            allowed = ", ".join(sorted(SIGNALOME_CANDIDATE_SCORING_POLICIES))
            raise ContractValidationError(
                "signalome workflow request config.clustering."
                "candidate_scoring_policy "
                f"must be one of: {allowed}"
            )
        if self.clustering_engine not in SIGNALOME_CLUSTERING_ENGINES:
            allowed = ", ".join(sorted(SIGNALOME_CLUSTERING_ENGINES))
            raise ContractValidationError(
                "signalome workflow request config.clustering.clustering_engine "
                f"must be one of: {allowed}"
            )


@dataclass(frozen=True, slots=True)
class SignalomeValidationConfig:
    """Validation and strictness policy for signalome inputs."""

    score_preconditioning_policy: SignalomeScorePreconditioningPolicy = (
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP
    )
    localisation_requirement: LocalisationRequirement = field(
        default_factory=LocalisationRequirement
    )
    allow_mixed_total_protein_quantitative_meaning: bool = (
        SIGNALOME_ALLOW_MIXED_TOTAL_PROTEIN_QUANTITATIVE_MEANING_DEFAULT
    )
    reference_context_compatibility_policy: ReferenceContextCompatibilityPolicy = (
        REFERENCE_CONTEXT_COMPATIBILITY_POLICY_REQUIRE_KNOWN_MATCH
    )

    def __post_init__(self) -> None:
        if not isinstance(self.allow_mixed_total_protein_quantitative_meaning, bool):
            raise ContractValidationError(
                "signalome workflow request config.validation."
                "allow_mixed_total_protein_quantitative_meaning must be a bool"
            )
        if (
            self.score_preconditioning_policy
            not in SIGNALOME_SCORE_PRECONDITIONING_POLICIES
        ):
            allowed_policies = ", ".join(
                sorted(SIGNALOME_SCORE_PRECONDITIONING_POLICIES)
            )
            raise ContractValidationError(
                "signalome workflow request config.validation."
                "score_preconditioning_policy "
                f"must be one of: {allowed_policies}"
            )
        if not isinstance(self.localisation_requirement, LocalisationRequirement):
            raise ContractValidationError(
                "signalome workflow request config.validation."
                "localisation_requirement must be LocalisationRequirement"
            )
        reference_context_compatibility_policy = coerce_policy_enum(
            ReferenceContextCompatibilityPolicy,
            self.reference_context_compatibility_policy,
            field_name=(
                "signalome workflow request config.validation."
                "reference_context_compatibility_policy"
            ),
            error_type=ContractValidationError,
        )
        object.__setattr__(
            self,
            "reference_context_compatibility_policy",
            reference_context_compatibility_policy,
        )


@dataclass(frozen=True, slots=True)
class SignalomeOutputConfig:
    """Output-shape and score-profile association edge settings."""

    network_correlation_threshold: float = 0.5
    network_policy: SignalomeKinaseNetworkPolicy = (
        SIGNALOME_KINASE_NETWORK_POLICY_SIGNED
    )
    network_min_paired_finite_observations: int | None = None

    def __post_init__(self) -> None:
        _require_real_between(
            self.network_correlation_threshold,
            field_name=(
                "signalome workflow request config.output.network_correlation_threshold"
            ),
            minimum=0.0,
            maximum=1.0,
            error_type=ContractValidationError,
        )
        if self.network_policy not in SIGNALOME_KINASE_NETWORK_POLICIES:
            allowed_policies = ", ".join(sorted(SIGNALOME_KINASE_NETWORK_POLICIES))
            raise ContractValidationError(
                "signalome workflow request config.output.network_policy "
                f"must be one of: {allowed_policies}"
            )
        if self.network_min_paired_finite_observations is not None:
            _require_int_at_least(
                self.network_min_paired_finite_observations,
                field_name=(
                    "signalome workflow request config.output."
                    "network_min_paired_finite_observations"
                ),
                minimum=SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_FLOOR,
                error_type=ContractValidationError,
            )


@dataclass(frozen=True, slots=True)
class SignalomePerformanceConfig:
    """Advanced clustering-scale guardrails for large signalome runs."""

    max_exact_tree_sites: int = SIGNALOME_MAX_EXACT_TREE_SITES_DEFAULT
    max_full_candidate_scoring_sites: int = (
        SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_DEFAULT
    )

    def __post_init__(self) -> None:
        _require_int_at_least(
            self.max_exact_tree_sites,
            field_name=(
                "signalome workflow request config.performance.max_exact_tree_sites"
            ),
            minimum=SIGNALOME_MAX_EXACT_TREE_SITES_FLOOR,
            error_type=ContractValidationError,
        )
        _require_int_at_least(
            self.max_full_candidate_scoring_sites,
            field_name=(
                "signalome workflow request config.performance."
                "max_full_candidate_scoring_sites"
            ),
            minimum=SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_FLOOR,
            error_type=ContractValidationError,
        )


@dataclass(frozen=True, slots=True)
class SignalomeConfig:
    """Public signalome workflow configuration grouped by user intent."""

    scientific: SignalomeScientificConfig = field(
        default_factory=SignalomeScientificConfig
    )
    clustering: SignalomeClusteringConfig = field(
        default_factory=SignalomeClusteringConfig
    )
    validation: SignalomeValidationConfig = field(
        default_factory=SignalomeValidationConfig
    )
    output: SignalomeOutputConfig = field(default_factory=SignalomeOutputConfig)
    performance: SignalomePerformanceConfig = field(
        default_factory=SignalomePerformanceConfig
    )

    def __post_init__(self) -> None:
        if not isinstance(self.scientific, SignalomeScientificConfig):
            raise ContractValidationError(
                "signalome workflow request config.scientific must be "
                "SignalomeScientificConfig"
            )
        if not isinstance(self.clustering, SignalomeClusteringConfig):
            raise ContractValidationError(
                "signalome workflow request config.clustering must be "
                "SignalomeClusteringConfig"
            )
        if not isinstance(self.validation, SignalomeValidationConfig):
            raise ContractValidationError(
                "signalome workflow request config.validation must be "
                "SignalomeValidationConfig"
            )
        if not isinstance(self.output, SignalomeOutputConfig):
            raise ContractValidationError(
                "signalome workflow request config.output must be SignalomeOutputConfig"
            )
        if not isinstance(self.performance, SignalomePerformanceConfig):
            raise ContractValidationError(
                "signalome workflow request config.performance must be "
                "SignalomePerformanceConfig"
            )

    @classmethod
    def strict(cls) -> SignalomeConfig:
        """Return strict score preconditioning for signalome workflows."""
        return cls(
            validation=SignalomeValidationConfig(
                score_preconditioning_policy=(
                    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP
                )
            )
        )

    @classmethod
    def permissive_missing_scores(cls) -> SignalomeConfig:
        """Allow and report all-missing score-row drops."""
        return cls(
            validation=SignalomeValidationConfig(
                score_preconditioning_policy=(
                    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
                )
            )
        )

    @classmethod
    def sampled_candidate_scoring(cls) -> SignalomeConfig:
        """Return config using sampled candidate module-count scoring while keeping scale guards enabled."""
        return cls(
            clustering=SignalomeClusteringConfig(
                candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
            )
        )

    @classmethod
    def production(cls) -> SignalomeConfig:
        """Return production signalome config with strict site-level localisation."""
        return cls(
            validation=SignalomeValidationConfig(
                localisation_requirement=LocalisationRequirement.production_site_level()
            )
        )


__all__ = [
    "SIGNALOME_ALLOW_MIXED_TOTAL_PROTEIN_QUANTITATIVE_MEANING_DEFAULT",
    "SIGNALOME_ASSIGNMENT_POLICIES",
    "SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY",
    "SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP",
    "SIGNALOME_CANDIDATE_SCORING_POLICIES",
    "SIGNALOME_CANDIDATE_SCORING_POLICY_FULL",
    "SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED",
    "SIGNALOME_CLUSTERING_ENGINES",
    "SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON",
    "SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL",
    "SIGNALOME_KINASE_NETWORK_POLICIES",
    "SIGNALOME_KINASE_NETWORK_POLICY_ABSOLUTE_THRESHOLD",
    "SIGNALOME_KINASE_NETWORK_POLICY_POSITIVE_ONLY",
    "SIGNALOME_KINASE_NETWORK_POLICY_SIGNED",
    "SIGNALOME_MAX_EXACT_TREE_SITES_DEFAULT",
    "SIGNALOME_MAX_EXACT_TREE_SITES_FLOOR",
    "SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_DEFAULT",
    "SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_FLOOR",
    "SIGNALOME_MODULE_COUNT_FLOOR",
    "SIGNALOME_MODULE_SELECTION_FALLBACK_THRESHOLD_DEFAULT",
    "SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_DEFAULT",
    "SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_FLOOR",
    "SIGNALOME_MODULE_SELECTION_PRIMARY_THRESHOLD_DEFAULT",
    "SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_DEFAULT",
    "SIGNALOME_NETWORK_MIN_PAIRED_FINITE_OBSERVATIONS_FLOOR",
    "SIGNALOME_SCORE_PRECONDITIONING_POLICIES",
    "SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT",
    "SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP",
    "SignalomeAssignmentPolicy",
    "SignalomeCandidateScoringPolicy",
    "SignalomeClusteringConfig",
    "SignalomeClusteringEngine",
    "SignalomeConfig",
    "SignalomeKinaseNetworkPolicy",
    "LocalisationRequirement",
    "SignalomeOutputConfig",
    "SignalomePerformanceConfig",
    "SignalomeScientificConfig",
    "SignalomeScorePreconditioningPolicy",
    "SignalomeValidationConfig",
    "ReferenceContextCompatibilityPolicy",
]
