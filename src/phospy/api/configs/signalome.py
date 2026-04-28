"""Public signalome workflow configuration models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from phospy.api.configs.common import _require_int_at_least, _require_real_between
from phospy.errors.validation import WorkflowValidationError

SIGNALOME_MODULE_COUNT_FLOOR = 1
SIGNALOME_MODULE_SELECTION_PRIMARY_THRESHOLD_DEFAULT = 0.5
SIGNALOME_MODULE_SELECTION_FALLBACK_THRESHOLD_DEFAULT = 0.1
SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_FLOOR = 1
SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_DEFAULT = 10

SIGNALOME_TREE_ENGINE_EXACT = "exact"
SignalomeTreeEngine = Literal["exact"]
SIGNALOME_TREE_ENGINES = frozenset({SIGNALOME_TREE_ENGINE_EXACT})

SIGNALOME_CANDIDATE_SCORING_POLICY_FULL = "full"
SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED = "sampled"
SignalomeCandidateScoringPolicy = Literal["full", "sampled"]
SIGNALOME_CANDIDATE_SCORING_POLICIES = frozenset(
    {
        SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    }
)

SIGNALOME_MAX_EXACT_TREE_SITES_FLOOR = 1
SIGNALOME_MAX_EXACT_TREE_SITES_DEFAULT = 2000
SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_FLOOR = 1
SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_DEFAULT = 2000

SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON = "exact_python"
SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL = "scipy_hierarchical"
SignalomeClusteringEngine = Literal["exact_python", "scipy_hierarchical"]
SIGNALOME_CLUSTERING_ENGINES = frozenset(
    {
        SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
        SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
    }
)

SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY = "cutoff_binary"
SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP = "weighted_top"
SignalomeAssignmentPolicy = Literal["cutoff_binary", "weighted_top"]
SIGNALOME_ASSIGNMENT_POLICIES = frozenset(
    {
        SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
        SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
    }
)

SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT = "allow_and_report"
SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP = "error_on_drop"
SignalomeScorePreconditioningPolicy = Literal["allow_and_report", "error_on_drop"]
SIGNALOME_SCORE_PRECONDITIONING_POLICIES = frozenset(
    {
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
    }
)

SIGNALOME_KINASE_NETWORK_POLICY_POSITIVE_ONLY = "positive_only"
SIGNALOME_KINASE_NETWORK_POLICY_ABSOLUTE_THRESHOLD = "absolute_threshold"
SIGNALOME_KINASE_NETWORK_POLICY_SIGNED = "signed"
SignalomeKinaseNetworkPolicy = Literal[
    "positive_only",
    "absolute_threshold",
    "signed",
]
SIGNALOME_KINASE_NETWORK_POLICIES = frozenset(
    {
        SIGNALOME_KINASE_NETWORK_POLICY_POSITIVE_ONLY,
        SIGNALOME_KINASE_NETWORK_POLICY_ABSOLUTE_THRESHOLD,
        SIGNALOME_KINASE_NETWORK_POLICY_SIGNED,
    }
)


@dataclass(frozen=True, slots=True)
class SignalomeConfig:
    """Public signalome workflow configuration."""

    substrate_support_cutoff: float = 0.5
    network_correlation_threshold: float = 0.5
    network_policy: SignalomeKinaseNetworkPolicy = (
        SIGNALOME_KINASE_NETWORK_POLICY_SIGNED
    )
    assignment_policy: SignalomeAssignmentPolicy = (
        SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY
    )
    score_preconditioning_policy: SignalomeScorePreconditioningPolicy = (
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
    )
    module_count: int | None = None
    module_selection_primary_correlation_threshold: float = (
        SIGNALOME_MODULE_SELECTION_PRIMARY_THRESHOLD_DEFAULT
    )
    module_selection_fallback_correlation_threshold: float = (
        SIGNALOME_MODULE_SELECTION_FALLBACK_THRESHOLD_DEFAULT
    )
    module_selection_max_clusters: int = SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_DEFAULT
    tree_engine: SignalomeTreeEngine = SIGNALOME_TREE_ENGINE_EXACT
    candidate_scoring_policy: SignalomeCandidateScoringPolicy = (
        SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    )
    max_exact_tree_sites: int = SIGNALOME_MAX_EXACT_TREE_SITES_DEFAULT
    max_full_candidate_scoring_sites: int = (
        SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_DEFAULT
    )
    clustering_engine: SignalomeClusteringEngine = (
        SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON
    )

    def __post_init__(self) -> None:
        _require_real_between(
            self.substrate_support_cutoff,
            field_name="signalome workflow request config.substrate_support_cutoff",
            minimum=0.0,
            maximum=1.0,
            error_type=WorkflowValidationError,
        )
        _require_real_between(
            self.network_correlation_threshold,
            field_name=(
                "signalome workflow request config.network_correlation_threshold"
            ),
            minimum=0.0,
            maximum=1.0,
            error_type=WorkflowValidationError,
        )
        if self.network_policy not in SIGNALOME_KINASE_NETWORK_POLICIES:
            allowed_policies = ", ".join(sorted(SIGNALOME_KINASE_NETWORK_POLICIES))
            raise WorkflowValidationError(
                "signalome workflow request config.network_policy "
                f"must be one of: {allowed_policies}"
            )
        if self.assignment_policy not in SIGNALOME_ASSIGNMENT_POLICIES:
            allowed_policies = ", ".join(sorted(SIGNALOME_ASSIGNMENT_POLICIES))
            raise WorkflowValidationError(
                "signalome workflow request config.assignment_policy "
                f"must be one of: {allowed_policies}"
            )
        if (
            self.score_preconditioning_policy
            not in SIGNALOME_SCORE_PRECONDITIONING_POLICIES
        ):
            allowed_policies = ", ".join(
                sorted(SIGNALOME_SCORE_PRECONDITIONING_POLICIES)
            )
            raise WorkflowValidationError(
                "signalome workflow request config.score_preconditioning_policy "
                f"must be one of: {allowed_policies}"
            )
        if self.module_count is not None:
            _require_int_at_least(
                self.module_count,
                field_name="signalome workflow request config.module_count",
                minimum=SIGNALOME_MODULE_COUNT_FLOOR,
                error_type=WorkflowValidationError,
            )
        _require_real_between(
            self.module_selection_primary_correlation_threshold,
            field_name=(
                "signalome workflow request config."
                "module_selection_primary_correlation_threshold"
            ),
            minimum=0.0,
            maximum=1.0,
            error_type=WorkflowValidationError,
        )
        _require_real_between(
            self.module_selection_fallback_correlation_threshold,
            field_name=(
                "signalome workflow request config."
                "module_selection_fallback_correlation_threshold"
            ),
            minimum=0.0,
            maximum=1.0,
            error_type=WorkflowValidationError,
        )
        _require_int_at_least(
            self.module_selection_max_clusters,
            field_name="signalome workflow request config.module_selection_max_clusters",
            minimum=SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_FLOOR,
            error_type=WorkflowValidationError,
        )
        if self.tree_engine not in SIGNALOME_TREE_ENGINES:
            allowed = ", ".join(sorted(SIGNALOME_TREE_ENGINES))
            raise WorkflowValidationError(
                "signalome workflow request config.tree_engine "
                f"must be one of: {allowed}"
            )
        if self.candidate_scoring_policy not in SIGNALOME_CANDIDATE_SCORING_POLICIES:
            allowed = ", ".join(sorted(SIGNALOME_CANDIDATE_SCORING_POLICIES))
            raise WorkflowValidationError(
                "signalome workflow request config.candidate_scoring_policy "
                f"must be one of: {allowed}"
            )
        _require_int_at_least(
            self.max_exact_tree_sites,
            field_name="signalome workflow request config.max_exact_tree_sites",
            minimum=SIGNALOME_MAX_EXACT_TREE_SITES_FLOOR,
            error_type=WorkflowValidationError,
        )
        _require_int_at_least(
            self.max_full_candidate_scoring_sites,
            field_name="signalome workflow request config.max_full_candidate_scoring_sites",
            minimum=SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_FLOOR,
            error_type=WorkflowValidationError,
        )
        if self.clustering_engine not in SIGNALOME_CLUSTERING_ENGINES:
            allowed = ", ".join(sorted(SIGNALOME_CLUSTERING_ENGINES))
            raise WorkflowValidationError(
                "signalome workflow request config.clustering_engine "
                f"must be one of: {allowed}"
            )


__all__ = [
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
    "SIGNALOME_SCORE_PRECONDITIONING_POLICIES",
    "SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT",
    "SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP",
    "SIGNALOME_TREE_ENGINES",
    "SIGNALOME_TREE_ENGINE_EXACT",
    "SignalomeAssignmentPolicy",
    "SignalomeCandidateScoringPolicy",
    "SignalomeClusteringEngine",
    "SignalomeConfig",
    "SignalomeKinaseNetworkPolicy",
    "SignalomeScorePreconditioningPolicy",
    "SignalomeTreeEngine",
]
