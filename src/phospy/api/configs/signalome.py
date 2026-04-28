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
SIGNALOME_CLUSTER_TREE_BACKEND_EXACT = "exact"
SignalomeClusterTreeBackend = Literal["exact"]
SIGNALOME_CLUSTER_TREE_BACKENDS = frozenset(
    {
        SIGNALOME_CLUSTER_TREE_BACKEND_EXACT,
    }
)
SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL = "full"
SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED = "sampled"
SignalomeCandidateScoringBackend = Literal["full", "sampled"]
SIGNALOME_CANDIDATE_SCORING_BACKENDS = frozenset(
    {
        SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
        SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED,
    }
)
SIGNALOME_MAX_EXACT_CLUSTER_TREE_SITES_FLOOR = 1
# This default matches the established full-correlation contract threshold.
SIGNALOME_MAX_EXACT_CLUSTER_TREE_SITES_DEFAULT = 2000
SIGNALOME_MAX_FULL_CORRELATION_SITES_FLOOR = 1
SIGNALOME_MAX_FULL_CORRELATION_SITES_DEFAULT = 2000
SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON = "exact_python"
SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL = "scipy_hierarchical"
SignalomeClusteringBackend = Literal["exact_python", "scipy_hierarchical"]
SIGNALOME_CLUSTERING_BACKENDS = frozenset(
    {
        SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON,
        SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL,
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
    """Public signalome workflow configuration.

    `network_policy` controls how score correlations are thresholded and encoded
    in `kinase_network.edges.correlation`:

    - `"positive_only"`: keep only positive correlations `>= threshold`.
    - `"absolute_threshold"`: keep correlations where `abs(correlation) >=
      threshold` and emit unsigned absolute correlation values.
    - `"signed"`: keep correlations where `abs(correlation) >= threshold` and
      emit signed correlation values.

    `assignment_policy` controls module-support attribution:

    - `"cutoff_binary"`: binary support per kinase from
      `substrate_support_cutoff`.
    - `"weighted_top"`: fractional support propagated from per-site
      `top_kinase_weights` ties.

    `score_preconditioning_policy` controls how all-missing downstream score
    rows are handled before score-driven signalome stages:

    - `"allow_and_report"` (default): drop all-missing rows and continue,
      reporting exact counts in diagnostics.
    - `"error_on_drop"`: fail signalome interpretation when any all-missing
      rows would be dropped.

    `cluster_tree_backend` controls site-cluster tree construction:

    - `"exact"`: build the exact Ward/agglomerative cluster tree.
    Alias: `tree_engine`.

    `candidate_scoring_backend` controls module-selection candidate scoring:

    - `"full"`: use full site-by-site correlations for candidate scoring.
    - `"sampled"`: use sampled within-cluster correlation estimates.
      This only affects candidate module-count evaluation. It does not make
      signalome clustering approximate.
    Alias: `candidate_scoring_policy`.

    `max_exact_cluster_tree_sites` hard-limits exact cluster-tree
    construction. This guard applies regardless of candidate scoring mode.
    Final module assignment still depends on the exact cluster tree.
    `max_full_correlation_sites` hard-limits full candidate-correlation scoring.

    `clustering_backend` controls which internal numerical backend performs
    low-level clustering mechanics:

    - `"exact_python"`: current exact Python implementation.
    - `"scipy_hierarchical"`: SciPy-backed hierarchical clustering backend.
    Alias: `clustering_engine`.

    `module_count`, when provided, is enforced strictly. It must be an integer
    between 1 and the number of available clustering sites for the current run.
    Invalid requests fail fast; the workflow does not clamp or rewrite the value.

    """

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
    cluster_tree_backend: SignalomeClusterTreeBackend = (
        SIGNALOME_CLUSTER_TREE_BACKEND_EXACT
    )
    tree_engine: SignalomeClusterTreeBackend | None = None
    candidate_scoring_backend: SignalomeCandidateScoringBackend = (
        SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL
    )
    candidate_scoring_policy: SignalomeCandidateScoringBackend | None = None
    max_exact_cluster_tree_sites: int = SIGNALOME_MAX_EXACT_CLUSTER_TREE_SITES_DEFAULT
    max_full_correlation_sites: int = SIGNALOME_MAX_FULL_CORRELATION_SITES_DEFAULT
    clustering_backend: SignalomeClusteringBackend = (
        SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON
    )
    clustering_engine: SignalomeClusteringBackend | None = None

    def __post_init__(self) -> None:
        resolved_cluster_tree_backend = _resolve_signalome_alias(
            legacy_field_name=(
                "signalome workflow request config.cluster_tree_backend"
            ),
            alias_field_name="signalome workflow request config.tree_engine",
            legacy_value=self.cluster_tree_backend,
            alias_value=self.tree_engine,
            legacy_default=SIGNALOME_CLUSTER_TREE_BACKEND_EXACT,
        )
        resolved_candidate_scoring_backend = _resolve_signalome_alias(
            legacy_field_name=(
                "signalome workflow request config.candidate_scoring_backend"
            ),
            alias_field_name=(
                "signalome workflow request config.candidate_scoring_policy"
            ),
            legacy_value=self.candidate_scoring_backend,
            alias_value=self.candidate_scoring_policy,
            legacy_default=SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
        )
        resolved_clustering_backend = _resolve_signalome_alias(
            legacy_field_name="signalome workflow request config.clustering_backend",
            alias_field_name="signalome workflow request config.clustering_engine",
            legacy_value=self.clustering_backend,
            alias_value=self.clustering_engine,
            legacy_default=SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON,
        )
        object.__setattr__(
            self,
            "cluster_tree_backend",
            resolved_cluster_tree_backend,
        )
        object.__setattr__(
            self,
            "tree_engine",
            resolved_cluster_tree_backend,
        )
        object.__setattr__(
            self,
            "candidate_scoring_backend",
            resolved_candidate_scoring_backend,
        )
        object.__setattr__(
            self,
            "candidate_scoring_policy",
            resolved_candidate_scoring_backend,
        )
        object.__setattr__(
            self,
            "clustering_backend",
            resolved_clustering_backend,
        )
        object.__setattr__(
            self,
            "clustering_engine",
            resolved_clustering_backend,
        )

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
        if self.cluster_tree_backend not in SIGNALOME_CLUSTER_TREE_BACKENDS:
            allowed_backends = ", ".join(sorted(SIGNALOME_CLUSTER_TREE_BACKENDS))
            raise WorkflowValidationError(
                "signalome workflow request config.cluster_tree_backend "
                f"must be one of: {allowed_backends}"
            )
        if self.candidate_scoring_backend not in SIGNALOME_CANDIDATE_SCORING_BACKENDS:
            allowed_backends = ", ".join(sorted(SIGNALOME_CANDIDATE_SCORING_BACKENDS))
            raise WorkflowValidationError(
                "signalome workflow request config.candidate_scoring_backend "
                f"must be one of: {allowed_backends}"
            )
        _require_int_at_least(
            self.max_exact_cluster_tree_sites,
            field_name="signalome workflow request config.max_exact_cluster_tree_sites",
            minimum=SIGNALOME_MAX_EXACT_CLUSTER_TREE_SITES_FLOOR,
            error_type=WorkflowValidationError,
        )
        _require_int_at_least(
            self.max_full_correlation_sites,
            field_name="signalome workflow request config.max_full_correlation_sites",
            minimum=SIGNALOME_MAX_FULL_CORRELATION_SITES_FLOOR,
            error_type=WorkflowValidationError,
        )
        if self.clustering_backend not in SIGNALOME_CLUSTERING_BACKENDS:
            allowed_backends = ", ".join(sorted(SIGNALOME_CLUSTERING_BACKENDS))
            raise WorkflowValidationError(
                "signalome workflow request config.clustering_backend "
                f"must be one of: {allowed_backends}"
            )


def _resolve_signalome_alias(
    *,
    legacy_field_name: str,
    alias_field_name: str,
    legacy_value: str,
    alias_value: str | None,
    legacy_default: str,
) -> str:
    if alias_value is None:
        return str(legacy_value)
    resolved_alias = str(alias_value)
    resolved_legacy = str(legacy_value)
    if resolved_legacy == str(legacy_default):
        return resolved_alias
    if resolved_alias != resolved_legacy:
        raise WorkflowValidationError(
            f"{legacy_field_name}={resolved_legacy!r} conflicts with "
            f"{alias_field_name}={resolved_alias!r}; provide only one value or use "
            "matching values."
        )
    return resolved_alias


__all__ = [
    "SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL",
    "SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED",
    "SIGNALOME_CANDIDATE_SCORING_BACKENDS",
    "SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON",
    "SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL",
    "SIGNALOME_CLUSTERING_BACKENDS",
    "SIGNALOME_CLUSTER_TREE_BACKEND_EXACT",
    "SIGNALOME_CLUSTER_TREE_BACKENDS",
    "SIGNALOME_MAX_EXACT_CLUSTER_TREE_SITES_DEFAULT",
    "SIGNALOME_MAX_EXACT_CLUSTER_TREE_SITES_FLOOR",
    "SIGNALOME_MAX_FULL_CORRELATION_SITES_DEFAULT",
    "SIGNALOME_MAX_FULL_CORRELATION_SITES_FLOOR",
    "SignalomeCandidateScoringBackend",
    "SignalomeClusterTreeBackend",
    "SIGNALOME_ASSIGNMENT_POLICIES",
    "SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY",
    "SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP",
    "SIGNALOME_KINASE_NETWORK_POLICIES",
    "SIGNALOME_KINASE_NETWORK_POLICY_ABSOLUTE_THRESHOLD",
    "SIGNALOME_KINASE_NETWORK_POLICY_POSITIVE_ONLY",
    "SIGNALOME_KINASE_NETWORK_POLICY_SIGNED",
    "SIGNALOME_MODULE_COUNT_FLOOR",
    "SIGNALOME_MODULE_SELECTION_FALLBACK_THRESHOLD_DEFAULT",
    "SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_DEFAULT",
    "SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_FLOOR",
    "SIGNALOME_MODULE_SELECTION_PRIMARY_THRESHOLD_DEFAULT",
    "SIGNALOME_SCORE_PRECONDITIONING_POLICIES",
    "SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT",
    "SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP",
    "SignalomeAssignmentPolicy",
    "SignalomeClusteringBackend",
    "SignalomeConfig",
    "SignalomeKinaseNetworkPolicy",
    "SignalomeScorePreconditioningPolicy",
]
