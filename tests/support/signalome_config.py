from __future__ import annotations

from phospy.api.configs import (
    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
    SIGNALOME_KINASE_NETWORK_POLICY_SIGNED,
    SIGNALOME_MAX_EXACT_TREE_SITES_DEFAULT,
    SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_DEFAULT,
    SIGNALOME_MODULE_SELECTION_FALLBACK_THRESHOLD_DEFAULT,
    SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_DEFAULT,
    SIGNALOME_MODULE_SELECTION_PRIMARY_THRESHOLD_DEFAULT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    SIGNALOME_TREE_ENGINE_EXACT,
    SignalomeAssignmentPolicy,
    SignalomeCandidateScoringPolicy,
    SignalomeClusteringConfig,
    SignalomeClusteringEngine,
    SignalomeConfig,
    SignalomeKinaseNetworkPolicy,
    SignalomeOutputConfig,
    SignalomePerformanceConfig,
    SignalomeScientificConfig,
    SignalomeScorePreconditioningPolicy,
    SignalomeTreeEngine,
    SignalomeValidationConfig,
)


def build_signalome_config(
    *,
    substrate_support_cutoff: float = 0.5,
    network_correlation_threshold: float = 0.5,
    network_policy: SignalomeKinaseNetworkPolicy = SIGNALOME_KINASE_NETWORK_POLICY_SIGNED,
    assignment_policy: SignalomeAssignmentPolicy = SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    score_preconditioning_policy: SignalomeScorePreconditioningPolicy = (
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
    ),
    allow_mixed_total_protein_quantitative_meaning: bool = False,
    module_count: int | None = None,
    module_selection_primary_correlation_threshold: float = (
        SIGNALOME_MODULE_SELECTION_PRIMARY_THRESHOLD_DEFAULT
    ),
    module_selection_fallback_correlation_threshold: float = (
        SIGNALOME_MODULE_SELECTION_FALLBACK_THRESHOLD_DEFAULT
    ),
    module_selection_max_clusters: int = SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_DEFAULT,
    tree_engine: SignalomeTreeEngine = SIGNALOME_TREE_ENGINE_EXACT,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy = (
        SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    ),
    max_exact_tree_sites: int = SIGNALOME_MAX_EXACT_TREE_SITES_DEFAULT,
    max_full_candidate_scoring_sites: int = (
        SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_DEFAULT
    ),
    clustering_engine: SignalomeClusteringEngine = (
        SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL
    ),
) -> SignalomeConfig:
    return SignalomeConfig(
        scientific=SignalomeScientificConfig(
            substrate_support_cutoff=substrate_support_cutoff,
            assignment_policy=assignment_policy,
        ),
        clustering=SignalomeClusteringConfig(
            module_count=module_count,
            module_selection_primary_correlation_threshold=(
                module_selection_primary_correlation_threshold
            ),
            module_selection_fallback_correlation_threshold=(
                module_selection_fallback_correlation_threshold
            ),
            module_selection_max_clusters=module_selection_max_clusters,
            tree_engine=tree_engine,
            candidate_scoring_policy=candidate_scoring_policy,
            clustering_engine=clustering_engine,
        ),
        validation=SignalomeValidationConfig(
            score_preconditioning_policy=score_preconditioning_policy,
            allow_mixed_total_protein_quantitative_meaning=(
                allow_mixed_total_protein_quantitative_meaning
            ),
        ),
        output=SignalomeOutputConfig(
            network_correlation_threshold=network_correlation_threshold,
            network_policy=network_policy,
        ),
        performance=SignalomePerformanceConfig(
            max_exact_tree_sites=max_exact_tree_sites,
            max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        ),
    )
