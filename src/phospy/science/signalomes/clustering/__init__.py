"""Signalome clustering public facade and backend boundary."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.science.signalomes.clustering import exact_python as _exact
from phospy.science.signalomes.clustering import selection as _selection
from phospy.science.signalomes.clustering.backend_dispatch import (
    available_clustering_engines,
    resolve_clustering_engine,
    run_clustering_engine,
)
from phospy.science.signalomes.clustering.exact_python import (
    MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER,
    NEAR_CONSTANT_PROFILE_VARIANCE_TOLERANCE,
    SIGNALOME_CANDIDATE_SCORING_APPLIES_TO,
    SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD,
    SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY,
    SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE,
    ClusterSitesResult,
    SignalomeCandidateScoringPolicy,
    SignalomeClusteringScoringMode,
    SignalomeTreeEngine,
    derive_protein_modules,
)
from phospy.science.signalomes.clustering.models import (
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION,
    SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
    SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL_VERSION,
    SignalomeClusteringEngineRequest,
    SignalomeClusteringEngineResult,
)
from phospy.science.signalomes.clustering.policies import (
    MAX_FULL_CORRELATION_SITE_COUNT,
    SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_APPLIES_TO,
    SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_DROP_FULLY_MISSING_THEN_COLUMN_MEDIAN_IMPUTE,
    SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_IMPUTED_VALUES_EXPOSED_IN_OUTPUT_TABLES,
    SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE,
    SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
    SIGNALOME_CLUSTERING_SCORING_MODE_EXACT,
    SIGNALOME_MODULE_SELECTION_STABILITY_ASSIGNMENT_METRIC,
    SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_MAX_SITES,
    SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_PERTURBATIONS,
    SIGNALOME_MODULE_SELECTION_STABILITY_PERTURBATION_SCALE,
    SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_CALLER_FIXED,
    SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_INPUT_DERIVED,
    SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_NOT_APPLICABLE,
    SIGNALOME_MODULE_SELECTION_STABILITY_THRESHOLD_DELTA,
    SIGNALOME_TREE_ENGINE_EXACT,
    SignalomeClusteringMissingValuePolicy,
)
from phospy.science.signalomes.clustering.protocol import SignalomeClusteringEngine
from phospy.science.signalomes.clustering.tree_building import (
    build_cluster_labels_from_tree,
)

# Backward-compatible private alias for older internal callers/tests.
_build_cluster_tree = _exact._build_cluster_tree


def run_signalome_clustering_engine(
    *,
    scoring_matrix: pd.DataFrame,
    site_to_protein_group_id: pd.Series | None = None,
    site_to_protein: pd.Series | None = None,
    requested_module_count: int | None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy | None = None,
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    max_full_candidate_scoring_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
    clustering_engine: str = SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
    module_selection_stability_perturbations: int = (
        SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_PERTURBATIONS
    ),
    module_selection_stability_seed: int | None = None,
    module_selection_stability_max_sites: int = (
        SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_MAX_SITES
    ),
) -> SignalomeClusteringEngineResult:
    """Run configured signalome clustering backend."""

    return run_clustering_engine(
        request=SignalomeClusteringEngineRequest(
            scoring_matrix=scoring_matrix,
            site_to_protein_group_id=site_to_protein_group_id,
            site_to_protein=site_to_protein,
            requested_module_count=requested_module_count,
            primary_threshold=primary_threshold,
            fallback_threshold=fallback_threshold,
            max_clusters=max_clusters,
            candidate_scoring_policy=candidate_scoring_policy,
            max_exact_tree_sites=max_exact_tree_sites,
            max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
            module_selection_stability_perturbations=(
                module_selection_stability_perturbations
            ),
            module_selection_stability_seed=module_selection_stability_seed,
            module_selection_stability_max_sites=(module_selection_stability_max_sites),
        ),
        clustering_engine=clustering_engine,
    )


def cluster_sites(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    tree_engine: SignalomeTreeEngine = SIGNALOME_TREE_ENGINE_EXACT,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy | None = None,
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    max_full_candidate_scoring_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
    module_selection_stability_perturbations: int = (
        SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_PERTURBATIONS
    ),
    module_selection_stability_seed: int | None = None,
    module_selection_stability_max_sites: int = (
        SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_MAX_SITES
    ),
) -> pd.Series:
    return _exact.cluster_sites(
        scoring_matrix=scoring_matrix,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        tree_engine=tree_engine,
        candidate_scoring_policy=candidate_scoring_policy,
        max_exact_tree_sites=max_exact_tree_sites,
        max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        module_selection_stability_perturbations=(
            module_selection_stability_perturbations
        ),
        module_selection_stability_seed=module_selection_stability_seed,
        module_selection_stability_max_sites=module_selection_stability_max_sites,
    )


def cluster_sites_with_diagnostics(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    tree_engine: SignalomeTreeEngine = SIGNALOME_TREE_ENGINE_EXACT,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy | None = None,
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    max_full_candidate_scoring_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
    module_selection_stability_perturbations: int = (
        SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_PERTURBATIONS
    ),
    module_selection_stability_seed: int | None = None,
    module_selection_stability_max_sites: int = (
        SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_MAX_SITES
    ),
) -> ClusterSitesResult:
    return _exact.cluster_sites_with_diagnostics(
        scoring_matrix=scoring_matrix,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        tree_engine=tree_engine,
        candidate_scoring_policy=candidate_scoring_policy,
        max_exact_tree_sites=max_exact_tree_sites,
        max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        module_selection_stability_perturbations=(
            module_selection_stability_perturbations
        ),
        module_selection_stability_seed=module_selection_stability_seed,
        module_selection_stability_max_sites=module_selection_stability_max_sites,
    )


def select_module_count(
    scoring_values: pd.DataFrame | np.ndarray,
    *,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    module_selection_stability_perturbations: int = (
        SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_PERTURBATIONS
    ),
    module_selection_stability_seed: int | None = None,
    module_selection_stability_max_sites: int = (
        SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_MAX_SITES
    ),
) -> int:
    return _selection.select_module_count(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        max_exact_tree_sites=max_exact_tree_sites,
        module_selection_stability_perturbations=(
            module_selection_stability_perturbations
        ),
        module_selection_stability_seed=module_selection_stability_seed,
        module_selection_stability_max_sites=module_selection_stability_max_sites,
    )


def select_module_count_with_diagnostics(
    *,
    scoring_values: pd.DataFrame | np.ndarray,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    module_selection_stability_perturbations: int = (
        SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_PERTURBATIONS
    ),
    module_selection_stability_seed: int | None = None,
    module_selection_stability_max_sites: int = (
        SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_MAX_SITES
    ),
):
    return _selection.select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        max_exact_tree_sites=max_exact_tree_sites,
        module_selection_stability_perturbations=(
            module_selection_stability_perturbations
        ),
        module_selection_stability_seed=module_selection_stability_seed,
        module_selection_stability_max_sites=module_selection_stability_max_sites,
    )


def fit_cluster_labels(
    scoring_values: np.ndarray,
    cluster_count: int,
    *,
    tree_engine: SignalomeTreeEngine = SIGNALOME_TREE_ENGINE_EXACT,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy = SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
) -> np.ndarray:
    return _exact.fit_cluster_labels(
        scoring_values=scoring_values,
        cluster_count=cluster_count,
        tree_engine=tree_engine,
        candidate_scoring_policy=candidate_scoring_policy,
        max_exact_tree_sites=max_exact_tree_sites,
    )


__all__ = [
    "ClusterSitesResult",
    "MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER",
    "MAX_FULL_CORRELATION_SITE_COUNT",
    "NEAR_CONSTANT_PROFILE_VARIANCE_TOLERANCE",
    "SIGNALOME_CANDIDATE_SCORING_APPLIES_TO",
    "SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED",
    "SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_APPLIES_TO",
    "SIGNALOME_CANDIDATE_SCORING_POLICY_FULL",
    "SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED",
    "SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD",
    "SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY",
    "SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT",
    "SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_DROP_FULLY_MISSING_THEN_COLUMN_MEDIAN_IMPUTE",
    "SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_IMPUTED_VALUES_EXPOSED_IN_OUTPUT_TABLES",
    "SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE",
    "SIGNALOME_CLUSTERING_SCORING_MODE_AUTO",
    "SIGNALOME_CLUSTERING_SCORING_MODE_EXACT",
    "SIGNALOME_MODULE_SELECTION_STABILITY_ASSIGNMENT_METRIC",
    "SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_MAX_SITES",
    "SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_PERTURBATIONS",
    "SIGNALOME_MODULE_SELECTION_STABILITY_PERTURBATION_SCALE",
    "SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_CALLER_FIXED",
    "SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_INPUT_DERIVED",
    "SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_NOT_APPLICABLE",
    "SIGNALOME_MODULE_SELECTION_STABILITY_THRESHOLD_DELTA",
    "SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON",
    "SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION",
    "SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL",
    "SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL_VERSION",
    "SIGNALOME_TREE_ENGINE_EXACT",
    "SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE",
    "SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE",
    "SignalomeClusteringEngine",
    "SignalomeClusteringMissingValuePolicy",
    "SignalomeClusteringEngineRequest",
    "SignalomeClusteringEngineResult",
    "available_clustering_engines",
    "build_cluster_labels_from_tree",
    "cluster_sites",
    "cluster_sites_with_diagnostics",
    "derive_protein_modules",
    "fit_cluster_labels",
    "resolve_clustering_engine",
    "select_module_count",
    "select_module_count_with_diagnostics",
    "run_clustering_engine",
    "run_signalome_clustering_engine",
]
