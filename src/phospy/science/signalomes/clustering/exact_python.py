"""Compatibility facade for the exact-Python clustering backend."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.science.signalomes.clustering.backends.exact_python import (
    ExactPythonTreeEngine,
    ExactWardClusterTree,
)
from phospy.science.signalomes.clustering.candidate_scoring import (
    build_correlation_exclusion_note,
    build_correlation_matrix_with_exclusions,
    cluster_median_correlation,
    cluster_median_correlation_approximate,
    summarize_profile_degeneracy,
)
from phospy.science.signalomes.clustering.candidate_selection import (
    filter_cluster_candidates,
)
from phospy.science.signalomes.clustering.diagnostic_schemas import (
    build_tree_engine_diagnostics,
)
from phospy.science.signalomes.clustering.models import (
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION,
    SignalomeClusteringEngineRequest,
    SignalomeClusteringEngineResult,
)
from phospy.science.signalomes.clustering.orchestration import (
    ClusterSitesResult,
    cluster_sites,
    cluster_sites_with_diagnostics,
    select_module_count,
    select_module_count_with_diagnostics,
)
from phospy.science.signalomes.clustering.policies import (
    MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER,
    MAX_FULL_CORRELATION_SITE_COUNT,
    NEAR_CONSTANT_PROFILE_VARIANCE_TOLERANCE,
    SIGNALOME_CANDIDATE_SCORING_APPLIES_TO,
    SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD,
    SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY,
    SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT,
    SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_APPLIES_TO,
    SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_DROP_FULLY_MISSING_THEN_COLUMN_MEDIAN_IMPUTE,
    SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_IMPUTED_VALUES_EXPOSED_IN_OUTPUT_TABLES,
    SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE,
    SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
    SIGNALOME_CLUSTERING_SCORING_MODE_EXACT,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE,
    SIGNALOME_MODULE_SELECTION_STABILITY_ASSIGNMENT_METRIC,
    SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_MAX_SITES,
    SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_PERTURBATIONS,
    SIGNALOME_MODULE_SELECTION_STABILITY_PERTURBATION_SCALE,
    SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_CALLER_FIXED,
    SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_INPUT_DERIVED,
    SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_NOT_APPLICABLE,
    SIGNALOME_MODULE_SELECTION_STABILITY_THRESHOLD_DELTA,
    SIGNALOME_TREE_ENGINE_EXACT,
    SignalomeCandidateScoringPolicy,
    SignalomeClusteringMissingValuePolicy,
    SignalomeClusteringScoringMode,
    SignalomeTreeEngine,
)
from phospy.science.signalomes.clustering.protein_modules import derive_protein_modules
from phospy.science.signalomes.clustering.tree_building import (
    ClusterTreeOperations,
    build_cluster_labels_from_tree,
    fit_cluster_labels,
)
from phospy.science.signalomes.clustering.tree_building import (
    build_cluster_tree as _compat_build_cluster_tree,
)
from phospy.science.signalomes.clustering.tree_engine_adapter import (
    run_clustering_with_tree_engine,
)

# Backward-compatible private aliases used by performance/benchmark contracts.
_WardClusterTree = ExactWardClusterTree
_build_cluster_tree = _compat_build_cluster_tree


@dataclass(frozen=True, slots=True)
class ExactPythonClusteringBackend:
    """Top-level exact-Python clustering backend (shared orchestration + engine)."""

    name: str = SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON
    version: str = SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION

    def run(
        self,
        request: SignalomeClusteringEngineRequest,
    ) -> SignalomeClusteringEngineResult:
        return run_clustering_with_tree_engine(
            request=request,
            tree_engine=ExactPythonTreeEngine(),
            clustering_engine=self.name,
            backend_version=self.version,
            backend_diagnostics=build_tree_engine_diagnostics(
                uses_scipy=False,
                linkage_method="ward",
                distance_metric="euclidean",
            ),
        )


__all__ = [
    "ClusterSitesResult",
    "ClusterTreeOperations",
    "ExactPythonClusteringBackend",
    "MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER",
    "MAX_FULL_CORRELATION_SITE_COUNT",
    "NEAR_CONSTANT_PROFILE_VARIANCE_TOLERANCE",
    "SIGNALOME_CANDIDATE_SCORING_APPLIES_TO",
    "SIGNALOME_CANDIDATE_SCORING_POLICY_FULL",
    "SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED",
    "SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED",
    "SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD",
    "SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY",
    "SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT",
    "SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_APPLIES_TO",
    "SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_DROP_FULLY_MISSING_THEN_COLUMN_MEDIAN_IMPUTE",
    "SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_IMPUTED_VALUES_EXPOSED_IN_OUTPUT_TABLES",
    "SIGNALOME_TREE_ENGINE_EXACT",
    "SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE",
    "SIGNALOME_CLUSTERING_SCORING_MODE_AUTO",
    "SIGNALOME_CLUSTERING_SCORING_MODE_EXACT",
    "SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE",
    "SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE",
    "SIGNALOME_MODULE_SELECTION_STABILITY_ASSIGNMENT_METRIC",
    "SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_MAX_SITES",
    "SIGNALOME_MODULE_SELECTION_STABILITY_DEFAULT_PERTURBATIONS",
    "SIGNALOME_MODULE_SELECTION_STABILITY_PERTURBATION_SCALE",
    "SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_CALLER_FIXED",
    "SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_INPUT_DERIVED",
    "SIGNALOME_MODULE_SELECTION_STABILITY_SEED_POLICY_NOT_APPLICABLE",
    "SIGNALOME_MODULE_SELECTION_STABILITY_THRESHOLD_DELTA",
    "SignalomeCandidateScoringPolicy",
    "SignalomeClusteringMissingValuePolicy",
    "SignalomeTreeEngine",
    "SignalomeClusteringScoringMode",
    "build_cluster_labels_from_tree",
    "build_correlation_exclusion_note",
    "build_correlation_matrix_with_exclusions",
    "cluster_median_correlation",
    "cluster_median_correlation_approximate",
    "cluster_sites",
    "cluster_sites_with_diagnostics",
    "derive_protein_modules",
    "filter_cluster_candidates",
    "fit_cluster_labels",
    "select_module_count",
    "select_module_count_with_diagnostics",
    "summarize_profile_degeneracy",
]
