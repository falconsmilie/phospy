"""Tree-engine adaptation and engine-result assembly."""

from __future__ import annotations

from phospy.science.signalomes.clustering.contracts import ClusterTreeEngine
from phospy.science.signalomes.clustering.diagnostic_schemas import (
    SignalomeTreeEngineDiagnostics,
)
from phospy.science.signalomes.clustering.diagnostics import (
    SignalomeEngineDiagnosticsBuilder,
)
from phospy.science.signalomes.clustering.models import (
    SignalomeClusteringEngineRequest,
    SignalomeClusteringEngineResult,
)
from phospy.science.signalomes.clustering.orchestration import (
    cluster_sites_with_diagnostics,
)
from phospy.science.signalomes.clustering.policies import SIGNALOME_TREE_ENGINE_EXACT
from phospy.science.signalomes.clustering.protein_modules import derive_protein_modules
from phospy.science.signalomes.clustering.tree_building import (
    ClusterTreeOperationsAdapter,
)

_DIAGNOSTICS_BUILDER = SignalomeEngineDiagnosticsBuilder()


def run_clustering_with_tree_engine(
    *,
    request: SignalomeClusteringEngineRequest,
    tree_engine: ClusterTreeEngine,
    clustering_engine: str,
    backend_version: str,
    backend_diagnostics: SignalomeTreeEngineDiagnostics,
) -> SignalomeClusteringEngineResult:
    """Run shared orchestration with an injected tree engine implementation."""

    clustering_result = cluster_sites_with_diagnostics(
        scoring_matrix=request.scoring_matrix,
        requested_module_count=request.requested_module_count,
        primary_threshold=request.primary_threshold,
        fallback_threshold=request.fallback_threshold,
        max_clusters=request.max_clusters,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=request.candidate_scoring_policy,
        max_exact_tree_sites=request.max_exact_tree_sites,
        max_full_candidate_scoring_sites=request.max_full_candidate_scoring_sites,
        cluster_tree_operations=ClusterTreeOperationsAdapter(engine=tree_engine),
    )
    protein_modules = derive_protein_modules(
        site_clusters=clustering_result.site_clusters,
        site_to_protein=request.site_to_protein,
    )
    selected_module_count = int(
        clustering_result.module_selection_diagnostics.selected_module_count
    )
    resolved_backend_diagnostics = _DIAGNOSTICS_BUILDER.backend_diagnostics(
        clustering_engine=clustering_engine,
        backend_version=backend_version,
        tree_engine=tree_engine,
        tree_engine_diagnostics=backend_diagnostics,
        selected_module_count=selected_module_count,
        input_site_count=int(request.scoring_matrix.shape[0]),
        exact_tree_path_used=bool(clustering_result.exact_cluster_tree_built),
    )
    threshold_metadata = _DIAGNOSTICS_BUILDER.threshold_metadata(
        primary_threshold=request.primary_threshold,
        fallback_threshold=request.fallback_threshold,
    )
    limit_metadata = _DIAGNOSTICS_BUILDER.limit_metadata(
        max_exact_tree_sites=request.max_exact_tree_sites,
        max_full_candidate_scoring_sites=request.max_full_candidate_scoring_sites,
        max_clusters=request.max_clusters,
    )
    return SignalomeClusteringEngineResult(
        site_clusters=clustering_result.site_clusters,
        protein_modules=protein_modules,
        selected_module_count=selected_module_count,
        module_selection_diagnostics=clustering_result.module_selection_diagnostics,
        clustering_preparation_diagnostics=(
            clustering_result.clustering_preparation_diagnostics
        ),
        backend_name=str(clustering_engine),
        backend_version=str(backend_version),
        approximation_used=bool(clustering_result.approximation_used),
        exact_cluster_tree_built=bool(clustering_result.exact_cluster_tree_built),
        tree_implementation=str(tree_engine.name),
        candidate_scoring_mode=str(clustering_result.candidate_scoring_mode),
        candidate_scoring_evaluated=bool(clustering_result.candidate_scoring_evaluated),
        candidate_scoring_skip_reason=clustering_result.candidate_scoring_skip_reason,
        candidate_scoring_sampling=clustering_result.candidate_scoring_sampling,
        backend_diagnostics=resolved_backend_diagnostics,
        threshold_metadata=threshold_metadata,
        limit_metadata=limit_metadata,
    )


__all__ = ["run_clustering_with_tree_engine"]
