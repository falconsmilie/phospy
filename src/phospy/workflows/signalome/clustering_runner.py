"""Clustering orchestration for signalome workflow execution."""

from __future__ import annotations

from collections.abc import Callable

import pandas as pd

from phospy.errors.workflows import WorkflowStageError
from phospy.signalomes.clustering import (
    ClusterSitesResult,
    cluster_sites_with_diagnostics,
    derive_protein_modules,
)
from phospy.workflows.signalome.component_helpers import (
    raise_boundary_error,
    requested_module_count_label,
)
from phospy.workflows.signalome.component_models import (
    SignalomeClusteringRunResult,
    SignalomeExecutionMetadata,
    SignalomeScaleGuardDecision,
)
from phospy.workflows.signalome.constants import (
    SIGNALOME_EXECUTOR_MODULE_CONSTRUCTION_SEAM,
)
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeExecutionConfig,
    ResolvedSignalomeWorkflowRequest,
)


class SignalomeClusteringRunner:
    """Run clustering and protein-module derivation for signalome execution."""

    def __init__(
        self,
        *,
        cluster_sites: Callable[..., ClusterSitesResult] = (
            cluster_sites_with_diagnostics
        ),
        derive_modules: Callable[..., pd.Series] = derive_protein_modules,
    ) -> None:
        self._cluster_sites = cluster_sites
        self._derive_modules = derive_modules

    def run(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        execution_metadata: SignalomeExecutionMetadata,
    ) -> SignalomeClusteringRunResult:
        try:
            clustering_result = self._cluster_sites(
                scoring_matrix=request.downstream_score_matrix,
                requested_module_count=config.requested_module_count,
                primary_threshold=config.module_selection_primary_threshold,
                fallback_threshold=config.module_selection_fallback_threshold,
                max_clusters=config.module_selection_max_clusters,
                cluster_tree_backend=config.cluster_tree_backend,
                candidate_scoring_backend=config.candidate_scoring_backend,
                max_exact_cluster_tree_sites=config.max_exact_cluster_tree_sites,
                max_full_correlation_sites=config.max_full_correlation_sites,
            )
            protein_modules = self._derive_modules(
                site_clusters=clustering_result.site_clusters,
                site_to_protein=request.site_to_protein,
            )
            return SignalomeClusteringRunResult(
                clustering_result=clustering_result,
                protein_modules=protein_modules,
            )
        except (WorkflowStageError, ValueError) as exc:
            raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "ensure downstream signalome scores are compatible with module "
                    "selection clustering and module-count policy settings"
                ),
                requested_module_count=requested_module_count_label(
                    config.requested_module_count
                ),
                module_selection_primary_correlation_threshold=config.module_selection_primary_threshold,
                module_selection_fallback_correlation_threshold=config.module_selection_fallback_threshold,
                module_selection_max_clusters=config.module_selection_max_clusters,
                cluster_tree_backend=config.cluster_tree_backend,
                candidate_scoring_backend=config.candidate_scoring_backend,
                max_exact_cluster_tree_sites=config.max_exact_cluster_tree_sites,
                max_full_correlation_sites=config.max_full_correlation_sites,
                downstream_score_sites=execution_metadata.downstream_score_sites,
                downstream_score_kinases=execution_metadata.downstream_score_kinases,
                stage_error=str(exc),
            )

    @staticmethod
    def collect_execution_metadata(
        request: ResolvedSignalomeWorkflowRequest,
    ) -> SignalomeExecutionMetadata:
        return SignalomeExecutionMetadata(
            prediction_sites=int(request.prediction_matrix.shape[0]),
            prediction_kinases=int(request.prediction_matrix.shape[1]),
            downstream_score_sites=int(request.downstream_score_matrix.shape[0]),
            downstream_score_kinases=int(request.downstream_score_matrix.shape[1]),
            downstream_score_source=request.downstream_score_source,
        )

    @staticmethod
    def summarize_scale_guard(
        *,
        config: ResolvedSignalomeExecutionConfig,
        site_count: int,
        clustering_result: ClusterSitesResult,
    ) -> SignalomeScaleGuardDecision:
        return SignalomeScaleGuardDecision(
            site_count=int(site_count),
            cluster_tree_backend=str(config.cluster_tree_backend),
            candidate_scoring_backend=str(config.candidate_scoring_backend),
            candidate_scoring_requested_backend=str(config.candidate_scoring_backend),
            max_exact_cluster_tree_sites=int(config.max_exact_cluster_tree_sites),
            max_full_correlation_sites=int(config.max_full_correlation_sites),
            exact_cluster_tree_built=bool(clustering_result.exact_cluster_tree_built),
            candidate_scoring_mode=str(clustering_result.candidate_scoring_mode),
            candidate_scoring_evaluated=bool(
                clustering_result.candidate_scoring_evaluated
            ),
            candidate_scoring_skip_reason=(
                None
                if clustering_result.candidate_scoring_skip_reason is None
                else str(clustering_result.candidate_scoring_skip_reason)
            ),
            candidate_scoring_sampling=clustering_result.candidate_scoring_sampling,
            scale_guard_passed=True,
        )


__all__ = ["SignalomeClusteringRunner"]
