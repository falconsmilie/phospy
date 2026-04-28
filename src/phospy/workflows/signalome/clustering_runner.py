"""Clustering orchestration for signalome workflow execution."""

from __future__ import annotations

from collections.abc import Callable
from typing import cast

import pandas as pd

from phospy.errors.workflows import (
    SignalomeModuleCountValidationError,
    WorkflowStageError,
)
from phospy.signalomes.clustering import (
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE,
    ClusterSitesResult,
    cluster_sites_with_diagnostics,
    derive_protein_modules,
    run_signalome_clustering_backend,
)
from phospy.signalomes.clustering.exact_python import (
    SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
    SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED,
    SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED,
    _CandidateScoringMode,
)
from phospy.signalomes.clustering.models import SignalomeClusteringBackendResult
from phospy.signalomes.clustering.validation import validate_requested_module_count
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
        run_backend_clustering: Callable[
            ..., SignalomeClusteringBackendResult
        ] = run_signalome_clustering_backend,
    ) -> None:
        self._cluster_sites = cluster_sites
        self._derive_modules = derive_modules
        self._run_backend_clustering = run_backend_clustering
        self._use_legacy_injected_cluster_functions = (
            cluster_sites is not cluster_sites_with_diagnostics
            or derive_modules is not derive_protein_modules
        )

    def run(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        execution_metadata: SignalomeExecutionMetadata,
    ) -> SignalomeClusteringRunResult:
        try:
            validated_requested_module_count = validate_requested_module_count(
                requested_module_count=config.requested_module_count,
                available_clustering_site_count=int(
                    request.downstream_score_matrix.shape[0]
                ),
                field_name="signalome workflow request config.module_count",
            )
        except SignalomeModuleCountValidationError as exc:
            raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "choose config.module_count between 1 and the number of "
                    "available clustering sites, or omit config.module_count for "
                    "automatic module-count selection"
                ),
                requested_module_count=requested_module_count_label(
                    config.requested_module_count
                ),
                available_clustering_site_count=int(
                    request.downstream_score_matrix.shape[0]
                ),
                affected_configuration_field="signalome workflow request config.module_count",
                validation_error=str(exc),
            )
        try:
            if self._use_legacy_injected_cluster_functions:
                clustering_result = self._cluster_sites(
                    scoring_matrix=request.downstream_score_matrix,
                    requested_module_count=validated_requested_module_count,
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
            else:
                backend_result = self._run_backend_clustering(
                    scoring_matrix=request.downstream_score_matrix,
                    site_to_protein=request.site_to_protein,
                    requested_module_count=validated_requested_module_count,
                    primary_threshold=config.module_selection_primary_threshold,
                    fallback_threshold=config.module_selection_fallback_threshold,
                    max_clusters=config.module_selection_max_clusters,
                    cluster_tree_backend=config.cluster_tree_backend,
                    candidate_scoring_backend=config.candidate_scoring_backend,
                    max_exact_cluster_tree_sites=config.max_exact_cluster_tree_sites,
                    max_full_correlation_sites=config.max_full_correlation_sites,
                    backend_name=config.clustering_backend,
                    tree_engine=config.tree_engine,
                    candidate_scoring_policy=config.candidate_scoring_policy,
                    clustering_engine=config.clustering_engine,
                )
                if backend_result.candidate_scoring_mode not in {
                    SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
                    SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED,
                    SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED,
                }:
                    raise ValueError(
                        "backend returned unsupported candidate_scoring_mode: "
                        f"{backend_result.candidate_scoring_mode!r}"
                    )
                clustering_result = ClusterSitesResult(
                    site_clusters=backend_result.site_clusters,
                    module_selection_diagnostics=backend_result.module_selection_diagnostics,
                    cluster_tree_backend=backend_result.cluster_tree_backend,
                    candidate_scoring_mode=cast(
                        _CandidateScoringMode,
                        backend_result.candidate_scoring_mode,
                    ),
                    exact_cluster_tree_built=backend_result.exact_cluster_tree_built,
                    candidate_scoring_sampling=backend_result.candidate_scoring_sampling,
                    candidate_scoring_evaluated=backend_result.candidate_scoring_evaluated,
                    candidate_scoring_skip_reason=backend_result.candidate_scoring_skip_reason,
                    backend_name=backend_result.backend_name,
                    backend_version=backend_result.backend_version,
                    approximation_used=backend_result.approximation_used,
                    backend_diagnostics=backend_result.backend_diagnostics,
                )
                protein_modules = backend_result.protein_modules
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
        selected_module_count = int(
            clustering_result.module_selection_diagnostics.selected_module_count
        )
        final_module_assignment_backend = (
            SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE
            if selected_module_count <= 1
            else SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE
        )
        return SignalomeScaleGuardDecision(
            site_count=int(site_count),
            selected_module_count=selected_module_count,
            clustering_backend=str(clustering_result.backend_name),
            clustering_backend_version=str(clustering_result.backend_version),
            backend_diagnostics=clustering_result.backend_diagnostics,
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
            final_module_assignment_backend=final_module_assignment_backend,
        )


__all__ = ["SignalomeClusteringRunner"]
