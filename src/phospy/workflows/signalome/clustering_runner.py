"""Clustering orchestration for signalome workflow execution."""

from __future__ import annotations

from collections.abc import Callable

import pandas as pd

from phospy.errors.workflows import (
    SignalomeModuleCountValidationError,
    WorkflowStageError,
)
from phospy.science.signalomes.clustering import (
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE,
    ClusterSitesResult,
    run_signalome_clustering_engine,
)
from phospy.science.signalomes.clustering.diagnostic_schemas import (
    SignalomeBackendDiagnostics,
    SignalomeCandidateScoringSamplingDiagnostics,
)
from phospy.science.signalomes.clustering.exact_python import (
    SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
)
from phospy.science.signalomes.clustering.models import SignalomeClusteringEngineResult
from phospy.science.signalomes.clustering.policies import (
    _CandidateScoringMode,
)
from phospy.science.signalomes.clustering.validation import (
    validate_requested_module_count,
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
        run_backend_clustering: Callable[
            ..., SignalomeClusteringEngineResult
        ] = run_signalome_clustering_engine,
    ) -> None:
        self._run_backend_clustering = run_backend_clustering

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
                field_name="signalome workflow request config.clustering.module_count",
            )
        except SignalomeModuleCountValidationError as exc:
            raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "choose config.clustering.module_count between 1 and the number "
                    "of available clustering sites, or omit "
                    "config.clustering.module_count for automatic module-count "
                    "selection"
                ),
                requested_module_count=requested_module_count_label(
                    config.requested_module_count
                ),
                available_clustering_site_count=int(
                    request.downstream_score_matrix.shape[0]
                ),
                affected_configuration_field=(
                    "signalome workflow request config.clustering.module_count"
                ),
                validation_error=str(exc),
            )
        try:
            backend_result = self._run_backend_clustering(
                scoring_matrix=request.downstream_score_matrix,
                site_to_protein=request.site_to_protein,
                requested_module_count=validated_requested_module_count,
                primary_threshold=config.module_selection_primary_threshold,
                fallback_threshold=config.module_selection_fallback_threshold,
                max_clusters=config.module_selection_max_clusters,
                candidate_scoring_policy=config.candidate_scoring_policy,
                max_exact_tree_sites=config.max_exact_tree_sites,
                max_full_candidate_scoring_sites=config.max_full_candidate_scoring_sites,
                clustering_engine=config.clustering_engine,
            )
            if backend_result.candidate_scoring_mode not in {
                SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
                SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
                SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED,
            }:
                raise ValueError(
                    "backend returned unsupported candidate_scoring_mode: "
                    f"{backend_result.candidate_scoring_mode!r}"
                )
            clustering_result = ClusterSitesResult(
                site_clusters=backend_result.site_clusters,
                module_selection_diagnostics=backend_result.module_selection_diagnostics,
                clustering_preparation_diagnostics=(
                    backend_result.clustering_preparation_diagnostics
                ),
                tree_engine=backend_result.tree_implementation,
                candidate_scoring_mode=_validated_candidate_scoring_mode(
                    backend_result.candidate_scoring_mode
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
                candidate_scoring_policy=config.candidate_scoring_policy,
                max_exact_tree_sites=config.max_exact_tree_sites,
                max_full_candidate_scoring_sites=config.max_full_candidate_scoring_sites,
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
            downstream_score_selection_policy=request.downstream_score_selection_policy,
        )

    @staticmethod
    def summarize_scale_guard(
        *,
        config: ResolvedSignalomeExecutionConfig,
        site_count: int,
        site_to_protein: pd.Series,
        downstream_score_kinases: int,
        clustering_result: ClusterSitesResult,
    ) -> SignalomeScaleGuardDecision:
        def _sampled_sites_total(
            payload: SignalomeCandidateScoringSamplingDiagnostics | None,
        ) -> int | None:
            if payload is None:
                return None
            return int(payload["per_cluster_sample_count_summary"]["total"])

        def _sampled_pairs_total(
            payload: SignalomeCandidateScoringSamplingDiagnostics | None,
        ) -> int | None:
            if payload is None:
                return None
            return int(payload["actual_sampled_pair_count"])

        def _tree_generation_backend_name(
            payload: SignalomeBackendDiagnostics | None,
        ) -> str:
            if payload is None:
                return "unknown_tree_backend"
            return str(payload["tree_implementation"])

        input_protein_count = int(
            pd.Index(site_to_protein.astype(str)).nunique(dropna=True)
        )
        input_kinase_count = int(downstream_score_kinases)
        candidate_scores = (
            clustering_result.module_selection_diagnostics.candidate_scores
        )
        candidate_module_counts_evaluated = int(len(candidate_scores))
        candidate_module_count_upper_bound = int(
            clustering_result.module_selection_diagnostics.max_clusters_evaluated
        )
        selected_module_count = int(
            clustering_result.module_selection_diagnostics.selected_module_count
        )
        final_module_assignment_backend = (
            SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE
            if selected_module_count <= 1
            else SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE
        )
        sampled_candidate_scoring = clustering_result.candidate_scoring_sampling
        candidate_scoring_is_approximate = bool(
            clustering_result.candidate_scoring_evaluated
            and str(clustering_result.candidate_scoring_mode)
            == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
        )
        tree_generation_backend = _tree_generation_backend_name(
            clustering_result.backend_diagnostics
        )
        return SignalomeScaleGuardDecision(
            site_count=int(site_count),
            input_protein_count=input_protein_count,
            input_kinase_count=input_kinase_count,
            selected_module_count=selected_module_count,
            candidate_module_counts_evaluated=candidate_module_counts_evaluated,
            candidate_module_count_upper_bound=candidate_module_count_upper_bound,
            clustering_engine=str(clustering_result.backend_name),
            clustering_engine_version=str(clustering_result.backend_version),
            backend_diagnostics=clustering_result.backend_diagnostics,
            tree_implementation=tree_generation_backend,
            tree_generation_backend=tree_generation_backend,
            tree_generation_mode="full_exact_tree_construction",
            tree_generation_is_approximate=False,
            tree_generation_scope="module_count_selection_and_final_assignment",
            tree_generation_guard_triggered=False,
            candidate_scoring_policy=str(config.candidate_scoring_policy),
            candidate_scoring_requested_policy=str(config.candidate_scoring_policy),
            candidate_scoring_strategy=str(clustering_result.candidate_scoring_mode),
            candidate_scoring_is_approximate=candidate_scoring_is_approximate,
            candidate_scoring_guard_triggered=False,
            candidate_scoring_sampled_site_total=_sampled_sites_total(
                sampled_candidate_scoring
            ),
            candidate_scoring_sampled_pair_count=_sampled_pairs_total(
                sampled_candidate_scoring
            ),
            max_exact_tree_sites=int(config.max_exact_tree_sites),
            max_full_candidate_scoring_sites=int(
                config.max_full_candidate_scoring_sites
            ),
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
            candidate_scoring_sampling=sampled_candidate_scoring,
            scale_guard_passed=True,
            final_module_assignment_backend=final_module_assignment_backend,
        )


def _validated_candidate_scoring_mode(value: str) -> _CandidateScoringMode:
    if value == SIGNALOME_CANDIDATE_SCORING_POLICY_FULL:
        return SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    if value == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED:
        return SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    return SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED


__all__ = ["SignalomeClusteringRunner"]
