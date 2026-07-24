"""Kinase network construction for signalome workflow execution."""

from __future__ import annotations

from collections.abc import Callable

import pandas as pd

from phospy.errors.workflows import WorkflowStageError
from phospy.science.signalomes.clustering import ClusterSitesResult
from phospy.science.signalomes.models import SignalomeNetworkCorrelationDiagnostics
from phospy.science.signalomes.science import build_kinase_network_with_diagnostics
from phospy.workflows.signalome.component_helpers import (
    module_selection_details,
    raise_boundary_error,
    score_variance_kinases,
    support_details,
)
from phospy.workflows.signalome.component_models import (
    SignalomeExecutionMetadata,
    SignalomeNetworkBuildResult,
    SignalomeSupportSummary,
)
from phospy.workflows.signalome.constants import SIGNALOME_EXECUTOR_NETWORK_SEAM
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeExecutionConfig,
    ResolvedSignalomeWorkflowRequest,
)


class SignalomeNetworkBuilder:
    """Build kinase score-profile association outputs and diagnostics."""

    def __init__(
        self,
        *,
        build_network: Callable[
            ...,
            tuple[
                pd.DataFrame,
                pd.DataFrame,
                pd.DataFrame,
                SignalomeNetworkCorrelationDiagnostics,
            ],
        ] = build_kinase_network_with_diagnostics,
    ) -> None:
        self._build_network = build_network

    def run(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        clustering_result: ClusterSitesResult,
        support_summary: SignalomeSupportSummary,
        execution_metadata: SignalomeExecutionMetadata,
    ) -> SignalomeNetworkBuildResult:
        score_variance_count = score_variance_kinases(request.downstream_score_matrix)
        try:
            (
                network_edges,
                network_nodes,
                candidate_correlations,
                correlation_diagnostics,
            ) = self._build_network(
                downstream_score_matrix=request.downstream_score_matrix,
                kinase_order=request.prediction_matrix.columns.astype(str).tolist(),
                kinase_substrates=support_summary.kinase_substrates,
                threshold=config.network_correlation_threshold,
                network_policy=config.network_policy,
                min_paired_observations=config.network_min_paired_finite_observations,
            )
        except WorkflowStageError as exc:
            raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_NETWORK_SEAM,
                next_action=(
                    "ensure interpreted score and prediction matrices share the same "
                    "kinase set and contain variable score profiles"
                ),
                shared_kinases=execution_metadata.prediction_kinases,
                **support_details(support_summary.support_counts),
                downstream_score_sites=execution_metadata.downstream_score_sites,
                downstream_score_source=execution_metadata.downstream_score_source,
                score_variance_kinases=score_variance_count,
                network_correlation_threshold=config.network_correlation_threshold,
                network_policy=config.network_policy,
                network_min_paired_finite_observations=(
                    config.network_min_paired_finite_observations
                ),
                stage_error=str(exc),
            )
        if network_edges.empty and not _empty_network_allowed_by_observation_policy(
            correlation_diagnostics
        ):
            next_action = (
                "lower network_correlation_threshold or provide more variable "
                "score profiles so kinase score-profile correlations can be estimated"
            )
            if correlation_diagnostics.edges_skipped_insufficient_paired_observations:
                next_action = (
                    "provide at least "
                    f"{int(config.network_min_paired_finite_observations)} paired "
                    "finite observations for candidate kinase score-profile "
                    "correlations, or explicitly set "
                    "config.output.network_min_paired_finite_observations to a "
                    "scientifically justified value of at least 3"
                )
            raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_NETWORK_SEAM,
                next_action=next_action,
                shared_kinases=execution_metadata.prediction_kinases,
                **support_details(support_summary.support_counts),
                downstream_score_sites=execution_metadata.downstream_score_sites,
                downstream_score_source=execution_metadata.downstream_score_source,
                score_variance_kinases=score_variance_count,
                network_correlation_threshold=config.network_correlation_threshold,
                network_policy=config.network_policy,
                network_min_paired_finite_observations=(
                    config.network_min_paired_finite_observations
                ),
                selected_module_count=module_selection_details(clustering_result)[
                    "selected_module_count"
                ],
                total_candidate_correlations=(
                    correlation_diagnostics.total_candidate_correlations
                ),
                undefined_correlations=correlation_diagnostics.undefined_correlations,
                constant_profile_correlations=(
                    correlation_diagnostics.constant_profile_correlations
                ),
                insufficient_observation_correlations=(
                    correlation_diagnostics.insufficient_observation_correlations
                ),
                missing_value_correlations=(
                    correlation_diagnostics.missing_value_correlations
                ),
                non_finite_value_correlations=(
                    correlation_diagnostics.non_finite_value_correlations
                ),
                edges_skipped_below_threshold=(
                    correlation_diagnostics.edges_skipped_below_threshold
                ),
                edges_skipped_insufficient_paired_observations=(
                    correlation_diagnostics.edges_skipped_insufficient_paired_observations
                ),
                edges_skipped_constant_profile=(
                    correlation_diagnostics.edges_skipped_constant_profile
                ),
                edges_skipped_missing_score=(
                    correlation_diagnostics.edges_skipped_missing_score
                ),
                edges_skipped_non_finite_score=(
                    correlation_diagnostics.edges_skipped_non_finite_score
                ),
                edges_skipped_undefined_correlation=(
                    correlation_diagnostics.edges_skipped_undefined_correlation
                ),
            )
        return SignalomeNetworkBuildResult(
            edges=network_edges,
            nodes=network_nodes,
            candidate_correlations=candidate_correlations,
            correlation_diagnostics=correlation_diagnostics,
        )


def _empty_network_allowed_by_observation_policy(
    diagnostics: SignalomeNetworkCorrelationDiagnostics,
) -> bool:
    """Allow zero accepted edges only when the minimum-observation policy owns it."""

    total_candidates = int(diagnostics.total_candidate_correlations)
    if total_candidates <= 0:
        return False
    return (
        int(diagnostics.edges_created) == 0
        and int(diagnostics.edges_skipped_insufficient_paired_observations)
        == total_candidates
    )


__all__ = ["SignalomeNetworkBuilder"]
