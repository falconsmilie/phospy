"""Internal executor for the signalome workflow."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.api.results import SignalomeWorkflowResult
from phospy.errors.workflows import WorkflowBoundaryError, WorkflowStageError
from phospy.signalomes.clustering import (
    ClusterSitesResult,
    cluster_sites_with_diagnostics,
    derive_protein_modules,
)
from phospy.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
)
from phospy.signalomes.science import (
    build_expanded_signalome_table,
    build_kinase_network,
    build_module_assignments,
    build_signalome_module_table,
    select_kinase_substrates,
)
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeExecutionConfig,
    ResolvedSignalomeWorkflowRequest,
)


@dataclass(frozen=True, slots=True)
class _ExecutionMetadata:
    prediction_sites: int
    prediction_kinases: int
    downstream_score_sites: int
    downstream_score_kinases: int
    downstream_score_source: str


@dataclass(frozen=True, slots=True)
class _ClusteringStage:
    clustering_result: ClusterSitesResult
    protein_modules: pd.Series


@dataclass(frozen=True, slots=True)
class _SupportSummary:
    kinase_substrates: dict[str, tuple[str, ...]]
    support_counts: dict[str, int]


@dataclass(frozen=True, slots=True)
class _SignalomeModuleStage:
    signalome_modules: pd.DataFrame
    module_count: int
    support_summary: _SupportSummary


class SignalomeWorkflowExecutor:
    """Run signalome stage logic and assemble `SignalomeWorkflowResult`."""

    _MODULE_ID_COLUMN = "module_id"
    _NETWORK_SEAM = "signalome.executor.network"
    _KINASE_SUPPORT_SEAM = "signalome.executor.kinase_support"
    _MODULE_CONSTRUCTION_SEAM = "signalome.executor.module_construction"
    _EXPANDED_SIGNALOME_SEAM = "signalome.executor.expanded_signalome"

    def run(self, request: ResolvedSignalomeWorkflowRequest) -> SignalomeWorkflowResult:
        config = request.execution_config
        execution_metadata = self._collect_execution_metadata(request)
        clustering_stage = self._run_clustering_and_module_derivation(
            request=request,
            config=config,
            execution_metadata=execution_metadata,
        )
        module_assignments = self._build_module_assignments(
            request=request,
            config=config,
            clustering_result=clustering_stage.clustering_result,
            protein_modules=clustering_stage.protein_modules,
            execution_metadata=execution_metadata,
        )
        signalome_module_stage = self._build_signalome_modules(
            request=request,
            config=config,
            clustering_result=clustering_stage.clustering_result,
            module_assignments=module_assignments,
            execution_metadata=execution_metadata,
        )
        network_edges, network_nodes = self._build_network(
            request=request,
            config=config,
            clustering_result=clustering_stage.clustering_result,
            support_summary=signalome_module_stage.support_summary,
            execution_metadata=execution_metadata,
        )
        expanded_signalome = self._build_expanded_signalome(
            request=request,
            config=config,
            module_assignments=module_assignments,
            signalome_modules=signalome_module_stage.signalome_modules,
            network_edges=network_edges,
            support_summary=signalome_module_stage.support_summary,
            module_count=signalome_module_stage.module_count,
            execution_metadata=execution_metadata,
        )

        return self._assemble_result(
            request=request,
            clustering_result=clustering_stage.clustering_result,
            module_assignments=module_assignments,
            signalome_modules=signalome_module_stage.signalome_modules,
            network_edges=network_edges,
            network_nodes=network_nodes,
            expanded_signalome=expanded_signalome,
        )

    @staticmethod
    def _assemble_result(
        *,
        request: ResolvedSignalomeWorkflowRequest,
        clustering_result: ClusterSitesResult,
        module_assignments: pd.DataFrame,
        signalome_modules: pd.DataFrame,
        network_edges: pd.DataFrame,
        network_nodes: pd.DataFrame,
        expanded_signalome: pd.DataFrame,
    ) -> SignalomeWorkflowResult:
        return SignalomeWorkflowResult._from_owned(
            dataset=request.dataset,
            kinase_result=request.kinase_result,
            module_assignments=SignalomeAssignments._from_owned(
                table=module_assignments
            ),
            signalome_modules=SignalomeModules._from_owned(table=signalome_modules),
            kinase_network=KinaseNetwork._from_owned(
                edges=network_edges,
                nodes=network_nodes,
            ),
            module_selection_diagnostics=clustering_result.module_selection_diagnostics,
            expanded_signalome=expanded_signalome,
        )

    @staticmethod
    def _collect_execution_metadata(
        request: ResolvedSignalomeWorkflowRequest,
    ) -> _ExecutionMetadata:
        return _ExecutionMetadata(
            prediction_sites=int(request.prediction_matrix.shape[0]),
            prediction_kinases=int(request.prediction_matrix.shape[1]),
            downstream_score_sites=int(request.downstream_score_matrix.shape[0]),
            downstream_score_kinases=int(request.downstream_score_matrix.shape[1]),
            downstream_score_source=request.downstream_score_source,
        )

    def _run_clustering_and_module_derivation(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        execution_metadata: _ExecutionMetadata,
    ) -> _ClusteringStage:
        try:
            clustering_result = cluster_sites_with_diagnostics(
                scoring_matrix=request.downstream_score_matrix,
                requested_module_count=config.requested_module_count,
                primary_threshold=config.module_selection_primary_threshold,
                fallback_threshold=config.module_selection_fallback_threshold,
                max_clusters=config.module_selection_max_clusters,
            )
            protein_modules = derive_protein_modules(
                site_clusters=clustering_result.site_clusters,
                site_to_protein=request.site_to_protein,
            )
            return _ClusteringStage(
                clustering_result=clustering_result,
                protein_modules=protein_modules,
            )
        except (WorkflowStageError, ValueError) as exc:
            self._raise_boundary_error(
                seam=self._MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "ensure downstream signalome scores are compatible with module "
                    "selection clustering and module-count policy settings"
                ),
                requested_module_count=self._requested_module_count_label(
                    config.requested_module_count
                ),
                module_selection_primary_correlation_threshold=config.module_selection_primary_threshold,
                module_selection_fallback_correlation_threshold=config.module_selection_fallback_threshold,
                module_selection_max_clusters=config.module_selection_max_clusters,
                downstream_score_sites=execution_metadata.downstream_score_sites,
                downstream_score_kinases=execution_metadata.downstream_score_kinases,
                stage_error=str(exc),
            )

    def _build_module_assignments(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        clustering_result: ClusterSitesResult,
        protein_modules: pd.Series,
        execution_metadata: _ExecutionMetadata,
    ) -> pd.DataFrame:
        try:
            return build_module_assignments(
                prediction_matrix=request.prediction_matrix,
                site_to_protein=request.site_to_protein,
                protein_modules=protein_modules,
            )
        except WorkflowStageError as exc:
            self._raise_boundary_error(
                seam=self._MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "ensure interpreted prediction inputs provide unique site IDs and "
                    "resolvable site-to-protein assignments"
                ),
                **self._prediction_shape_details(execution_metadata),
                substrate_support_cutoff=config.substrate_support_cutoff,
                network_correlation_threshold=config.network_correlation_threshold,
                network_policy=config.network_policy,
                **self._module_selection_details(clustering_result),
                stage_error=str(exc),
            )

    def _build_signalome_modules(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        clustering_result: ClusterSitesResult,
        module_assignments: pd.DataFrame,
        execution_metadata: _ExecutionMetadata,
    ) -> _SignalomeModuleStage:
        support_summary = self._select_supported_substrates(
            request=request,
            config=config,
            execution_metadata=execution_metadata,
        )
        support_counts = support_summary.support_counts
        try:
            signalome_modules = build_signalome_module_table(
                module_assignments=module_assignments,
                kinase_substrates=support_summary.kinase_substrates,
                kinase_order=request.prediction_matrix.columns.astype(str).tolist(),
                assignment_policy=config.assignment_policy,
            )
        except WorkflowStageError as exc:
            self._raise_boundary_error(
                seam=self._MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "ensure supported substrates map to interpreted proteins so module "
                    "aggregation can be computed"
                ),
                **self._support_details(support_counts),
                **self._prediction_shape_details(execution_metadata),
                substrate_support_cutoff=config.substrate_support_cutoff,
                network_correlation_threshold=config.network_correlation_threshold,
                network_policy=config.network_policy,
                stage_error=str(exc),
            )
        module_count = int(
            module_assignments.loc[
                module_assignments.loc[:, self._MODULE_ID_COLUMN].astype("int64") > 0,
                self._MODULE_ID_COLUMN,
            ].nunique(dropna=True)
        )
        module_rows_with_support = (
            int((signalome_modules.sum(axis=1) > 0.0).sum())
            if not signalome_modules.empty
            else 0
        )
        if (
            signalome_modules.empty
            or module_rows_with_support == 0
            or (module_count < 2 and support_counts["supported_kinases"] < 2)
        ):
            self._raise_boundary_error(
                seam=self._MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "increase interpreted signal diversity and ensure at least two "
                    "supported kinases contribute non-zero module signal"
                ),
                module_count=module_count,
                module_rows_with_support=module_rows_with_support,
                **self._support_details(support_counts),
                **self._prediction_shape_details(execution_metadata),
                substrate_support_cutoff=config.substrate_support_cutoff,
                network_correlation_threshold=config.network_correlation_threshold,
                selected_module_count=self._module_selection_details(clustering_result)[
                    "selected_module_count"
                ],
            )
        return _SignalomeModuleStage(
            signalome_modules=signalome_modules,
            module_count=module_count,
            support_summary=support_summary,
        )

    def _select_supported_substrates(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        execution_metadata: _ExecutionMetadata,
    ) -> _SupportSummary:
        kinase_substrates = select_kinase_substrates(
            prediction_matrix=request.prediction_matrix,
            cutoff=config.substrate_support_cutoff,
        )
        support_counts = self._summarize_support(kinase_substrates)
        if support_counts["supported_kinases"] == 0:
            self._raise_boundary_error(
                seam=self._KINASE_SUPPORT_SEAM,
                next_action=(
                    "lower substrate_support_cutoff or ensure prediction scores "
                    "include kinase-site support above the configured threshold"
                ),
                **self._prediction_shape_details(execution_metadata),
                **self._support_details(support_counts),
                substrate_support_cutoff=config.substrate_support_cutoff,
            )
        return _SupportSummary(
            kinase_substrates=kinase_substrates,
            support_counts=support_counts,
        )

    def _build_network(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        clustering_result: ClusterSitesResult,
        support_summary: _SupportSummary,
        execution_metadata: _ExecutionMetadata,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        score_variance_kinases = self._score_variance_kinases(
            request.downstream_score_matrix
        )
        try:
            network_edges, network_nodes = build_kinase_network(
                downstream_score_matrix=request.downstream_score_matrix,
                kinase_order=request.prediction_matrix.columns.astype(str).tolist(),
                kinase_substrates=support_summary.kinase_substrates,
                threshold=config.network_correlation_threshold,
                network_policy=config.network_policy,
            )
        except WorkflowStageError as exc:
            self._raise_boundary_error(
                seam=self._NETWORK_SEAM,
                next_action=(
                    "ensure interpreted score and prediction matrices share the same "
                    "kinase set and contain variable score signal"
                ),
                shared_kinases=execution_metadata.prediction_kinases,
                **self._support_details(support_summary.support_counts),
                downstream_score_sites=execution_metadata.downstream_score_sites,
                downstream_score_source=execution_metadata.downstream_score_source,
                score_variance_kinases=score_variance_kinases,
                network_correlation_threshold=config.network_correlation_threshold,
                network_policy=config.network_policy,
                stage_error=str(exc),
            )
        if network_edges.empty:
            self._raise_boundary_error(
                seam=self._NETWORK_SEAM,
                next_action=(
                    "lower network_correlation_threshold or provide more variable "
                    "score profiles so kinase correlations can be estimated"
                ),
                shared_kinases=execution_metadata.prediction_kinases,
                **self._support_details(support_summary.support_counts),
                downstream_score_sites=execution_metadata.downstream_score_sites,
                downstream_score_source=execution_metadata.downstream_score_source,
                score_variance_kinases=score_variance_kinases,
                network_correlation_threshold=config.network_correlation_threshold,
                network_policy=config.network_policy,
                selected_module_count=self._module_selection_details(clustering_result)[
                    "selected_module_count"
                ],
            )
        return network_edges, network_nodes

    def _build_expanded_signalome(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        module_assignments: pd.DataFrame,
        signalome_modules: pd.DataFrame,
        network_edges: pd.DataFrame,
        support_summary: _SupportSummary,
        module_count: int,
        execution_metadata: _ExecutionMetadata,
    ) -> pd.DataFrame:
        try:
            return build_expanded_signalome_table(
                module_assignments=module_assignments,
                signalome_modules=signalome_modules,
                kinase_network_edges=network_edges,
                kinase_substrates=support_summary.kinase_substrates,
                assignment_policy=config.assignment_policy,
            )
        except WorkflowStageError as exc:
            self._raise_boundary_error(
                seam=self._EXPANDED_SIGNALOME_SEAM,
                next_action=(
                    "ensure module assignments, module signal, and network topology "
                    "are mutually consistent for expanded signalome output"
                ),
                assignment_policy=config.assignment_policy,
                module_count=module_count,
                **self._support_details(support_summary.support_counts),
                **self._prediction_shape_details(execution_metadata),
                stage_error=str(exc),
            )

    @staticmethod
    def _summarize_support(
        kinase_substrates: dict[str, tuple[str, ...]],
    ) -> dict[str, int]:
        supported_site_ids: set[str] = set()
        supported_kinases = 0
        for substrates in kinase_substrates.values():
            resolved_sites = tuple(str(site_id) for site_id in substrates)
            if not resolved_sites:
                continue
            supported_kinases += 1
            supported_site_ids.update(resolved_sites)
        return {
            "supported_sites": int(len(supported_site_ids)),
            "supported_kinases": int(supported_kinases),
        }

    @staticmethod
    def _score_variance_kinases(downstream_score_matrix: pd.DataFrame) -> int:
        if downstream_score_matrix.empty:
            return 0
        variances = downstream_score_matrix.astype(float).var(axis=0, ddof=0)
        return int((variances > 0.0).sum())

    @staticmethod
    def _prediction_shape_details(
        execution_metadata: _ExecutionMetadata,
    ) -> dict[str, int]:
        return {
            "prediction_sites": execution_metadata.prediction_sites,
            "prediction_kinases": execution_metadata.prediction_kinases,
        }

    @staticmethod
    def _support_details(support_counts: dict[str, int]) -> dict[str, int]:
        return {
            "supported_sites": int(support_counts["supported_sites"]),
            "supported_kinases": int(support_counts["supported_kinases"]),
        }

    @staticmethod
    def _module_selection_details(
        clustering_result: ClusterSitesResult,
    ) -> dict[str, int | str]:
        diagnostics = clustering_result.module_selection_diagnostics
        return {
            "selected_module_count": int(diagnostics.selected_module_count),
            "requested_module_count": SignalomeWorkflowExecutor._requested_module_count_label(
                diagnostics.requested_module_count
            ),
        }

    @staticmethod
    def _requested_module_count_label(
        requested_module_count: int | None,
    ) -> int | str:
        if requested_module_count is None:
            return "auto"
        return int(requested_module_count)

    @staticmethod
    def _raise_boundary_error(
        *,
        seam: str,
        next_action: str,
        **details: int | float | str,
    ) -> None:
        details_text = ", ".join(f"{key}={value}" for key, value in details.items())
        raise WorkflowBoundaryError(
            "signalome workflow boundary validation failed at "
            f"seam={seam}; {details_text}; next_action={next_action}"
        )
