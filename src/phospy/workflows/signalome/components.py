"""Focused execution components for the signalome workflow."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import asdict, dataclass

import pandas as pd

from phospy.errors.workflows import WorkflowBoundaryError, WorkflowStageError
from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.hashing import fingerprint_optional_table
from phospy.provenance.models import (
    PreprocessingStageProvenance,
    RunProvenance,
    TableFingerprint,
)
from phospy.provenance.serialization import to_payload as provenance_to_payload
from phospy.signalomes.clustering import (
    MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER,
    SIGNALOME_CANDIDATE_SCORING_APPLIES_TO,
    SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED,
    SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD,
    SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE,
    ClusterSitesResult,
    cluster_sites_with_diagnostics,
    derive_protein_modules,
)
from phospy.signalomes.constants import MODULE_ID_COLUMN
from phospy.signalomes.context import (
    build_protein_site_context_table,
    build_site_membership_table,
)
from phospy.signalomes.models import SignalomeNetworkCorrelationDiagnostics
from phospy.signalomes.science import (
    build_kinase_network_with_diagnostics,
    build_module_assignments,
    build_signalome_module_table,
    select_kinase_substrates,
)
from phospy.workflows.signalome.constants import (
    SIGNALOME_EXECUTOR_CONTEXT_TABLES_SEAM,
    SIGNALOME_EXECUTOR_KINASE_SUPPORT_SEAM,
    SIGNALOME_EXECUTOR_MODULE_CONSTRUCTION_SEAM,
    SIGNALOME_EXECUTOR_NETWORK_SEAM,
    SIGNALOME_WORKFLOW_BOUNDARY_MESSAGE_PREFIX,
)
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeExecutionConfig,
    ResolvedSignalomeWorkflowRequest,
)


@dataclass(frozen=True, slots=True)
class SignalomeExecutionMetadata:
    prediction_sites: int
    prediction_kinases: int
    downstream_score_sites: int
    downstream_score_kinases: int
    downstream_score_source: str


@dataclass(frozen=True, slots=True)
class SignalomeScaleGuardDecision:
    site_count: int
    cluster_tree_backend: str
    candidate_scoring_backend: str
    max_exact_cluster_tree_sites: int
    max_full_correlation_sites: int
    exact_cluster_tree_built: bool
    candidate_scoring_mode: str
    candidate_scoring_sampling: dict[str, object] | None
    scale_guard_passed: bool
    candidate_scoring_evaluated: bool = False
    candidate_scoring_skip_reason: str | None = None
    candidate_scoring_applies_to: str = SIGNALOME_CANDIDATE_SCORING_APPLIES_TO
    final_module_assignment_backend: str = (
        SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE
    )
    final_module_assignment_uses_candidate_scoring: bool = False


@dataclass(frozen=True, slots=True)
class SignalomeClusteringRunResult:
    clustering_result: ClusterSitesResult
    protein_modules: pd.Series


@dataclass(frozen=True, slots=True)
class SignalomeSupportSummary:
    kinase_substrates: dict[str, tuple[str, ...]]
    support_counts: dict[str, int]


@dataclass(frozen=True, slots=True)
class SignalomeModuleTableBuildResult:
    module_assignments: pd.DataFrame
    signalome_modules: pd.DataFrame
    module_count: int
    support_summary: SignalomeSupportSummary


@dataclass(frozen=True, slots=True)
class SignalomeNetworkBuildResult:
    edges: pd.DataFrame
    nodes: pd.DataFrame
    candidate_correlations: pd.DataFrame
    correlation_diagnostics: SignalomeNetworkCorrelationDiagnostics


@dataclass(frozen=True, slots=True)
class SignalomeContextTableBuildResult:
    site_membership: pd.DataFrame
    protein_site_context: pd.DataFrame


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
            _raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "ensure downstream signalome scores are compatible with module "
                    "selection clustering and module-count policy settings"
                ),
                requested_module_count=_requested_module_count_label(
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
        candidate_scoring_sampling = clustering_result.candidate_scoring_sampling
        if (
            candidate_scoring_sampling is None
            and str(config.candidate_scoring_backend)
            == SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED
        ):
            candidate_scoring_sampling = {
                "sampling_cap": int(MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER),
                "sampling_method": SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD,
                "deterministic_seed_policy": (
                    SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY
                ),
                "actual_sampled_pair_count": 0,
                "per_cluster_sample_count_summary": {
                    "min": 0,
                    "max": 0,
                    "mean": 0.0,
                    "total": 0,
                },
            }
        return SignalomeScaleGuardDecision(
            site_count=int(site_count),
            cluster_tree_backend=str(config.cluster_tree_backend),
            candidate_scoring_backend=str(config.candidate_scoring_backend),
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
            candidate_scoring_sampling=candidate_scoring_sampling,
            candidate_scoring_applies_to=SIGNALOME_CANDIDATE_SCORING_APPLIES_TO,
            final_module_assignment_backend=(
                SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE
            ),
            final_module_assignment_uses_candidate_scoring=False,
            scale_guard_passed=True,
        )


class SignalomeModuleTableBuilder:
    """Build module-assignment and module summary tables for signalome execution."""

    _MODULE_ID_COLUMN = MODULE_ID_COLUMN

    def __init__(
        self,
        *,
        build_assignments: Callable[..., pd.DataFrame] = build_module_assignments,
        select_substrates: Callable[..., dict[str, tuple[str, ...]]] = (
            select_kinase_substrates
        ),
        build_modules: Callable[..., pd.DataFrame] = build_signalome_module_table,
    ) -> None:
        self._build_assignments = build_assignments
        self._select_substrates = select_substrates
        self._build_modules = build_modules

    def run(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        clustering_result: ClusterSitesResult,
        protein_modules: pd.Series,
        execution_metadata: SignalomeExecutionMetadata,
    ) -> SignalomeModuleTableBuildResult:
        module_assignments = self._build_module_assignments(
            request=request,
            config=config,
            clustering_result=clustering_result,
            protein_modules=protein_modules,
            execution_metadata=execution_metadata,
        )
        support_summary = self._select_supported_substrates(
            request=request,
            config=config,
            execution_metadata=execution_metadata,
        )
        support_counts = support_summary.support_counts
        try:
            signalome_modules = self._build_modules(
                module_assignments=module_assignments,
                kinase_substrates=support_summary.kinase_substrates,
                kinase_order=request.prediction_matrix.columns.astype(str).tolist(),
                assignment_policy=config.assignment_policy,
            )
        except WorkflowStageError as exc:
            _raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "ensure supported substrates map to interpreted proteins so module "
                    "aggregation can be computed"
                ),
                **_support_details(support_counts),
                **_prediction_shape_details(execution_metadata),
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
            _raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "increase interpreted signal diversity and ensure at least two "
                    "supported kinases contribute non-zero module signal"
                ),
                module_count=module_count,
                module_rows_with_support=module_rows_with_support,
                **_support_details(support_counts),
                **_prediction_shape_details(execution_metadata),
                substrate_support_cutoff=config.substrate_support_cutoff,
                network_correlation_threshold=config.network_correlation_threshold,
                selected_module_count=_module_selection_details(clustering_result)[
                    "selected_module_count"
                ],
            )
        return SignalomeModuleTableBuildResult(
            module_assignments=module_assignments,
            signalome_modules=signalome_modules,
            module_count=module_count,
            support_summary=support_summary,
        )

    def _build_module_assignments(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        clustering_result: ClusterSitesResult,
        protein_modules: pd.Series,
        execution_metadata: SignalomeExecutionMetadata,
    ) -> pd.DataFrame:
        try:
            return self._build_assignments(
                prediction_matrix=request.prediction_matrix,
                site_to_protein=request.site_to_protein,
                protein_modules=protein_modules,
            )
        except WorkflowStageError as exc:
            _raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "ensure interpreted prediction inputs provide unique site IDs and "
                    "resolvable site-to-protein assignments"
                ),
                **_prediction_shape_details(execution_metadata),
                substrate_support_cutoff=config.substrate_support_cutoff,
                network_correlation_threshold=config.network_correlation_threshold,
                network_policy=config.network_policy,
                **_module_selection_details(clustering_result),
                stage_error=str(exc),
            )

    def _select_supported_substrates(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        execution_metadata: SignalomeExecutionMetadata,
    ) -> SignalomeSupportSummary:
        kinase_substrates = self._select_substrates(
            prediction_matrix=request.prediction_matrix,
            cutoff=config.substrate_support_cutoff,
        )
        support_counts = _summarize_support(kinase_substrates)
        if support_counts["supported_kinases"] == 0:
            _raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_KINASE_SUPPORT_SEAM,
                next_action=(
                    "lower substrate_support_cutoff or ensure prediction scores "
                    "include kinase-site support above the configured threshold"
                ),
                **_prediction_shape_details(execution_metadata),
                **_support_details(support_counts),
                substrate_support_cutoff=config.substrate_support_cutoff,
            )
        return SignalomeSupportSummary(
            kinase_substrates=kinase_substrates,
            support_counts=support_counts,
        )


class SignalomeNetworkBuilder:
    """Build kinase network outputs and diagnostics for signalome execution."""

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
        score_variance_kinases = _score_variance_kinases(
            request.downstream_score_matrix
        )
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
            )
        except WorkflowStageError as exc:
            _raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_NETWORK_SEAM,
                next_action=(
                    "ensure interpreted score and prediction matrices share the same "
                    "kinase set and contain variable score signal"
                ),
                shared_kinases=execution_metadata.prediction_kinases,
                **_support_details(support_summary.support_counts),
                downstream_score_sites=execution_metadata.downstream_score_sites,
                downstream_score_source=execution_metadata.downstream_score_source,
                score_variance_kinases=score_variance_kinases,
                network_correlation_threshold=config.network_correlation_threshold,
                network_policy=config.network_policy,
                stage_error=str(exc),
            )
        if network_edges.empty:
            _raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_NETWORK_SEAM,
                next_action=(
                    "lower network_correlation_threshold or provide more variable "
                    "score profiles so kinase correlations can be estimated"
                ),
                shared_kinases=execution_metadata.prediction_kinases,
                **_support_details(support_summary.support_counts),
                downstream_score_sites=execution_metadata.downstream_score_sites,
                downstream_score_source=execution_metadata.downstream_score_source,
                score_variance_kinases=score_variance_kinases,
                network_correlation_threshold=config.network_correlation_threshold,
                network_policy=config.network_policy,
                selected_module_count=_module_selection_details(clustering_result)[
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
            )
        return SignalomeNetworkBuildResult(
            edges=network_edges,
            nodes=network_nodes,
            candidate_correlations=candidate_correlations,
            correlation_diagnostics=correlation_diagnostics,
        )


class SignalomeContextTableBuilder:
    """Build site and protein context tables for signalome execution."""

    def __init__(
        self,
        *,
        build_site_membership: Callable[
            ..., pd.DataFrame
        ] = build_site_membership_table,
        build_protein_context: Callable[..., pd.DataFrame] = (
            build_protein_site_context_table
        ),
    ) -> None:
        self._build_site_membership = build_site_membership
        self._build_protein_context = build_protein_context

    def run(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        clustering_result: ClusterSitesResult,
        module_assignments: pd.DataFrame,
        support_summary: SignalomeSupportSummary,
        execution_metadata: SignalomeExecutionMetadata,
    ) -> SignalomeContextTableBuildResult:
        try:
            site_membership = self._build_site_membership(
                module_assignments=module_assignments,
                site_clusters=clustering_result.site_clusters,
                site_metadata=request.dataset.site_metadata,
                prediction_matrix=request.prediction_matrix,
                kinase_substrates=support_summary.kinase_substrates,
                substrate_support_cutoff=config.substrate_support_cutoff,
                assignment_policy=config.assignment_policy,
            )
            protein_site_context = self._build_protein_context(
                site_membership=site_membership
            )
            return SignalomeContextTableBuildResult(
                site_membership=site_membership,
                protein_site_context=protein_site_context,
            )
        except (WorkflowStageError, ValueError, TypeError, KeyError) as exc:
            _raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_CONTEXT_TABLES_SEAM,
                next_action=(
                    "ensure module assignments, site clusters, site metadata, "
                    "prediction scores, and supported substrate mappings are "
                    "mutually consistent before building signalome context tables"
                ),
                **_prediction_shape_details(execution_metadata),
                module_assignment_rows=int(module_assignments.shape[0]),
                module_assignment_columns=int(module_assignments.shape[1]),
                site_cluster_count=int(clustering_result.site_clusters.nunique()),
                **_support_details(support_summary.support_counts),
                substrate_support_cutoff=config.substrate_support_cutoff,
                assignment_policy=config.assignment_policy,
                stage_error=str(exc),
            )


class SignalomeProvenanceBuilder:
    """Build workflow-level provenance for signalome execution."""

    def __init__(
        self,
        *,
        collect_environment: Callable[[], object] = collect_environment_provenance,
    ) -> None:
        self._collect_environment = collect_environment

    def build(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        clustering_result: ClusterSitesResult,
        module_assignments: pd.DataFrame,
        signalome_modules: pd.DataFrame,
        network_edges: pd.DataFrame,
        network_nodes: pd.DataFrame,
        candidate_correlations: pd.DataFrame,
        network_correlation_diagnostics: SignalomeNetworkCorrelationDiagnostics,
        expanded_signalome: pd.DataFrame,
        site_membership: pd.DataFrame,
        protein_site_context: pd.DataFrame,
        scale_guard_decision: SignalomeScaleGuardDecision,
    ) -> RunProvenance:
        input_tables = self._collect_fingerprints(
            (
                ("dataset.phospho", request.dataset.phospho),
                ("dataset.site_metadata", request.dataset.site_metadata),
                ("dataset.sample_metadata", request.dataset.sample_metadata),
                ("dataset.total", request.dataset.total),
                ("dataset.comparisons", request.dataset.comparisons),
                ("upstream.prediction.pred_mat", request.prediction_matrix),
                (
                    "upstream.scoring.downstream_score_matrix",
                    request.downstream_score_matrix,
                ),
            )
        )
        output_tables = self._collect_fingerprints(
            (
                ("outputs.signalome.module_assignments", module_assignments),
                ("outputs.signalome.signalome_modules", signalome_modules),
                ("outputs.signalome.kinase_network.edges", network_edges),
                ("outputs.signalome.kinase_network.nodes", network_nodes),
                (
                    "outputs.signalome.kinase_network.candidate_correlations",
                    candidate_correlations,
                ),
                ("outputs.signalome.expanded_signalome", expanded_signalome),
                ("outputs.signalome.site_membership", site_membership),
                ("outputs.signalome.protein_site_context", protein_site_context),
            )
        )
        upstream_provenance = request.kinase_result.provenance
        return RunProvenance(
            environment=self._collect_environment(),
            input_tables=input_tables,
            preprocessing_stages=self._dataset_preprocessing_stages(request),
            reference=request.kinase_result.references.provenance,
            workflow_name="signalome_workflow",
            workflow_parameters={
                "signalome_config": {
                    "substrate_support_cutoff": float(config.substrate_support_cutoff),
                    "network_correlation_threshold": float(
                        config.network_correlation_threshold
                    ),
                    "network_policy": str(config.network_policy),
                    "assignment_policy": str(config.assignment_policy),
                    "score_preconditioning_policy": str(
                        config.score_preconditioning_policy
                    ),
                    "module_selection_primary_correlation_threshold": float(
                        config.module_selection_primary_threshold
                    ),
                    "module_selection_fallback_correlation_threshold": float(
                        config.module_selection_fallback_threshold
                    ),
                    "module_selection_max_clusters": int(
                        config.module_selection_max_clusters
                    ),
                    "cluster_tree_backend": str(config.cluster_tree_backend),
                    "candidate_scoring_backend": str(config.candidate_scoring_backend),
                    "max_exact_cluster_tree_sites": int(
                        config.max_exact_cluster_tree_sites
                    ),
                    "max_full_correlation_sites": int(
                        config.max_full_correlation_sites
                    ),
                    "deprecated_clustering_backend_alias": str(
                        config.clustering_backend_alias
                    ),
                    "deprecated_max_exact_clustering_sites_alias": int(
                        config.max_exact_cluster_tree_sites
                    ),
                    "module_count": (
                        None
                        if config.requested_module_count is None
                        else int(config.requested_module_count)
                    ),
                },
                "scale_guard": {
                    "site_count": int(scale_guard_decision.site_count),
                    "cluster_tree_backend": str(
                        scale_guard_decision.cluster_tree_backend
                    ),
                    "candidate_scoring_backend": str(
                        scale_guard_decision.candidate_scoring_backend
                    ),
                    "max_exact_cluster_tree_sites": int(
                        scale_guard_decision.max_exact_cluster_tree_sites
                    ),
                    "max_full_correlation_sites": int(
                        scale_guard_decision.max_full_correlation_sites
                    ),
                    "exact_cluster_tree_built": bool(
                        scale_guard_decision.exact_cluster_tree_built
                    ),
                    "candidate_scoring_mode": str(
                        scale_guard_decision.candidate_scoring_mode
                    ),
                    "candidate_scoring_evaluated": bool(
                        scale_guard_decision.candidate_scoring_evaluated
                    ),
                    "candidate_scoring_skip_reason": (
                        None
                        if scale_guard_decision.candidate_scoring_skip_reason is None
                        else str(scale_guard_decision.candidate_scoring_skip_reason)
                    ),
                    "candidate_scoring_sampling": (
                        None
                        if scale_guard_decision.candidate_scoring_sampling is None
                        else dict(scale_guard_decision.candidate_scoring_sampling)
                    ),
                    "candidate_scoring_applies_to": str(
                        scale_guard_decision.candidate_scoring_applies_to
                    ),
                    "final_module_assignment_backend": str(
                        scale_guard_decision.final_module_assignment_backend
                    ),
                    "final_module_assignment_uses_candidate_scoring": bool(
                        scale_guard_decision.final_module_assignment_uses_candidate_scoring
                    ),
                    "scale_guard_passed": bool(scale_guard_decision.scale_guard_passed),
                },
                "module_selection_diagnostics": asdict(
                    clustering_result.module_selection_diagnostics
                ),
                "score_preconditioning_diagnostics": asdict(
                    request.score_preconditioning_diagnostics
                ),
                "network_correlation_diagnostics": asdict(
                    network_correlation_diagnostics
                ),
                "upstream_kinase_provenance": (
                    None
                    if upstream_provenance is None
                    else provenance_to_payload(upstream_provenance)
                ),
            },
            random_state=None,
            random_seed_policy=None,
            output_tables=output_tables,
        )

    @staticmethod
    def _dataset_preprocessing_stages(
        request: ResolvedSignalomeWorkflowRequest,
    ) -> tuple[PreprocessingStageProvenance, ...]:
        provenance = request.dataset.provenance
        if provenance is None:
            return ()
        return tuple(provenance.preprocessing_stages)

    @staticmethod
    def _collect_fingerprints(
        entries: tuple[tuple[str, pd.DataFrame | None], ...],
    ) -> tuple[TableFingerprint, ...]:
        fingerprints: list[TableFingerprint] = []
        for name, table in entries:
            fingerprint = fingerprint_optional_table(table, name=name)
            if fingerprint is None:
                continue
            fingerprints.append(fingerprint)
        return tuple(fingerprints)


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


def _score_variance_kinases(downstream_score_matrix: pd.DataFrame) -> int:
    if downstream_score_matrix.empty:
        return 0
    variances = downstream_score_matrix.astype(float).var(axis=0, ddof=0)
    return int((variances > 0.0).sum())


def _prediction_shape_details(
    execution_metadata: SignalomeExecutionMetadata,
) -> dict[str, int]:
    return {
        "prediction_sites": execution_metadata.prediction_sites,
        "prediction_kinases": execution_metadata.prediction_kinases,
    }


def _support_details(support_counts: dict[str, int]) -> dict[str, int]:
    return {
        "supported_sites": int(support_counts["supported_sites"]),
        "supported_kinases": int(support_counts["supported_kinases"]),
    }


def _module_selection_details(
    clustering_result: ClusterSitesResult,
) -> dict[str, int | str]:
    diagnostics = clustering_result.module_selection_diagnostics
    return {
        "selected_module_count": int(diagnostics.selected_module_count),
        "requested_module_count": _requested_module_count_label(
            diagnostics.requested_module_count
        ),
    }


def _requested_module_count_label(requested_module_count: int | None) -> int | str:
    if requested_module_count is None:
        return "auto"
    return int(requested_module_count)


def _raise_boundary_error(
    *,
    seam: str,
    next_action: str,
    **details: object,
) -> None:
    raise WorkflowBoundaryError(
        seam=seam,
        next_action=next_action,
        details=details,
        message_prefix=SIGNALOME_WORKFLOW_BOUNDARY_MESSAGE_PREFIX,
    )


__all__ = [
    "SignalomeClusteringRunResult",
    "SignalomeClusteringRunner",
    "SignalomeContextTableBuildResult",
    "SignalomeContextTableBuilder",
    "SignalomeExecutionMetadata",
    "SignalomeModuleTableBuildResult",
    "SignalomeModuleTableBuilder",
    "SignalomeNetworkBuildResult",
    "SignalomeNetworkBuilder",
    "SignalomeProvenanceBuilder",
    "SignalomeScaleGuardDecision",
    "SignalomeSupportSummary",
]
