"""Signalome provenance and diagnostics assembly."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import asdict

import pandas as pd

from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.hashing import fingerprint_optional_table
from phospy.provenance.models import (
    PreprocessingStageProvenance,
    RunProvenance,
    TableFingerprint,
)
from phospy.provenance.serialization import to_payload as provenance_to_payload
from phospy.signalomes.clustering import ClusterSitesResult
from phospy.signalomes.models import SignalomeNetworkCorrelationDiagnostics
from phospy.workflows.signalome.component_models import SignalomeScaleGuardDecision
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeExecutionConfig,
    ResolvedSignalomeWorkflowRequest,
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
                    "clustering_backend": str(config.clustering_backend),
                    "module_count": (
                        None
                        if config.requested_module_count is None
                        else int(config.requested_module_count)
                    ),
                },
                "scale_guard": {
                    "site_count": int(scale_guard_decision.site_count),
                    "selected_module_count": int(
                        scale_guard_decision.selected_module_count
                    ),
                    "clustering_backend": str(scale_guard_decision.clustering_backend),
                    "clustering_backend_version": str(
                        scale_guard_decision.clustering_backend_version
                    ),
                    "backend_diagnostics": (
                        None
                        if scale_guard_decision.backend_diagnostics is None
                        else dict(scale_guard_decision.backend_diagnostics)
                    ),
                    "cluster_tree_backend": str(
                        scale_guard_decision.cluster_tree_backend
                    ),
                    "candidate_scoring_backend": str(
                        scale_guard_decision.candidate_scoring_backend
                    ),
                    "candidate_scoring_requested_backend": str(
                        scale_guard_decision.candidate_scoring_requested_backend
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


__all__ = ["SignalomeProvenanceBuilder"]
