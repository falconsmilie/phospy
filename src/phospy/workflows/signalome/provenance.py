"""Signalome provenance and diagnostics assembly."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import asdict

import pandas as pd

from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.hashing import fingerprint_optional_table
from phospy.provenance.models import (
    EnvironmentProvenance,
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
        collect_environment: Callable[
            [], EnvironmentProvenance
        ] = collect_environment_provenance,
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
                    "tree_engine": str(config.tree_engine),
                    "candidate_scoring_policy": str(config.candidate_scoring_policy),
                    "max_exact_tree_sites": int(config.max_exact_tree_sites),
                    "max_full_candidate_scoring_sites": int(
                        config.max_full_candidate_scoring_sites
                    ),
                    "clustering_engine": str(config.clustering_engine),
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
                    "clustering_engine": str(scale_guard_decision.clustering_engine),
                    "clustering_engine_version": str(
                        scale_guard_decision.clustering_engine_version
                    ),
                    "backend_diagnostics": (
                        None
                        if scale_guard_decision.backend_diagnostics is None
                        else dict(scale_guard_decision.backend_diagnostics)
                    ),
                    "tree_engine": str(scale_guard_decision.tree_engine),
                    "candidate_scoring_policy": str(
                        scale_guard_decision.candidate_scoring_policy
                    ),
                    "candidate_scoring_requested_policy": str(
                        scale_guard_decision.candidate_scoring_requested_policy
                    ),
                    "max_exact_tree_sites": int(
                        scale_guard_decision.max_exact_tree_sites
                    ),
                    "max_full_candidate_scoring_sites": int(
                        scale_guard_decision.max_full_candidate_scoring_sites
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
                "signalome_score_semantics": _build_signalome_score_semantics(
                    request=request,
                    config=config,
                    clustering_result=clustering_result,
                    network_correlation_diagnostics=network_correlation_diagnostics,
                    scale_guard_decision=scale_guard_decision,
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


def _build_signalome_score_semantics(
    *,
    request: ResolvedSignalomeWorkflowRequest,
    config: ResolvedSignalomeExecutionConfig,
    clustering_result: ClusterSitesResult,
    network_correlation_diagnostics: SignalomeNetworkCorrelationDiagnostics,
    scale_guard_decision: SignalomeScaleGuardDecision,
) -> dict[str, object]:
    downstream_score_source = str(request.downstream_score_source)
    module_selection_diagnostics = clustering_result.module_selection_diagnostics
    preconditioning = request.score_preconditioning_diagnostics
    return {
        "downstream_score_source": downstream_score_source,
        "downstream_score_meaning": _downstream_score_meaning(
            downstream_score_source=downstream_score_source
        ),
        "module_selection_score_meaning": (
            "module-count selection uses within-cluster correlation summaries over "
            "downstream score profiles across candidate module counts"
        ),
        "candidate_scoring_mode": str(scale_guard_decision.candidate_scoring_mode),
        "candidate_scoring_scope": str(
            scale_guard_decision.candidate_scoring_applies_to
        ),
        "network_correlation_meaning": (
            "network candidate and edge scores are pairwise correlations between "
            "downstream kinase score profiles, computed on finite paired observations"
        ),
        "network_policy": str(config.network_policy),
        "negative_correlation_handling": _negative_correlation_handling(
            network_policy=str(config.network_policy),
            network_correlation_threshold=float(config.network_correlation_threshold),
        ),
        "missing_profile_handling": {
            "all_missing_rows_before_execution": {
                "policy": str(preconditioning.policy),
                "dropped_row_count": int(preconditioning.dropped_all_missing_row_count),
                "retained_row_count": int(preconditioning.retained_row_count),
            },
            "correlation_calculation": (
                "partially missing rows are retained for pairwise-complete "
                "correlation; pairs with missing/non-finite issues are labeled via "
                "correlation_status and excluded from finite-edge creation"
            ),
            "network_missing_value_correlations": int(
                network_correlation_diagnostics.missing_value_correlations
            ),
            "network_non_finite_value_correlations": int(
                network_correlation_diagnostics.non_finite_value_correlations
            ),
            "network_insufficient_observation_correlations": int(
                network_correlation_diagnostics.insufficient_observation_correlations
            ),
        },
        "constant_profile_handling": {
            "module_selection_zero_variance_profile_count": int(
                module_selection_diagnostics.zero_variance_profile_count
            ),
            "module_selection_near_constant_profile_count": int(
                module_selection_diagnostics.near_constant_profile_count
            ),
            "module_selection_excluded_from_correlation_count": int(
                module_selection_diagnostics.excluded_from_correlation_count
            ),
            "network_constant_profile_correlations": int(
                network_correlation_diagnostics.constant_profile_correlations
            ),
            "meaning": (
                "constant/near-constant profiles have weak or undefined "
                "correlation signal and are tracked explicitly in diagnostics"
            ),
        },
        "thresholds_and_limits": {
            "substrate_support_cutoff": float(config.substrate_support_cutoff),
            "module_selection_primary_correlation_threshold": float(
                config.module_selection_primary_threshold
            ),
            "module_selection_fallback_correlation_threshold": float(
                config.module_selection_fallback_threshold
            ),
            "module_selection_max_clusters": int(config.module_selection_max_clusters),
            "network_correlation_threshold": float(
                config.network_correlation_threshold
            ),
            "max_exact_tree_sites": int(config.max_exact_tree_sites),
            "max_full_candidate_scoring_sites": int(
                config.max_full_candidate_scoring_sites
            ),
        },
        "clustering_engine": str(scale_guard_decision.clustering_engine),
        "scientific_interpretation_limits": (
            "signalome module assignments, module scores, and kinase-network "
            "correlations are derived summary statistics for this dataset and "
            "configuration; they are not probabilities, calibrated confidence "
            "values, or direct evidence of causal biological regulation"
        ),
    }


def _downstream_score_meaning(*, downstream_score_source: str) -> str:
    if downstream_score_source == "rank_weighted_fusion_scores":
        return (
            "rank-weighted fusion of upstream downstream-score lanes; larger values "
            "indicate stronger relative downstream support within the run"
        )
    if downstream_score_source == "profile_scores":
        return (
            "upstream downstream profile-score lane; larger values indicate stronger "
            "relative downstream support within the run"
        )
    return (
        "upstream downstream score lane used for signalome construction; values are "
        "relative support scores within the run"
    )


def _negative_correlation_handling(
    *,
    network_policy: str,
    network_correlation_threshold: float,
) -> str:
    if network_policy == "positive_only":
        return (
            "negative correlations are never included as edges; only correlations "
            f">= {network_correlation_threshold} are retained"
        )
    if network_policy == "absolute_threshold":
        return (
            "negative correlations can pass eligibility via absolute magnitude; "
            "edge correlation values are stored as unsigned absolute magnitudes"
        )
    if network_policy == "signed":
        return (
            "negative correlations can pass eligibility when absolute magnitude "
            "meets threshold; edge correlation values retain sign"
        )
    return "negative-correlation handling depends on configured network_policy"


__all__ = ["SignalomeProvenanceBuilder"]
