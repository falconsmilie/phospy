"""Signalome provenance and diagnostics assembly."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import asdict

import pandas as pd

from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.hashing import _fingerprint_optional_table_with_normalized_axes
from phospy.provenance.models import (
    EnvironmentProvenance,
    PreprocessingStageProvenance,
    RunProvenance,
    TableFingerprint,
)
from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)
from phospy.provenance.serialization import to_payload as provenance_to_payload
from phospy.science.scoring.policy_models import DownstreamScoreSource
from phospy.science.signalomes.clustering import ClusterSitesResult
from phospy.science.signalomes.clustering.diagnostic_schemas import (
    backend_diagnostics_to_payload,
    candidate_scoring_sampling_diagnostics_to_payload,
)
from phospy.science.signalomes.clustering.policies import (
    SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_APPLIES_TO,
    SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_IMPUTED_VALUES_EXPOSED_IN_OUTPUT_TABLES,
)
from phospy.science.signalomes.clustering.scientific_policies import (
    PROTEIN_MODULE_FROM_SITE_MEMBERSHIP_POLICY,
    SignalomeMissingValueClusteringPolicy,
    build_signalome_module_candidate_score_policy,
)
from phospy.science.signalomes.clustering.tree_building import (
    SignalomeClusteringMissingValueDiagnostics,
    summarize_clustering_missing_value_diagnostics,
)
from phospy.science.signalomes.models import SignalomeNetworkCorrelationDiagnostics
from phospy.workflows.intensity_scale_evidence import (
    input_intensity_scale_evidence_payload,
)
from phospy.workflows.signalome.component_models import SignalomeScaleGuardDecision
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeExecutionConfig,
    ResolvedSignalomeWorkflowRequest,
)
from phospy.workflows.signalome.row_attrition import (
    build_signalome_row_attrition_provenance,
)
from phospy.workflows.signalome.scientific_policies import ScorePreconditioningPolicy


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

    def run(
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
        input_tables = _build_input_table_fingerprints(request)
        output_tables = _build_output_table_fingerprints(
            module_assignments=module_assignments,
            signalome_modules=signalome_modules,
            network_edges=network_edges,
            network_nodes=network_nodes,
            candidate_correlations=candidate_correlations,
            expanded_signalome=expanded_signalome,
            site_membership=site_membership,
            protein_site_context=protein_site_context,
        )
        upstream_provenance = request.kinase_result.provenance
        clustering_missing_value_diagnostics = (
            summarize_clustering_missing_value_diagnostics(
                request.downstream_score_matrix.to_numpy(dtype=float, copy=False)
            )
        )
        scientific_policies = _build_scientific_policy_records(
            request=request,
            config=config,
            clustering_result=clustering_result,
            scale_guard_decision=scale_guard_decision,
            clustering_missing_value_diagnostics=clustering_missing_value_diagnostics,
        )
        workflow_parameters = _build_workflow_parameters(
            request=request,
            config=config,
            clustering_result=clustering_result,
            network_correlation_diagnostics=network_correlation_diagnostics,
            scale_guard_decision=scale_guard_decision,
            clustering_missing_value_diagnostics=clustering_missing_value_diagnostics,
            upstream_provenance=upstream_provenance,
        )
        return RunProvenance(
            environment=self._collect_environment(),
            input_tables=input_tables,
            preprocessing_stages=self._dataset_preprocessing_stages(request),
            reference=request.kinase_result.references.provenance,
            workflow_name="signalome_workflow",
            workflow_parameters=workflow_parameters,
            random_state=None,
            random_seed_policy=None,
            output_tables=output_tables,
            scientific_policies=scientific_policies,
            reference_context=request.dataset.reference_context,
        )

    @staticmethod
    def _dataset_preprocessing_stages(
        request: ResolvedSignalomeWorkflowRequest,
    ) -> tuple[PreprocessingStageProvenance, ...]:
        provenance = request.dataset.provenance
        if provenance is None:
            return ()
        return tuple(provenance.preprocessing_stages)


def _build_input_table_fingerprints(
    request: ResolvedSignalomeWorkflowRequest,
) -> tuple[TableFingerprint, ...]:
    return _collect_fingerprints(
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


def _build_output_table_fingerprints(
    *,
    module_assignments: pd.DataFrame,
    signalome_modules: pd.DataFrame,
    network_edges: pd.DataFrame,
    network_nodes: pd.DataFrame,
    candidate_correlations: pd.DataFrame,
    expanded_signalome: pd.DataFrame,
    site_membership: pd.DataFrame,
    protein_site_context: pd.DataFrame,
) -> tuple[TableFingerprint, ...]:
    return _collect_fingerprints(
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


def _collect_fingerprints(
    entries: tuple[tuple[str, pd.DataFrame | None], ...],
) -> tuple[TableFingerprint, ...]:
    fingerprints: list[TableFingerprint] = []
    for name, table in entries:
        fingerprint = _fingerprint_optional_table_with_normalized_axes(table, name=name)
        if fingerprint is None:
            continue
        fingerprints.append(fingerprint)
    return tuple(fingerprints)


def _build_workflow_parameters(
    *,
    request: ResolvedSignalomeWorkflowRequest,
    config: ResolvedSignalomeExecutionConfig,
    clustering_result: ClusterSitesResult,
    network_correlation_diagnostics: SignalomeNetworkCorrelationDiagnostics,
    scale_guard_decision: SignalomeScaleGuardDecision,
    clustering_missing_value_diagnostics: SignalomeClusteringMissingValueDiagnostics,
    upstream_provenance: RunProvenance | None,
) -> dict[str, object]:
    payload = input_intensity_scale_evidence_payload(request.dataset)
    row_attrition = build_signalome_row_attrition_provenance(request)
    payload.update(
        {
            "site_token_validation": _build_site_token_validation_payload(request),
            "signalome_config": _build_signalome_config_payload(
                config=config,
                clustering_missing_value_diagnostics=(
                    clustering_missing_value_diagnostics
                ),
            ),
            "scale_guard": _build_scale_guard_payload(scale_guard_decision),
            "module_selection_diagnostics": asdict(
                clustering_result.module_selection_diagnostics
            ),
            "score_preconditioning_diagnostics": asdict(
                request.score_preconditioning_diagnostics
            ),
            "alignment_diagnostics": asdict(request.alignment_diagnostics),
            "network_correlation_diagnostics": asdict(network_correlation_diagnostics),
            "signalome_score_semantics": _build_signalome_score_semantics(
                request=request,
                config=config,
                clustering_result=clustering_result,
                network_correlation_diagnostics=network_correlation_diagnostics,
                scale_guard_decision=scale_guard_decision,
                clustering_missing_value_diagnostics=(
                    clustering_missing_value_diagnostics
                ),
            ),
            "upstream_kinase_provenance": (
                None
                if upstream_provenance is None
                else provenance_to_payload(upstream_provenance)
            ),
            **row_attrition.to_workflow_parameters(),
        }
    )
    return payload


def _build_site_token_validation_payload(
    request: ResolvedSignalomeWorkflowRequest,
) -> dict[str, object]:
    return {
        "mode": (
            "opaque_opt_in"
            if request.dataset.opaque_site_values_allowed
            else "strict_sty_residue_position"
        )
    }


def _build_signalome_config_payload(
    *,
    config: ResolvedSignalomeExecutionConfig,
    clustering_missing_value_diagnostics: SignalomeClusteringMissingValueDiagnostics,
) -> dict[str, object]:
    return {
        "scientific": {
            "substrate_support_cutoff": float(config.substrate_support_cutoff),
            "assignment_policy": str(config.assignment_policy),
        },
        "clustering": {
            "module_selection_primary_correlation_threshold": float(
                config.module_selection_primary_threshold
            ),
            "module_selection_fallback_correlation_threshold": float(
                config.module_selection_fallback_threshold
            ),
            "module_selection_max_clusters": int(config.module_selection_max_clusters),
            "candidate_scoring_policy": str(config.candidate_scoring_policy),
            "missing_value_policy": str(clustering_missing_value_diagnostics.policy),
            "clustering_engine": str(config.clustering_engine),
            "module_count": (
                None
                if config.requested_module_count is None
                else int(config.requested_module_count)
            ),
        },
        "validation": {
            "score_preconditioning_policy": str(config.score_preconditioning_policy),
            "allow_mixed_total_protein_quantitative_meaning": bool(
                config.allow_mixed_total_protein_quantitative_meaning
            ),
            "reference_context_compatibility_policy": str(
                config.reference_context_compatibility_policy
            ),
        },
        "output": {
            "network_correlation_threshold": float(
                config.network_correlation_threshold
            ),
            "network_policy": str(config.network_policy),
            "network_min_paired_finite_observations": int(
                config.network_min_paired_finite_observations
            ),
        },
        "performance": {
            "max_exact_tree_sites": int(config.max_exact_tree_sites),
            "max_full_candidate_scoring_sites": int(
                config.max_full_candidate_scoring_sites
            ),
        },
    }


def _build_scale_guard_payload(
    scale_guard_decision: SignalomeScaleGuardDecision,
) -> dict[str, object]:
    return {
        "site_count": int(scale_guard_decision.site_count),
        "input_protein_count": int(scale_guard_decision.input_protein_count),
        "input_kinase_count": int(scale_guard_decision.input_kinase_count),
        "selected_module_count": int(scale_guard_decision.selected_module_count),
        "candidate_module_counts_evaluated": int(
            scale_guard_decision.candidate_module_counts_evaluated
        ),
        "candidate_module_count_upper_bound": int(
            scale_guard_decision.candidate_module_count_upper_bound
        ),
        "clustering_engine": str(scale_guard_decision.clustering_engine),
        "clustering_engine_version": str(
            scale_guard_decision.clustering_engine_version
        ),
        "backend_diagnostics": (
            None
            if scale_guard_decision.backend_diagnostics is None
            else backend_diagnostics_to_payload(
                scale_guard_decision.backend_diagnostics
            )
        ),
        "tree_implementation": str(scale_guard_decision.tree_implementation),
        "tree_generation_backend": str(scale_guard_decision.tree_generation_backend),
        "tree_generation_mode": str(scale_guard_decision.tree_generation_mode),
        "tree_generation_is_approximate": bool(
            scale_guard_decision.tree_generation_is_approximate
        ),
        "tree_generation_scope": str(scale_guard_decision.tree_generation_scope),
        "tree_generation_guard_triggered": bool(
            scale_guard_decision.tree_generation_guard_triggered
        ),
        "candidate_scoring_policy": str(scale_guard_decision.candidate_scoring_policy),
        "candidate_scoring_requested_policy": str(
            scale_guard_decision.candidate_scoring_requested_policy
        ),
        "candidate_scoring_strategy": str(
            scale_guard_decision.candidate_scoring_strategy
        ),
        "candidate_scoring_is_approximate": bool(
            scale_guard_decision.candidate_scoring_is_approximate
        ),
        "candidate_scoring_guard_triggered": bool(
            scale_guard_decision.candidate_scoring_guard_triggered
        ),
        "candidate_scoring_sampled_site_total": (
            None
            if scale_guard_decision.candidate_scoring_sampled_site_total is None
            else int(scale_guard_decision.candidate_scoring_sampled_site_total)
        ),
        "candidate_scoring_sampled_pair_count": (
            None
            if scale_guard_decision.candidate_scoring_sampled_pair_count is None
            else int(scale_guard_decision.candidate_scoring_sampled_pair_count)
        ),
        "max_exact_tree_sites": int(scale_guard_decision.max_exact_tree_sites),
        "max_full_candidate_scoring_sites": int(
            scale_guard_decision.max_full_candidate_scoring_sites
        ),
        "exact_cluster_tree_built": bool(scale_guard_decision.exact_cluster_tree_built),
        "candidate_scoring_mode": str(scale_guard_decision.candidate_scoring_mode),
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
            else candidate_scoring_sampling_diagnostics_to_payload(
                scale_guard_decision.candidate_scoring_sampling
            )
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
    }


def _build_scientific_policy_records(
    *,
    request: ResolvedSignalomeWorkflowRequest,
    config: ResolvedSignalomeExecutionConfig,
    clustering_result: ClusterSitesResult,
    scale_guard_decision: SignalomeScaleGuardDecision,
    clustering_missing_value_diagnostics: SignalomeClusteringMissingValueDiagnostics,
) -> tuple[ScientificPolicyRecord, ...]:
    scientific_policies: list[ScientificPolicyRecord] = [
        build_signalome_module_candidate_score_policy(
            requested_policy=str(
                scale_guard_decision.candidate_scoring_requested_policy
            ),
            candidate_scoring_policy=str(config.candidate_scoring_policy),
            candidate_scoring_mode=str(clustering_result.candidate_scoring_mode),
            max_exact_tree_sites=(
                None
                if config.max_exact_tree_sites is None
                else int(config.max_exact_tree_sites)
            ),
            max_full_candidate_scoring_sites=int(
                config.max_full_candidate_scoring_sites
            ),
            candidate_scoring_evaluated=bool(
                clustering_result.candidate_scoring_evaluated
            ),
            candidate_scoring_skip_reason=(
                None
                if clustering_result.candidate_scoring_skip_reason is None
                else str(clustering_result.candidate_scoring_skip_reason)
            ),
            candidate_scoring_scope=str(
                scale_guard_decision.candidate_scoring_applies_to
            ),
            tree_generation_mode=str(scale_guard_decision.tree_generation_mode),
            tree_generation_is_approximate=bool(
                scale_guard_decision.tree_generation_is_approximate
            ),
        ),
        SignalomeMissingValueClusteringPolicy(
            missing_value_policy=str(clustering_missing_value_diagnostics.policy),
            applies_to=str(SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_APPLIES_TO),
            imputed_values_exposed_in_output_tables=bool(
                SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_IMPUTED_VALUES_EXPOSED_IN_OUTPUT_TABLES
            ),
        ).record,
        _build_signalome_assignment_policy_record(
            assignment_policy=str(config.assignment_policy)
        ),
        _build_signalome_network_policy_record(
            network_policy=str(config.network_policy),
            network_correlation_threshold=float(config.network_correlation_threshold),
            network_min_paired_finite_observations=int(
                config.network_min_paired_finite_observations
            ),
        ),
        ScorePreconditioningPolicy(
            policy=str(request.score_preconditioning_diagnostics.policy),
            input_row_count=int(
                request.score_preconditioning_diagnostics.input_row_count
            ),
            dropped_all_missing_row_count=int(
                request.score_preconditioning_diagnostics.dropped_all_missing_row_count
            ),
            retained_row_count=int(
                request.score_preconditioning_diagnostics.retained_row_count
            ),
        ).record,
        PROTEIN_MODULE_FROM_SITE_MEMBERSHIP_POLICY,
    ]
    if request.downstream_score_selection_policy is not None:
        scientific_policies.append(request.downstream_score_selection_policy.record)
    if config.candidate_scoring_policy_definition is not None:
        scientific_policies.append(config.candidate_scoring_policy_definition.record)
    return tuple(scientific_policies)


def _build_signalome_score_semantics(
    *,
    request: ResolvedSignalomeWorkflowRequest,
    config: ResolvedSignalomeExecutionConfig,
    clustering_result: ClusterSitesResult,
    network_correlation_diagnostics: SignalomeNetworkCorrelationDiagnostics,
    scale_guard_decision: SignalomeScaleGuardDecision,
    clustering_missing_value_diagnostics: SignalomeClusteringMissingValueDiagnostics,
) -> dict[str, object]:
    downstream_score_source = request.downstream_score_source
    module_selection_diagnostics = clustering_result.module_selection_diagnostics
    preconditioning = request.score_preconditioning_diagnostics
    return {
        "downstream_score_source": downstream_score_source.value,
        "downstream_score_meaning": _downstream_score_meaning(
            downstream_score_source=downstream_score_source
        ),
        "module_selection_score_meaning": (
            "module-count selection uses within-cluster correlation summaries over "
            "downstream score profiles across candidate module counts"
        ),
        "assignment_semantics": {
            "top_kinase": (
                "top supported kinase candidate by site-level prediction score; "
                "ties are reported with the configured selection policy"
            ),
            "module_top_kinase": (
                "top supported kinase candidate summarized across sites in the "
                "candidate module; this is a label, not a causal mechanism claim"
            ),
            "module_id": (
                "score-derived candidate kinase-supported module identifier for "
                "the current dataset and configuration"
            ),
        },
        "candidate_scoring_mode": str(scale_guard_decision.candidate_scoring_mode),
        "candidate_scoring_is_approximate": bool(
            scale_guard_decision.candidate_scoring_is_approximate
        ),
        "candidate_scoring_sampled_site_total": (
            None
            if scale_guard_decision.candidate_scoring_sampled_site_total is None
            else int(scale_guard_decision.candidate_scoring_sampled_site_total)
        ),
        "candidate_scoring_sampled_pair_count": (
            None
            if scale_guard_decision.candidate_scoring_sampled_pair_count is None
            else int(scale_guard_decision.candidate_scoring_sampled_pair_count)
        ),
        "candidate_scoring_scope": str(
            scale_guard_decision.candidate_scoring_applies_to
        ),
        "tree_generation_mode": str(scale_guard_decision.tree_generation_mode),
        "tree_generation_is_approximate": bool(
            scale_guard_decision.tree_generation_is_approximate
        ),
        "tree_generation_scope": str(scale_guard_decision.tree_generation_scope),
        "tree_generation_backend": str(scale_guard_decision.tree_generation_backend),
        "input_sizes": {
            "site_count": int(scale_guard_decision.site_count),
            "protein_count": int(scale_guard_decision.input_protein_count),
            "kinase_count": int(scale_guard_decision.input_kinase_count),
            "candidate_module_counts_evaluated": int(
                scale_guard_decision.candidate_module_counts_evaluated
            ),
            "candidate_module_count_upper_bound": int(
                scale_guard_decision.candidate_module_count_upper_bound
            ),
        },
        "scale_guard_status": {
            "exact_tree_guard_triggered": bool(
                scale_guard_decision.tree_generation_guard_triggered
            ),
            "candidate_scoring_guard_triggered": bool(
                scale_guard_decision.candidate_scoring_guard_triggered
            ),
            "passed": bool(scale_guard_decision.scale_guard_passed),
        },
        "network_correlation_meaning": (
            "network candidate and edge scores are score-profile associations: "
            "pairwise correlations between downstream kinase score profiles, "
            "computed on finite paired observations"
        ),
        "network_edge_semantics": {
            "nodes": (
                "kinases retained in the aligned prediction and downstream score "
                "matrices"
            ),
            "edges": (
                "correlation edge between kinase score profiles that passed the "
                "configured threshold and network policy"
            ),
            "direction": (
                "not inferred; source_kinase and target_kinase are deterministic "
                "table labels for an undirected score-profile association"
            ),
            "weight": (
                "correlation column derived from pairwise finite downstream score "
                "profile correlations; signed and absolute-threshold policies "
                "determine whether the stored value keeps sign or stores magnitude"
            ),
            "threshold_policy": {
                "network_policy": str(config.network_policy),
                "network_correlation_threshold": float(
                    config.network_correlation_threshold
                ),
                "network_min_paired_finite_observations": int(
                    config.network_min_paired_finite_observations
                ),
            },
            "edge_diagnostics": {
                "retained_edges": int(network_correlation_diagnostics.edges_created),
                "skipped_below_threshold": int(
                    network_correlation_diagnostics.edges_skipped_below_threshold
                ),
                "skipped_insufficient_paired_observations": int(
                    network_correlation_diagnostics.edges_skipped_insufficient_paired_observations
                ),
                "skipped_constant_profile": int(
                    network_correlation_diagnostics.edges_skipped_constant_profile
                ),
                "skipped_missing_score": int(
                    network_correlation_diagnostics.edges_skipped_missing_score
                ),
                "skipped_non_finite_score": int(
                    network_correlation_diagnostics.edges_skipped_non_finite_score
                ),
                "skipped_undefined_correlation": int(
                    network_correlation_diagnostics.edges_skipped_undefined_correlation
                ),
            },
            "interpretation_limit": (
                "correlations are not causal evidence and do not prove signalling "
                "relationships"
            ),
        },
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
            "clustering_distance_input": {
                "policy": str(clustering_missing_value_diagnostics.policy),
                "applies_to": SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_APPLIES_TO,
                "non_finite_handling": (
                    "non-finite values are treated as missing before imputation"
                ),
                "partial_missingness_handling": (
                    "missing entries are imputed with the median of the same column"
                ),
                "fully_missing_column_handling": (
                    "columns with all values missing are imputed with 0.0"
                ),
                "non_finite_input_value_count": int(
                    clustering_missing_value_diagnostics.non_finite_input_value_count
                ),
                "missing_after_non_finite_normalization_count": int(
                    clustering_missing_value_diagnostics.missing_after_non_finite_normalization_count
                ),
                "imputed_value_count": int(
                    clustering_missing_value_diagnostics.imputed_value_count
                ),
                "fully_missing_column_count": int(
                    clustering_missing_value_diagnostics.fully_missing_column_count
                ),
                "output_tables_include_imputed_values": bool(
                    SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_IMPUTED_VALUES_EXPOSED_IN_OUTPUT_TABLES
                ),
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
            "network_min_paired_finite_observations": int(
                config.network_min_paired_finite_observations
            ),
            "max_exact_tree_sites": int(config.max_exact_tree_sites),
            "max_full_candidate_scoring_sites": int(
                config.max_full_candidate_scoring_sites
            ),
        },
        "clustering_engine": str(scale_guard_decision.clustering_engine),
        "scientific_interpretation_limits": (
            "signalome module assignments, module scores, and kinase score-profile "
            "association edges are derived summary statistics for this dataset and "
            "configuration; they are not probabilities, calibrated confidence "
            "values, experimental validation of signalling relationships, or "
            "causal evidence"
        ),
    }


def _downstream_score_meaning(
    *, downstream_score_source: DownstreamScoreSource | str
) -> str:
    resolved_source = DownstreamScoreSource.parse(
        downstream_score_source,
        field_name="signalome provenance downstream_score_source",
    )
    if resolved_source is DownstreamScoreSource.RANK_WEIGHTED_FUSION_SCORES:
        return (
            "rank-weighted fusion of upstream downstream-score lanes; larger values "
            "indicate stronger relative downstream support within the run"
        )
    if resolved_source is DownstreamScoreSource.PROFILE_SCORES:
        return (
            "upstream downstream profile-score lane; larger values indicate stronger "
            "relative downstream support within the run"
        )
    if resolved_source is DownstreamScoreSource.KINASE_LIBRARY_MOTIF_SCORES:
        return (
            "upstream Kinase Library motif-score lane normalized to within-run "
            "unit support; larger values indicate stronger relative motif support"
        )
    if resolved_source is DownstreamScoreSource.COMBINED_PROFILE_MOTIF_SCORES:
        return (
            "upstream combined profile and Kinase Library motif lane; larger "
            "values indicate stronger relative downstream support within the run"
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


def _build_signalome_assignment_policy_record(
    *,
    assignment_policy: str,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.SIGNALOME_ASSIGNMENT_POLICY,
        name=f"signalome_assignment_policy_{assignment_policy}_v1",
        version="1",
        description=(
            "Assignment policy used to summarize site support as module-level "
            "top supported kinase candidate labels."
        ),
        parameters={"assignment_policy": assignment_policy},
        assumptions=(
            "Assignment policy affects top supported kinase candidate labels and tie handling.",
            "Labels summarize score support and do not infer causal regulation.",
        ),
        output_scale="Module assignment labels and top supported kinase summaries.",
        quantitative_meaning="signalome_assignment_rule",
    )


def _build_signalome_network_policy_record(
    *,
    network_policy: str,
    network_correlation_threshold: float,
    network_min_paired_finite_observations: int,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.SIGNALOME_NETWORK_POLICY,
        name=f"signalome_network_policy_{network_policy}_v1",
        version="1",
        description=(
            "Policy controlling exploratory kinase score-profile association edge "
            "eligibility."
        ),
        parameters={
            "network_policy": network_policy,
            "network_correlation_threshold": float(network_correlation_threshold),
            "network_min_paired_finite_observations": int(
                network_min_paired_finite_observations
            ),
        },
        assumptions=(
            "Network policy defines how sign and magnitude thresholds are applied.",
            "Edges represent score-profile associations, not inferred direction or causality.",
        ),
        output_scale="Signalome kinase score-profile association edge table.",
        quantitative_meaning="network_edge_eligibility_rule",
    )


__all__ = ["SignalomeProvenanceBuilder"]
