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
from phospy.workflows.signalome.component_models import SignalomeScaleGuardDecision
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeExecutionConfig,
    ResolvedSignalomeWorkflowRequest,
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
        clustering_missing_value_diagnostics = (
            summarize_clustering_missing_value_diagnostics(
                request.downstream_score_matrix.to_numpy(dtype=float, copy=False)
            )
        )
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
                network_correlation_threshold=float(
                    config.network_correlation_threshold
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
            scientific_policies.append(
                config.candidate_scoring_policy_definition.record
            )
        return RunProvenance(
            environment=self._collect_environment(),
            input_tables=input_tables,
            preprocessing_stages=self._dataset_preprocessing_stages(request),
            reference=request.kinase_result.references.provenance,
            workflow_name="signalome_workflow",
            workflow_parameters={
                "signalome_config": {
                    "scientific": {
                        "substrate_support_cutoff": float(
                            config.substrate_support_cutoff
                        ),
                        "assignment_policy": str(config.assignment_policy),
                    },
                    "clustering": {
                        "module_selection_primary_correlation_threshold": float(
                            config.module_selection_primary_threshold
                        ),
                        "module_selection_fallback_correlation_threshold": float(
                            config.module_selection_fallback_threshold
                        ),
                        "module_selection_max_clusters": int(
                            config.module_selection_max_clusters
                        ),
                        "candidate_scoring_policy": str(
                            config.candidate_scoring_policy
                        ),
                        "missing_value_policy": str(
                            clustering_missing_value_diagnostics.policy
                        ),
                        "clustering_engine": str(config.clustering_engine),
                        "module_count": (
                            None
                            if config.requested_module_count is None
                            else int(config.requested_module_count)
                        ),
                    },
                    "validation": {
                        "score_preconditioning_policy": str(
                            config.score_preconditioning_policy
                        ),
                    },
                    "output": {
                        "network_correlation_threshold": float(
                            config.network_correlation_threshold
                        ),
                        "network_policy": str(config.network_policy),
                    },
                    "performance": {
                        "max_exact_tree_sites": int(config.max_exact_tree_sites),
                        "max_full_candidate_scoring_sites": int(
                            config.max_full_candidate_scoring_sites
                        ),
                    },
                },
                "scale_guard": {
                    "site_count": int(scale_guard_decision.site_count),
                    "input_protein_count": int(
                        scale_guard_decision.input_protein_count
                    ),
                    "input_kinase_count": int(scale_guard_decision.input_kinase_count),
                    "selected_module_count": int(
                        scale_guard_decision.selected_module_count
                    ),
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
                    "tree_implementation": str(
                        scale_guard_decision.tree_implementation
                    ),
                    "tree_generation_backend": str(
                        scale_guard_decision.tree_generation_backend
                    ),
                    "tree_generation_mode": str(
                        scale_guard_decision.tree_generation_mode
                    ),
                    "tree_generation_is_approximate": bool(
                        scale_guard_decision.tree_generation_is_approximate
                    ),
                    "tree_generation_scope": str(
                        scale_guard_decision.tree_generation_scope
                    ),
                    "tree_generation_guard_triggered": bool(
                        scale_guard_decision.tree_generation_guard_triggered
                    ),
                    "candidate_scoring_policy": str(
                        scale_guard_decision.candidate_scoring_policy
                    ),
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
                        if scale_guard_decision.candidate_scoring_sampled_site_total
                        is None
                        else int(
                            scale_guard_decision.candidate_scoring_sampled_site_total
                        )
                    ),
                    "candidate_scoring_sampled_pair_count": (
                        None
                        if scale_guard_decision.candidate_scoring_sampled_pair_count
                        is None
                        else int(
                            scale_guard_decision.candidate_scoring_sampled_pair_count
                        )
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
                },
                "module_selection_diagnostics": asdict(
                    clustering_result.module_selection_diagnostics
                ),
                "score_preconditioning_diagnostics": asdict(
                    request.score_preconditioning_diagnostics
                ),
                "alignment_diagnostics": asdict(request.alignment_diagnostics),
                "network_correlation_diagnostics": asdict(
                    network_correlation_diagnostics
                ),
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
            },
            random_state=None,
            random_seed_policy=None,
            output_tables=output_tables,
            scientific_policies=tuple(scientific_policies),
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
            canonical_table = _canonicalise_for_provenance_fingerprint(table)
            fingerprint = fingerprint_optional_table(canonical_table, name=name)
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


def _canonicalise_for_provenance_fingerprint(
    table: pd.DataFrame | None,
) -> pd.DataFrame | None:
    if table is None:
        return None
    canonical = table
    try:
        canonical = canonical.sort_index(axis=0, kind="mergesort")
    except Exception:
        pass
    try:
        canonical = canonical.sort_index(axis=1, kind="mergesort")
    except Exception:
        pass
    return canonical


def _build_signalome_assignment_policy_record(
    *,
    assignment_policy: str,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.SIGNALOME_ASSIGNMENT_POLICY,
        name=f"signalome_assignment_policy_{assignment_policy}_v1",
        version="1",
        description="Assignment policy used to map site support to module-level labels.",
        parameters={"assignment_policy": assignment_policy},
        assumptions=("Assignment policy affects top-kinase labels and tie handling.",),
        output_scale="Module assignment labels and top-kinase summaries.",
        quantitative_meaning="signalome_assignment_rule",
    )


def _build_signalome_network_policy_record(
    *,
    network_policy: str,
    network_correlation_threshold: float,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.SIGNALOME_NETWORK_POLICY,
        name=f"signalome_network_policy_{network_policy}_v1",
        version="1",
        description="Policy controlling signalome kinase-network edge eligibility.",
        parameters={
            "network_policy": network_policy,
            "network_correlation_threshold": float(network_correlation_threshold),
        },
        assumptions=(
            "Network policy defines how sign and magnitude thresholds are applied.",
        ),
        output_scale="Signalome kinase-network edge table.",
        quantitative_meaning="network_edge_eligibility_rule",
    )


__all__ = ["SignalomeProvenanceBuilder"]
