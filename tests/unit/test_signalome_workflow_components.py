from __future__ import annotations

import pandas as pd
import pandas.testing as pdt

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    KinaseWorkflowResult,
    Organism,
    ReferenceBundle,
    SignalomeWorkflowRequest,
)
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    SignalomeWorkflowResult,
)
from phospy.provenance.scientific_policy_models import ScientificPolicyId
from phospy.science.signalomes.clustering import (
    ClusterSitesResult,
)
from phospy.science.signalomes.clustering.models import SignalomeClusteringEngineResult
from phospy.science.signalomes.constants import SITE_CLUSTER_COLUMN
from phospy.science.signalomes.context import (
    build_protein_site_context_table,
    build_site_membership_table,
)
from phospy.science.signalomes.models import (
    SignalomeModuleSelectionDiagnostics,
)
from phospy.science.signalomes.science import (
    build_expanded_signalome_table,
    build_kinase_network_with_diagnostics,
    build_module_assignments,
    build_signalome_module_table,
    select_kinase_substrates,
)
from phospy.workflows.signalome.clustering_runner import SignalomeClusteringRunner
from phospy.workflows.signalome.context_tables import SignalomeContextTableBuilder
from phospy.workflows.signalome.interpreter import SignalomeWorkflowInterpreter
from phospy.workflows.signalome.module_tables import SignalomeModuleTableBuilder
from phospy.workflows.signalome.network_builder import SignalomeNetworkBuilder
from phospy.workflows.signalome.provenance import SignalomeProvenanceBuilder
from phospy.workflows.signalome.result_assembly import SignalomeResultAssembler
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.signalome_config import build_signalome_config


def _dataset() -> AnalysisReadyPhosphoDataset:
    site_ids = ["P1;S1;", "P1;S2;", "P2;S3;", "P3;S4;"]
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0, 4.0],
            "sample_b": [1.2, 2.2, 3.2, 4.2],
        },
        index=site_ids,
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["P1", "P1", "P2", "P3"],
            "site": ["S1", "S2", "S3", "S4"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["S1", "S2", "S3", "S4"]
            ],
            "protein_id": ["P1", "P1", "P2", "P3"],
        },
        index=site_ids,
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _bundle(site_ids: list[str]) -> ReferenceBundle:
    unique_sites = pd.Index([str(site_id) for site_id in site_ids], name="site_id")
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1", "K2"],
                "substrate_site": [unique_sites[0], unique_sites[1]],
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    ("A" * 15)
                    + str(site_id).split(";")[1].strip().upper()[0]
                    + ("A" * 15)
                    for site_id in unique_sites
                ]
            },
            index=unique_sites,
        ),
    )


def _matrix(
    *,
    values: list[list[float]],
    site_ids: list[str],
    kinases: list[str],
) -> pd.DataFrame:
    return pd.DataFrame(
        values,
        index=pd.Index(site_ids, name="site_id"),
        columns=pd.Index(kinases, name="kinase"),
        dtype=float,
    )


def _kinase_result(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    prediction_matrix: pd.DataFrame,
    score_matrix: pd.DataFrame,
) -> KinaseWorkflowResult:
    return KinaseWorkflowResult(
        dataset=dataset,
        references=_bundle(site_ids=dataset.phospho.index.astype(str).tolist()),
        scoring_result=KinaseScoringResult(
            profile_scores=score_matrix,
            rank_weighted_fusion_scores=score_matrix,
        ),
        prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
        activity_result=None,
    )


def _resolved_request():
    dataset = _dataset()
    site_ids = dataset.phospho.index.astype(str).tolist()
    kinases = ["K1", "K2"]
    prediction_matrix = _matrix(
        values=[
            [0.9, 0.1],
            [0.1, 0.9],
            [0.8, 0.2],
            [0.7, 0.3],
        ],
        site_ids=site_ids,
        kinases=kinases,
    )
    score_matrix = _matrix(
        values=[
            [1.0, 0.0],
            [0.0, 1.0],
            [0.9, 0.1],
            [float("nan"), float("nan")],
        ],
        site_ids=site_ids,
        kinases=kinases,
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=build_signalome_config(
            substrate_support_cutoff=0.5,
            network_correlation_threshold=0.3,
            module_count=2,
            score_preconditioning_policy="allow_and_report",
        ),
    )
    return SignalomeWorkflowInterpreter().run(request)


def test_signalome_clustering_runner_returns_expected_diagnostics() -> None:
    resolved = _resolved_request()
    metadata = SignalomeClusteringRunner.collect_execution_metadata(resolved)
    observed_backend_kwargs: dict[str, object] = {}
    expected_cluster_result = ClusterSitesResult(
        site_clusters=pd.Series(
            [1, 2, 2],
            index=resolved.downstream_score_matrix.index.copy(),
            dtype="int64",
            name=SITE_CLUSTER_COLUMN,
        ),
        module_selection_diagnostics=SignalomeModuleSelectionDiagnostics(
            strategy="correlation_thresholds",
            selected_module_count=2,
            requested_module_count=2,
            threshold_used=0.5,
            max_clusters_evaluated=2,
            candidate_scores={},
            reason="component test",
        ),
        tree_engine="exact",
        candidate_scoring_mode="full",
        exact_cluster_tree_built=True,
    )
    expected_protein_modules = pd.Series(
        {"P1": 1, "P2": 2, "P3": 2},
        dtype="int64",
        name="module_id",
    )

    def _run_backend(**kwargs: object) -> SignalomeClusteringEngineResult:
        observed_backend_kwargs.update(kwargs)
        return SignalomeClusteringEngineResult(
            site_clusters=expected_cluster_result.site_clusters,
            protein_modules=expected_protein_modules,
            selected_module_count=2,
            module_selection_diagnostics=expected_cluster_result.module_selection_diagnostics,
            backend_name="exact_python",
            backend_version="1",
            approximation_used=False,
            exact_cluster_tree_built=True,
            tree_implementation="exact",
            candidate_scoring_mode="full",
            candidate_scoring_evaluated=True,
            candidate_scoring_skip_reason=None,
            candidate_scoring_sampling=None,
            backend_diagnostics=None,
            threshold_metadata={
                "primary_threshold": float(
                    resolved.execution_config.module_selection_primary_threshold
                ),
                "fallback_threshold": float(
                    resolved.execution_config.module_selection_fallback_threshold
                ),
            },
            limit_metadata={
                "max_exact_tree_sites": int(
                    resolved.execution_config.max_exact_tree_sites
                ),
                "max_full_candidate_scoring_sites": int(
                    resolved.execution_config.max_full_candidate_scoring_sites
                ),
                "max_clusters": int(
                    resolved.execution_config.module_selection_max_clusters
                ),
            },
        )

    runner = SignalomeClusteringRunner(run_backend_clustering=_run_backend)
    output = runner.run(
        request=resolved,
        config=resolved.execution_config,
        execution_metadata=metadata,
    )

    pdt.assert_series_equal(
        output.clustering_result.site_clusters,
        expected_cluster_result.site_clusters,
    )
    pdt.assert_series_equal(output.protein_modules, expected_protein_modules)
    assert output.clustering_result.module_selection_diagnostics == (
        expected_cluster_result.module_selection_diagnostics
    )
    assert observed_backend_kwargs["requested_module_count"] == 2
    assert "tree_engine" not in observed_backend_kwargs
    assert observed_backend_kwargs["candidate_scoring_policy"] == "full"
    assert observed_backend_kwargs["site_to_protein"] is resolved.site_to_protein


def test_signalome_module_table_builder_preserves_module_summary_shape() -> None:
    resolved = _resolved_request()
    metadata = SignalomeClusteringRunner.collect_execution_metadata(resolved)
    clustering = SignalomeClusteringRunner().run(
        request=resolved,
        config=resolved.execution_config,
        execution_metadata=metadata,
    )
    observed = SignalomeModuleTableBuilder().run(
        request=resolved,
        config=resolved.execution_config,
        clustering_result=clustering.clustering_result,
        protein_modules=clustering.protein_modules,
        execution_metadata=metadata,
    )

    expected_assignments = build_module_assignments(
        prediction_matrix=resolved.prediction_matrix,
        site_to_protein=resolved.site_to_protein,
        protein_modules=clustering.protein_modules,
    )
    expected_substrates = select_kinase_substrates(
        prediction_matrix=resolved.prediction_matrix,
        cutoff=resolved.execution_config.substrate_support_cutoff,
    )
    expected_modules = build_signalome_module_table(
        module_assignments=expected_assignments,
        kinase_substrates=expected_substrates,
        kinase_order=resolved.prediction_matrix.columns.astype(str).tolist(),
        assignment_policy=resolved.execution_config.assignment_policy,
    )

    pdt.assert_frame_equal(observed.module_assignments, expected_assignments)
    pdt.assert_frame_equal(observed.signalome_modules, expected_modules)
    assert observed.signalome_modules.shape == expected_modules.shape


def test_signalome_network_builder_preserves_edge_schema() -> None:
    resolved = _resolved_request()
    metadata = SignalomeClusteringRunner.collect_execution_metadata(resolved)
    clustering = SignalomeClusteringRunner().run(
        request=resolved,
        config=resolved.execution_config,
        execution_metadata=metadata,
    )
    module_stage = SignalomeModuleTableBuilder().run(
        request=resolved,
        config=resolved.execution_config,
        clustering_result=clustering.clustering_result,
        protein_modules=clustering.protein_modules,
        execution_metadata=metadata,
    )
    observed = SignalomeNetworkBuilder().run(
        request=resolved,
        config=resolved.execution_config,
        clustering_result=clustering.clustering_result,
        support_summary=module_stage.support_summary,
        execution_metadata=metadata,
    )
    expected_edges, expected_nodes, expected_candidates, expected_diagnostics = (
        build_kinase_network_with_diagnostics(
            downstream_score_matrix=resolved.downstream_score_matrix,
            kinase_order=resolved.prediction_matrix.columns.astype(str).tolist(),
            kinase_substrates=module_stage.support_summary.kinase_substrates,
            threshold=resolved.execution_config.network_correlation_threshold,
            network_policy=resolved.execution_config.network_policy,
        )
    )

    pdt.assert_frame_equal(observed.edges, expected_edges)
    pdt.assert_frame_equal(observed.nodes, expected_nodes)
    pdt.assert_frame_equal(observed.candidate_correlations, expected_candidates)
    assert observed.correlation_diagnostics == expected_diagnostics
    assert observed.edges.shape == expected_edges.shape


def test_signalome_context_table_builder_flags_multisite_proteins() -> None:
    resolved = _resolved_request()
    metadata = SignalomeClusteringRunner.collect_execution_metadata(resolved)
    clustering = SignalomeClusteringRunner().run(
        request=resolved,
        config=resolved.execution_config,
        execution_metadata=metadata,
    )
    module_stage = SignalomeModuleTableBuilder().run(
        request=resolved,
        config=resolved.execution_config,
        clustering_result=clustering.clustering_result,
        protein_modules=clustering.protein_modules,
        execution_metadata=metadata,
    )
    observed = SignalomeContextTableBuilder().run(
        request=resolved,
        config=resolved.execution_config,
        clustering_result=clustering.clustering_result,
        module_assignments=module_stage.module_assignments,
        support_summary=module_stage.support_summary,
        execution_metadata=metadata,
    )
    expected_site_membership = build_site_membership_table(
        module_assignments=module_stage.module_assignments,
        site_clusters=clustering.clustering_result.site_clusters,
        site_metadata=resolved.dataset.site_metadata,
        prediction_matrix=resolved.prediction_matrix,
        kinase_substrates=module_stage.support_summary.kinase_substrates,
        substrate_support_cutoff=resolved.execution_config.substrate_support_cutoff,
        assignment_policy=resolved.execution_config.assignment_policy,
    )
    expected_protein_context = build_protein_site_context_table(
        site_membership=expected_site_membership
    )

    pdt.assert_frame_equal(observed.site_membership, expected_site_membership)
    pdt.assert_frame_equal(observed.protein_site_context, expected_protein_context)
    context = observed.protein_site_context.set_index("protein_id")
    assert bool(context.loc["P1", "multi_site_protein"])
    assert bool(context.loc["P1", "ambiguous_module_context"])


def test_signalome_provenance_builder_records_scale_and_backend_fields() -> None:
    resolved = _resolved_request()
    metadata = SignalomeClusteringRunner.collect_execution_metadata(resolved)
    clustering = SignalomeClusteringRunner().run(
        request=resolved,
        config=resolved.execution_config,
        execution_metadata=metadata,
    )
    module_stage = SignalomeModuleTableBuilder().run(
        request=resolved,
        config=resolved.execution_config,
        clustering_result=clustering.clustering_result,
        protein_modules=clustering.protein_modules,
        execution_metadata=metadata,
    )
    network_stage = SignalomeNetworkBuilder().run(
        request=resolved,
        config=resolved.execution_config,
        clustering_result=clustering.clustering_result,
        support_summary=module_stage.support_summary,
        execution_metadata=metadata,
    )
    context_stage = SignalomeContextTableBuilder().run(
        request=resolved,
        config=resolved.execution_config,
        clustering_result=clustering.clustering_result,
        module_assignments=module_stage.module_assignments,
        support_summary=module_stage.support_summary,
        execution_metadata=metadata,
    )
    expanded_signalome = build_expanded_signalome_table(
        module_assignments=module_stage.module_assignments,
        signalome_modules=module_stage.signalome_modules,
        kinase_network_edges=network_stage.edges,
        kinase_substrates=module_stage.support_summary.kinase_substrates,
        assignment_policy=resolved.execution_config.assignment_policy,
    )
    scale_guard = SignalomeClusteringRunner.summarize_scale_guard(
        config=resolved.execution_config,
        site_count=metadata.downstream_score_sites,
        site_to_protein=resolved.site_to_protein,
        downstream_score_kinases=metadata.downstream_score_kinases,
        clustering_result=clustering.clustering_result,
    )
    provenance = SignalomeProvenanceBuilder().build(
        request=resolved,
        config=resolved.execution_config,
        clustering_result=clustering.clustering_result,
        module_assignments=module_stage.module_assignments,
        signalome_modules=module_stage.signalome_modules,
        network_edges=network_stage.edges,
        network_nodes=network_stage.nodes,
        candidate_correlations=network_stage.candidate_correlations,
        network_correlation_diagnostics=network_stage.correlation_diagnostics,
        expanded_signalome=expanded_signalome,
        site_membership=context_stage.site_membership,
        protein_site_context=context_stage.protein_site_context,
        scale_guard_decision=scale_guard,
    )

    assert provenance.workflow_name == "signalome_workflow"
    assert "signalome_config" in provenance.workflow_parameters
    assert "scale_guard" in provenance.workflow_parameters
    assert "module_selection_diagnostics" in provenance.workflow_parameters
    assert "alignment_diagnostics" in provenance.workflow_parameters
    assert "network_correlation_diagnostics" in provenance.workflow_parameters
    signalome_config = provenance.workflow_parameters["signalome_config"]
    assert "scientific" in signalome_config
    assert "clustering" in signalome_config
    assert "validation" in signalome_config
    assert "output" in signalome_config
    assert "performance" in signalome_config
    assert "tree_engine" not in signalome_config["clustering"]
    assert signalome_config["clustering"]["candidate_scoring_policy"] == "full"
    assert signalome_config["clustering"]["missing_value_policy"] == (
        "column_median_imputation_with_zero_for_all_missing_columns"
    )
    assert signalome_config["performance"]["max_exact_tree_sites"] == 2000
    assert signalome_config["clustering"]["module_count"] == 2
    scale_guard = provenance.workflow_parameters["scale_guard"]
    assert "tree_implementation" in scale_guard
    assert "candidate_scoring_policy" in scale_guard
    assert "candidate_scoring_requested_policy" in scale_guard
    assert "candidate_scoring_strategy" in scale_guard
    assert "candidate_scoring_is_approximate" in scale_guard
    assert "candidate_scoring_guard_triggered" in scale_guard
    assert "candidate_scoring_sampled_site_total" in scale_guard
    assert "candidate_scoring_sampled_pair_count" in scale_guard
    assert "exact_cluster_tree_built" in scale_guard
    assert "candidate_scoring_mode" in scale_guard
    assert "tree_generation_backend" in scale_guard
    assert "tree_generation_mode" in scale_guard
    assert "tree_generation_is_approximate" in scale_guard
    assert "tree_generation_scope" in scale_guard
    assert "tree_generation_guard_triggered" in scale_guard
    assert "input_protein_count" in scale_guard
    assert "input_kinase_count" in scale_guard
    assert "candidate_module_counts_evaluated" in scale_guard
    assert "candidate_module_count_upper_bound" in scale_guard
    output_names = {fingerprint.name for fingerprint in provenance.output_tables}
    assert "outputs.signalome.module_assignments" in output_names
    assert "outputs.signalome.signalome_modules" in output_names
    assert "outputs.signalome.site_membership" in output_names
    policy_ids = {policy.id for policy in provenance.scientific_policies}
    assert ScientificPolicyId.SIGNALOME_MODULE_CANDIDATE_SCORE in policy_ids
    assert ScientificPolicyId.SIGNALOME_MISSING_VALUE_CLUSTERING in policy_ids
    assert ScientificPolicyId.SIGNALOME_SCORE_PRECONDITIONING in policy_ids
    assert ScientificPolicyId.PROTEIN_MODULE_FROM_SITE_MEMBERSHIP in policy_ids


def test_signalome_result_assembly_preserves_public_result_shape() -> None:
    resolved = _resolved_request()
    metadata = SignalomeClusteringRunner.collect_execution_metadata(resolved)
    clustering = SignalomeClusteringRunner().run(
        request=resolved,
        config=resolved.execution_config,
        execution_metadata=metadata,
    )
    module_stage = SignalomeModuleTableBuilder().run(
        request=resolved,
        config=resolved.execution_config,
        clustering_result=clustering.clustering_result,
        protein_modules=clustering.protein_modules,
        execution_metadata=metadata,
    )
    network_stage = SignalomeNetworkBuilder().run(
        request=resolved,
        config=resolved.execution_config,
        clustering_result=clustering.clustering_result,
        support_summary=module_stage.support_summary,
        execution_metadata=metadata,
    )
    context_stage = SignalomeContextTableBuilder().run(
        request=resolved,
        config=resolved.execution_config,
        clustering_result=clustering.clustering_result,
        module_assignments=module_stage.module_assignments,
        support_summary=module_stage.support_summary,
        execution_metadata=metadata,
    )
    assembler = SignalomeResultAssembler()
    expanded_signalome = assembler.build_expanded_signalome(
        request=resolved,
        config=resolved.execution_config,
        module_assignments=module_stage.module_assignments,
        signalome_modules=module_stage.signalome_modules,
        network_edges=network_stage.edges,
        support_summary=module_stage.support_summary,
        module_count=module_stage.module_count,
        execution_metadata=metadata,
    )
    scale_guard = SignalomeClusteringRunner.summarize_scale_guard(
        config=resolved.execution_config,
        site_count=metadata.downstream_score_sites,
        site_to_protein=resolved.site_to_protein,
        downstream_score_kinases=metadata.downstream_score_kinases,
        clustering_result=clustering.clustering_result,
    )
    provenance = SignalomeProvenanceBuilder().build(
        request=resolved,
        config=resolved.execution_config,
        clustering_result=clustering.clustering_result,
        module_assignments=module_stage.module_assignments,
        signalome_modules=module_stage.signalome_modules,
        network_edges=network_stage.edges,
        network_nodes=network_stage.nodes,
        candidate_correlations=network_stage.candidate_correlations,
        network_correlation_diagnostics=network_stage.correlation_diagnostics,
        expanded_signalome=expanded_signalome,
        site_membership=context_stage.site_membership,
        protein_site_context=context_stage.protein_site_context,
        scale_guard_decision=scale_guard,
    )
    result = assembler.assemble_result(
        request=resolved,
        clustering_result=clustering.clustering_result,
        module_assignments=module_stage.module_assignments,
        signalome_modules=module_stage.signalome_modules,
        network_edges=network_stage.edges,
        network_nodes=network_stage.nodes,
        candidate_correlations=network_stage.candidate_correlations,
        network_correlation_diagnostics=network_stage.correlation_diagnostics,
        expanded_signalome=expanded_signalome,
        site_membership=context_stage.site_membership,
        protein_site_context=context_stage.protein_site_context,
        provenance=provenance,
    )

    assert isinstance(result, SignalomeWorkflowResult)
    assert not result.module_assignments.table.empty
    assert not result.signalome_modules.table.empty
    assert not result.kinase_network.edges.empty
    assert result.kinase_network.nodes is not None
    assert result.kinase_network.candidate_correlations is not None
    assert result.site_membership is not None
    assert result.protein_site_context is not None
    assert result.provenance.workflow_name == "signalome_workflow"
