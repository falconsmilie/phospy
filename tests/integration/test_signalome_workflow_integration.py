from __future__ import annotations

from collections.abc import Mapping

import pandas.testing as pdt
import pytest
from pandas.api.types import (
    is_bool_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_object_dtype,
    is_string_dtype,
)

import phospy.workflows.signalome.interpreter as signalome_interpreter
from phospy import (
    AnalysisReadyPhosphoDataset,
    KinaseWorkflow,
    SignalomeWorkflow,
)
from phospy.api import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    KinaseWorkflowResult,
    ReferencePreset,
    SignalomeWorkflowRequest,
)
from phospy.api.configs import SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL
from phospy.api.results import KinaseScoringResult
from phospy.errors import WorkflowBoundaryError, WorkflowValidationError
from phospy.signalomes.constants import (
    EXPANDED_SIGNALOME_ROW_KIND_COLUMN,
    EXPANDED_SIGNALOME_ROW_KIND_SITE,
)
from phospy.workflows.signalome.constants import (
    SIGNALOME_INTERPRETER_SCORE_PRECONDITIONING_SEAM,
)
from tests.support.rewrite_fixture_data import (
    build_rat_l6_dataset,
    load_signalome_l6_provenance_golden,
)
from tests.support.signalome_config import build_signalome_config

pytestmark = pytest.mark.integration


def _is_text_dtype(values: object) -> bool:
    return is_object_dtype(values) or is_string_dtype(values)


def _fingerprints_by_name(
    fingerprints: tuple[object, ...],
) -> dict[str, Mapping[str, object]]:
    return {
        str(item.name): {
            "rows": int(item.rows),
            "columns": int(item.columns),
            "hash_algorithm": str(item.hash_algorithm),
            "hash_value": str(item.hash_value),
        }
        for item in fingerprints
    }


def _assert_expected_fingerprint_map(
    *,
    observed: Mapping[str, Mapping[str, object]],
    expected: Mapping[str, object],
) -> None:
    expected_map = {
        str(name): values
        for name, values in expected.items()
        if isinstance(values, Mapping)
    }
    assert set(observed) == set(expected_map)
    for table_name, table_expected in expected_map.items():
        table_observed = observed[table_name]
        assert table_observed == {
            "rows": int(table_expected["rows"]),
            "columns": int(table_expected["columns"]),
            "hash_algorithm": str(table_expected["hash_algorithm"]),
            "hash_value": str(table_expected["hash_value"]),
        }


def test_signalome_workflow_runs_dataset_to_kinase_to_signalome_path() -> None:
    dataset = build_rat_l6_dataset(n_sites=260)
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                deterministic_max_selected_kinases=12,
                adaptive_ensemble_runs=12,
            ),
            activity_config=None,
        )
    )
    result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=build_signalome_config(substrate_support_cutoff=0.5),
        )
    )

    assignments = result.module_assignments.table
    assert not assignments.empty
    assert assignments.index.name == "site_id"
    assert {
        "protein_id",
        "module_id",
        "top_kinase",
        "top_score",
        "top_kinase_candidates",
        "top_kinase_weights",
        "top_kinase_tie_count",
        "top_kinase_is_ambiguous",
        "top_kinase_selection_policy",
        "module_top_kinase",
        "module_top_kinase_candidates",
        "module_top_kinase_tie_count",
        "module_top_kinase_is_ambiguous",
        "module_top_kinase_selection_policy",
    }.issubset(set(assignments.columns))
    assert _is_text_dtype(assignments.loc[:, "protein_id"])
    assert is_integer_dtype(assignments.loc[:, "module_id"])
    assert _is_text_dtype(assignments.loc[:, "top_kinase"])
    assert is_float_dtype(assignments.loc[:, "top_score"])
    assert is_integer_dtype(assignments.loc[:, "top_kinase_tie_count"])
    assert is_bool_dtype(assignments.loc[:, "top_kinase_is_ambiguous"])
    assert _is_text_dtype(assignments.loc[:, "top_kinase_selection_policy"])
    assert _is_text_dtype(assignments.loc[:, "module_top_kinase"])
    assert is_integer_dtype(assignments.loc[:, "module_top_kinase_tie_count"])
    assert is_bool_dtype(assignments.loc[:, "module_top_kinase_is_ambiguous"])
    assert _is_text_dtype(assignments.loc[:, "module_top_kinase_selection_policy"])

    modules = result.signalome_modules.table
    assert not modules.empty
    assert modules.index.name == "module_id"
    assert modules.columns.name == "kinase"
    assert is_float_dtype(modules.to_numpy(dtype=float))

    network_nodes = result.kinase_network.nodes
    assert network_nodes is not None
    assert not network_nodes.empty
    assert network_nodes.index.name == "kinase"
    assert {"degree", "n_substrates"} == set(network_nodes.columns)
    assert is_integer_dtype(network_nodes.loc[:, "degree"])
    assert is_integer_dtype(network_nodes.loc[:, "n_substrates"])

    network_edges = result.kinase_network.edges
    assert not network_edges.empty
    assert {"source_kinase", "target_kinase", "correlation"} == set(
        network_edges.columns
    )
    assert _is_text_dtype(network_edges.loc[:, "source_kinase"])
    assert _is_text_dtype(network_edges.loc[:, "target_kinase"])
    assert is_float_dtype(network_edges.loc[:, "correlation"])
    candidate_correlations = result.kinase_network.candidate_correlations
    assert candidate_correlations is not None
    assert not candidate_correlations.empty
    assert {
        "source_kinase",
        "target_kinase",
        "correlation",
        "correlation_status",
        "valid_observations",
        "correlation_reason",
    } == set(candidate_correlations.columns)
    assert _is_text_dtype(candidate_correlations.loc[:, "source_kinase"])
    assert _is_text_dtype(candidate_correlations.loc[:, "target_kinase"])
    assert is_float_dtype(candidate_correlations.loc[:, "correlation"])
    assert _is_text_dtype(candidate_correlations.loc[:, "correlation_status"])
    assert is_integer_dtype(candidate_correlations.loc[:, "valid_observations"])
    assert (
        result.kinase_network.correlation_diagnostics.total_candidate_correlations
        == int(candidate_correlations.shape[0])
    )
    assert result.kinase_network.correlation_diagnostics.finite_correlations >= int(
        network_edges.shape[0]
    )

    expanded = result.expanded_signalome
    assert expanded is not None
    assert not expanded.empty
    assert {
        "kinase",
        "row_kind",
        "assignment_policy",
        "linked_kinases",
        "regulated_module_ids",
        "site_id",
        "site_order",
        "protein_id",
        "module_id",
        "support_kinases",
        "support_weight",
        "top_kinase",
        "top_score",
    } == set(expanded.columns)
    assert _is_text_dtype(expanded.loc[:, "kinase"])
    assert _is_text_dtype(expanded.loc[:, "row_kind"])
    assert is_integer_dtype(expanded.loc[:, "site_order"])
    assert is_integer_dtype(expanded.loc[:, "module_id"])
    assert is_float_dtype(expanded.loc[:, "support_weight"])
    assert is_float_dtype(expanded.loc[:, "top_score"])
    assert (
        expanded.loc[:, EXPANDED_SIGNALOME_ROW_KIND_COLUMN]
        == EXPANDED_SIGNALOME_ROW_KIND_SITE
    ).any()


def test_signalome_workflow_runs_with_scipy_clustering_engine() -> None:
    dataset = build_rat_l6_dataset(n_sites=260)
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                deterministic_max_selected_kinases=12,
                adaptive_ensemble_runs=12,
            ),
            activity_config=None,
        )
    )
    result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=build_signalome_config(
                substrate_support_cutoff=0.5,
                clustering_engine=SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
            ),
        )
    )

    assert not result.module_assignments.table.empty
    assert not result.signalome_modules.table.empty
    assert result.provenance is not None
    signalome_config = result.provenance.workflow_parameters["signalome_config"]
    assert (
        signalome_config["clustering"]["clustering_engine"]
        == SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL
    )
    scale_guard = result.provenance.workflow_parameters["scale_guard"]
    assert (
        scale_guard["clustering_engine"]
        == SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL
    )
    backend_diagnostics = scale_guard["backend_diagnostics"]
    assert isinstance(backend_diagnostics, dict)
    assert backend_diagnostics["uses_scipy"] is True
    score_semantics = result.provenance.workflow_parameters["signalome_score_semantics"]
    assert score_semantics["clustering_engine"] == (
        SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL
    )
    assert score_semantics["candidate_scoring_scope"] == (
        "candidate_module_count_evaluation_only"
    )
    assert (
        score_semantics["network_policy"]
        == signalome_config["output"]["network_policy"]
    )


def test_signalome_workflow_requires_explicit_dataset_site_metadata_protein_id() -> (
    None
):
    base_dataset = build_rat_l6_dataset(n_sites=260)
    dataset_without_protein = AnalysisReadyPhosphoDataset(
        phospho=base_dataset.phospho,
        site_metadata=base_dataset.site_metadata.drop(columns=["protein_id"]),
        sample_metadata=base_dataset.sample_metadata,
        total=base_dataset.total,
        organism=base_dataset.organism,
        intensity_scale_state=base_dataset.intensity_scale_state,
        processing_state=base_dataset.processing_state,
    )
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset_without_protein,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                deterministic_max_selected_kinases=12,
                adaptive_ensemble_runs=12,
            ),
            activity_config=None,
        )
    )

    with pytest.raises(
        WorkflowValidationError,
        match="site_metadata is missing required columns: protein_id",
    ):
        SignalomeWorkflow().run(
            SignalomeWorkflowRequest(
                kinase_result=kinase_result,
                config=build_signalome_config(substrate_support_cutoff=0.5),
            )
        )


def test_signalome_workflow_uses_explicit_dataset_protein_identity_when_present() -> (
    None
):
    base_dataset = build_rat_l6_dataset(n_sites=260)
    site_metadata = base_dataset.site_metadata.copy(deep=True)
    site_metadata.loc[:, "protein_id"] = [
        f"PROT_{position:05d}" for position in range(site_metadata.shape[0])
    ]
    dataset = AnalysisReadyPhosphoDataset(
        phospho=base_dataset.phospho,
        site_metadata=site_metadata,
        sample_metadata=base_dataset.sample_metadata,
        total=base_dataset.total,
        organism=base_dataset.organism,
        intensity_scale_state=base_dataset.intensity_scale_state,
        processing_state=base_dataset.processing_state,
    )
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                deterministic_max_selected_kinases=12,
                adaptive_ensemble_runs=12,
            ),
            activity_config=None,
        )
    )
    result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=build_signalome_config(substrate_support_cutoff=0.5),
        )
    )

    assignments = result.module_assignments.table
    assert not assignments.empty
    assert assignments.loc[:, "protein_id"].astype(str).str.startswith("PROT_").all()


def test_signalome_threshold_knobs_do_not_cross_couple_unrelated_outputs() -> None:
    dataset = build_rat_l6_dataset(n_sites=260)
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                deterministic_max_selected_kinases=12,
                adaptive_ensemble_runs=12,
            ),
            activity_config=None,
        )
    )

    baseline = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=build_signalome_config(
                substrate_support_cutoff=0.5,
                network_correlation_threshold=0.5,
            ),
        )
    )
    support_shifted = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=build_signalome_config(
                substrate_support_cutoff=0.65,
                network_correlation_threshold=0.5,
            ),
        )
    )
    network_shifted = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=build_signalome_config(
                substrate_support_cutoff=0.5,
                network_correlation_threshold=0.8,
            ),
        )
    )

    pdt.assert_frame_equal(
        baseline.kinase_network.edges,
        support_shifted.kinase_network.edges,
        check_dtype=False,
    )
    pdt.assert_frame_equal(
        baseline.module_assignments.table,
        network_shifted.module_assignments.table,
        check_dtype=False,
    )
    pdt.assert_frame_equal(
        baseline.signalome_modules.table,
        network_shifted.signalome_modules.table,
        check_dtype=False,
    )


def test_signalome_network_uses_rank_weighted_fusion_downstream_scores_when_available(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset = build_rat_l6_dataset(n_sites=260)
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                deterministic_max_selected_kinases=12,
                adaptive_ensemble_runs=12,
            ),
            activity_config=None,
        )
    )
    request = SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=build_signalome_config(substrate_support_cutoff=0.5),
    )

    rank_weighted_fusion_lane = SignalomeWorkflow().run(request)

    def _force_profile_lane(*, profile_scores, rank_weighted_fusion_scores):
        _ = rank_weighted_fusion_scores
        return profile_scores, "profile_scores"

    monkeypatch.setattr(
        signalome_interpreter,
        "select_downstream_score_matrix",
        _force_profile_lane,
    )
    profile_lane = SignalomeWorkflow().run(request)

    assert not rank_weighted_fusion_lane.kinase_network.edges.equals(
        profile_lane.kinase_network.edges
    )
    assert not rank_weighted_fusion_lane.module_assignments.table.equals(
        profile_lane.module_assignments.table
    )
    pdt.assert_series_equal(
        rank_weighted_fusion_lane.module_assignments.table.loc[:, "top_kinase"],
        profile_lane.module_assignments.table.loc[:, "top_kinase"],
        check_dtype=False,
    )


def test_signalome_workflow_accepts_sparse_missing_rank_weighted_fusion_score_rows() -> (
    None
):
    dataset = build_rat_l6_dataset(n_sites=260)
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                deterministic_max_selected_kinases=12,
                adaptive_ensemble_runs=12,
            ),
            activity_config=None,
        )
    )
    rank_weighted_fusion_scores = (
        kinase_result.scoring_result.rank_weighted_fusion_scores
    )
    assert rank_weighted_fusion_scores is not None
    sparse_rank_weighted_fusion_scores = rank_weighted_fusion_scores.copy(deep=True)
    sparse_rank_weighted_fusion_scores.iloc[:5, :] = float("nan")
    sparse_kinase_result = KinaseWorkflowResult(
        dataset=kinase_result.dataset,
        references=kinase_result.references,
        scoring_result=KinaseScoringResult(
            profile_scores=kinase_result.scoring_result.profile_scores,
            motif_scores=kinase_result.scoring_result.motif_scores,
            rank_weighted_fusion_scores=sparse_rank_weighted_fusion_scores,
            score_fusion_weights=kinase_result.scoring_result.score_fusion_weights,
        ),
        prediction_result=kinase_result.prediction_result,
        activity_result=kinase_result.activity_result,
    )

    result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=sparse_kinase_result,
            config=build_signalome_config(
                substrate_support_cutoff=0.5,
                network_correlation_threshold=0.2,
                score_preconditioning_policy="allow_and_report",
            ),
        )
    )

    assert not result.module_assignments.table.empty
    assert not result.signalome_modules.table.empty
    expected_input_rows = int(kinase_result.prediction_result.pred_mat.shape[0])
    expected_retained_rows = expected_input_rows - 5
    assert result.module_assignments.table.shape[0] == expected_retained_rows
    assert result.score_preconditioning_diagnostics.input_row_count == int(
        kinase_result.prediction_result.pred_mat.shape[0]
    )
    assert result.score_preconditioning_diagnostics.dropped_all_missing_row_count == 5
    assert result.score_preconditioning_diagnostics.retained_row_count == (
        expected_retained_rows
    )
    assert result.score_preconditioning_diagnostics.policy == "allow_and_report"


def test_signalome_workflow_rejects_sparse_missing_rank_weighted_fusion_rows_under_default_strict_policy() -> (
    None
):
    dataset = build_rat_l6_dataset(n_sites=260)
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                deterministic_max_selected_kinases=12,
                adaptive_ensemble_runs=12,
            ),
            activity_config=None,
        )
    )
    rank_weighted_fusion_scores = (
        kinase_result.scoring_result.rank_weighted_fusion_scores
    )
    assert rank_weighted_fusion_scores is not None
    sparse_rank_weighted_fusion_scores = rank_weighted_fusion_scores.copy(deep=True)
    sparse_rank_weighted_fusion_scores.iloc[:5, :] = float("nan")
    sparse_kinase_result = KinaseWorkflowResult(
        dataset=kinase_result.dataset,
        references=kinase_result.references,
        scoring_result=KinaseScoringResult(
            profile_scores=kinase_result.scoring_result.profile_scores,
            motif_scores=kinase_result.scoring_result.motif_scores,
            rank_weighted_fusion_scores=sparse_rank_weighted_fusion_scores,
            score_fusion_weights=kinase_result.scoring_result.score_fusion_weights,
        ),
        prediction_result=kinase_result.prediction_result,
        activity_result=kinase_result.activity_result,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SignalomeWorkflow().run(
            SignalomeWorkflowRequest(
                kinase_result=sparse_kinase_result,
                config=build_signalome_config(
                    substrate_support_cutoff=0.5,
                    network_correlation_threshold=0.2,
                ),
            )
        )

    error = exc_info.value
    message = str(error)
    assert error.seam == SIGNALOME_INTERPRETER_SCORE_PRECONDITIONING_SEAM
    assert error.details["dropped_all_missing_row_count"] == 5
    assert error.details["score_preconditioning_policy"] == "error_on_drop"
    assert "dropped_all_missing_row_count=5" in message
    assert "score_preconditioning_policy=error_on_drop" in message


def test_signalome_l6_provenance_matches_golden_contract() -> None:
    dataset = build_rat_l6_dataset(n_sites=260)
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                deterministic_max_selected_kinases=12,
                adaptive_ensemble_runs=12,
            ),
            activity_config=None,
        )
    )
    result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=build_signalome_config(substrate_support_cutoff=0.5),
        )
    )
    provenance = result.provenance
    assert provenance is not None
    golden = load_signalome_l6_provenance_golden()

    assert provenance.workflow_name == golden["workflow_name"]
    assert provenance.random_state is None
    assert provenance.random_seed_policy is None
    assert (
        provenance.workflow_parameters["signalome_config"]
        == (golden["workflow_parameters"]["signalome_config"])
    )
    assert (
        provenance.workflow_parameters["score_preconditioning_diagnostics"]
        == (golden["workflow_parameters"]["score_preconditioning_diagnostics"])
    )
    upstream = provenance.workflow_parameters["upstream_kinase_provenance"]
    assert isinstance(upstream, dict)
    assert (
        upstream["workflow_name"]
        == golden["upstream_kinase_provenance"]["workflow_name"]
    )
    assert (
        upstream["random_state"] == golden["upstream_kinase_provenance"]["random_state"]
    )
    assert provenance.reference is not None
    assert provenance.reference.source_type == golden["reference"]["source_type"]
    assert provenance.reference.organism == golden["reference"]["organism"]
    assert provenance.reference.bundle_id == golden["reference"]["bundle_id"]

    _assert_expected_fingerprint_map(
        observed=_fingerprints_by_name(provenance.input_tables),
        expected=golden["input_tables"],
    )
    _assert_expected_fingerprint_map(
        observed=_fingerprints_by_name(provenance.output_tables),
        expected=golden["output_tables"],
    )
    _assert_expected_fingerprint_map(
        observed=_fingerprints_by_name(provenance.reference.table_fingerprints),
        expected=golden["reference"]["table_fingerprints"],
    )
    assert [
        {"id": item.id.value, "name": item.name, "version": item.version}
        for item in provenance.scientific_policies
    ] == golden["scientific_policies"]
