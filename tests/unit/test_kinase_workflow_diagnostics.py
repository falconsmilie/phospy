from __future__ import annotations

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    KinaseWorkflow,
)
from phospy.api import (
    DatasetBuildRequest,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferencePreset,
)
from phospy.errors import WorkflowBoundaryError
from phospy.science.activities.models import KinaseActivityInputs, PredMatOverlapSummary
from phospy.science.prediction.models import KinasePredictionResult
from phospy.science.references.resolution import ReferenceResolver
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseActivityExecutionConfig,
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.executor import KinaseWorkflowExecutor
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.rewrite_fixture_data import build_rat_l6_dataset
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)


def _dataset(
    *,
    site_ids: list[str],
    sample_names: list[str],
) -> AnalysisReadyPhosphoDataset:
    site_index = site_key_index_from_display_ids(site_ids)
    phospho = pd.DataFrame(
        {
            sample: [
                float((index + 1) * (sample_position + 1))
                for index in range(len(site_ids))
            ]
            for sample_position, sample in enumerate(sample_names)
        },
        index=site_index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": site_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": [site.split(";", 1)[0] for site in site_ids],
            "protein_id": [site.split(";", 1)[0] for site in site_ids],
            "site": [site.split(";")[1] for site in site_ids],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in [site.split(";")[1] for site in site_ids]
            ],
            "localisation_confidence": [0.95 for _ in site_ids],
        },
        index=site_index.copy(),
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


def _bundle(kinase_substrate_map: pd.DataFrame) -> ReferenceBundle:
    unique_sites = pd.Index(
        kinase_substrate_map.loc[:, "substrate_site"].astype(str).str.strip().unique(),
        name="site_id",
    )
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=kinase_substrate_map,
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


def _resolved_request(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    references: ReferenceBundle,
    min_substrates: int = 2,
    top_k: int = 3,
    deterministic_max_selected_kinases: int = 2,
    adaptive_ensemble_runs: int = 2,
    threshold: float = 0.7,
    activity_min_substrates: int = 2,
    activity_top_n_substrates: int = 3,
) -> ResolvedKinaseWorkflowRequest:
    projected_map, projected_sequences, scoring_site_index, site_identity_map = (
        _project_reference_inputs(dataset=dataset, references=references)
    )
    execution_config = ResolvedKinaseExecutionConfig(
        scoring_min_substrates=int(min_substrates),
        include_diagnostic_scoring_tables=False,
        profile_missing_value_strategy="strict",
        prediction_top_k=int(top_k),
        prediction_deterministic_max_selected_kinases=int(
            deterministic_max_selected_kinases
        ),
        prediction_adaptive_ensemble_runs=int(adaptive_ensemble_runs),
        prediction_mode="deterministic_ranking",
        prediction_adaptive_policy="stable",
        prediction_n_iterations=5,
        prediction_random_state=None,
        activity=ResolvedKinaseActivityExecutionConfig(
            method="simplified_weighted_substrate_activity",
            threshold=float(threshold),
            min_substrates=int(activity_min_substrates),
            top_n_substrates=int(activity_top_n_substrates),
            ksea_min_substrates=5,
            ksea_evidence_threshold=float(threshold),
            ksea_p_value_method="normal_approximation",
            ksea_adjust_p_values=True,
        ),
    )
    return ResolvedKinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        kinase_substrate_map=projected_map,
        site_sequences=projected_sequences,
        site_identity_map=site_identity_map,
        scoring_site_index=scoring_site_index,
        activity_phospho_matrix=dataset.phospho.loc[scoring_site_index, :].copy(
            deep=True
        ),
        execution_config=execution_config,
    )


def _project_reference_inputs(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    references: ReferenceBundle,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Index, pd.DataFrame]:
    site_metadata = dataset.site_metadata
    display_to_site_key = {
        str(display_id): str(site_key)
        for site_key, display_id in site_metadata.loc[
            :, ["site_key", "display_id"]
        ].itertuples(index=False)
    }
    scoring_display_ids = [
        str(display_id)
        for display_id in site_metadata.loc[:, "display_id"].astype(str).tolist()
        if str(display_id) in set(references.site_sequences.index.astype(str))
    ]
    scoring_site_index = pd.Index(
        [display_to_site_key[display_id] for display_id in scoring_display_ids],
        name=dataset.phospho.index.name,
    )
    site_sequences = references.site_sequences.reindex(scoring_display_ids).copy(
        deep=True
    )
    site_sequences.loc[:, "display_id"] = scoring_display_ids
    site_sequences.index = scoring_site_index.copy()
    projected_rows: list[dict[str, str]] = []
    for kinase, display_id in references.kinase_substrate_map.loc[
        :, ["kinase", "substrate_site"]
    ].itertuples(index=False):
        site_key = display_to_site_key.get(str(display_id))
        if site_key is None:
            continue
        projected_rows.append(
            {
                "kinase": str(kinase),
                "substrate_site": site_key,
                "display_id": str(display_id),
            }
        )
    site_identity_map = site_metadata.loc[
        scoring_site_index, ["site_key", "display_id"]
    ].copy(deep=True)
    return (
        pd.DataFrame(projected_rows),
        site_sequences,
        scoring_site_index,
        site_identity_map,
    )


def test_interpreter_resolves_execution_config_defaults_for_executor() -> None:
    dataset = _dataset(
        site_ids=["MAPK14;Y182;", "GSK3B;S9;"],
        sample_names=["sample_a", "sample_b"],
    )
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=_bundle(
            pd.DataFrame(
                {
                    "kinase": ["MAP2K6", "MAP2K6"],
                    "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
                }
            )
        ),
    )

    interpreted = KinaseWorkflowInterpreter().run(request)

    assert interpreted.execution_config.scoring_min_substrates == 2
    assert interpreted.execution_config.include_diagnostic_scoring_tables is False
    assert interpreted.execution_config.profile_missing_value_strategy == "strict"
    assert interpreted.execution_config.prediction_mode == "deterministic_ranking"
    assert interpreted.execution_config.prediction_top_k == 30
    assert (
        interpreted.execution_config.prediction_deterministic_max_selected_kinases == 10
    )
    assert interpreted.execution_config.prediction_adaptive_ensemble_runs == 10
    assert interpreted.execution_config.prediction_adaptive_policy == "stable"
    assert interpreted.execution_config.prediction_n_iterations == 5
    assert interpreted.execution_config.prediction_random_state is None
    assert interpreted.execution_config.activity is not None
    assert (
        interpreted.execution_config.activity.method
        == "simplified_weighted_substrate_activity"
    )
    assert interpreted.execution_config.activity.threshold == pytest.approx(0.6)
    assert interpreted.execution_config.activity.min_substrates == 3
    assert interpreted.execution_config.activity.top_n_substrates == 20
    assert interpreted.execution_config.activity.ksea_min_substrates == 5
    assert (
        interpreted.execution_config.activity.ksea_evidence_threshold
        == pytest.approx(0.6)
    )


def test_interpreter_merges_dataset_site_sequences_without_mutating_references() -> (
    None
):
    dataset = _dataset(
        site_ids=["MAPK14;Y182;", "GSK3B;S9;", "EXTRA;S1;"],
        sample_names=["sample_a", "sample_b"],
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6"],
                "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31, "A" * 31]},
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
        ),
    )
    original_reference_sequences = references.site_sequences.copy(deep=True)

    interpreted = KinaseWorkflowInterpreter().run(
        KinaseWorkflowRequest(dataset=dataset, references=references)
    )

    assert set(interpreted.site_sequences.index.astype(str)) == set(
        dataset.phospho.index.astype(str)
    )
    pd.testing.assert_frame_equal(
        references.site_sequences,
        original_reference_sequences,
    )
    assert interpreted.site_sequence_merge_diagnostics["dataset_sequences_added"] == 1


def test_scoring_results_are_reference_input_form_invariant_for_equivalent_content() -> (
    None
):
    dataset = build_rat_l6_dataset(n_sites=220)
    explicit_bundle = ReferenceResolver().run(
        ReferencePreset.AUTO,
        dataset_organism=dataset.organism,
    )
    preset_request = KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=6,
            deterministic_max_selected_kinases=8,
            adaptive_ensemble_runs=8,
        ),
        activity_config=None,
    )
    bundle_request = KinaseWorkflowRequest(
        dataset=dataset,
        references=explicit_bundle,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=6,
            deterministic_max_selected_kinases=8,
            adaptive_ensemble_runs=8,
        ),
        activity_config=None,
    )
    from_preset = KinaseWorkflow().run(preset_request)
    from_bundle = KinaseWorkflow().run(bundle_request)

    pd.testing.assert_frame_equal(
        from_preset.scoring_result.profile_scores,
        from_bundle.scoring_result.profile_scores,
        check_dtype=False,
    )
    assert from_preset.scoring_result.rank_weighted_fusion_scores is not None
    assert from_bundle.scoring_result.rank_weighted_fusion_scores is not None
    pd.testing.assert_frame_equal(
        from_preset.scoring_result.rank_weighted_fusion_scores,
        from_bundle.scoring_result.rank_weighted_fusion_scores,
        check_dtype=False,
    )
    pd.testing.assert_frame_equal(
        from_preset.prediction_result.pred_mat,
        from_bundle.prediction_result.pred_mat,
        check_dtype=False,
    )


def test_workflow_uses_execution_time_merged_site_sequences() -> None:
    dataset = _dataset(
        site_ids=["MAPK14;Y182;", "GSK3B;S9;", "EXTRA;S1;"],
        sample_names=["sample_a", "sample_b"],
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6"],
                "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31, "A" * 31]},
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
        ),
    )
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=2,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=None,
        )
    )

    expected_site_keys = dataset.phospho.index.astype(str).tolist()
    assert list(result.scoring_result.profile_scores.index) == expected_site_keys
    assert list(result.prediction_result.pred_mat.index) == expected_site_keys


def test_boundary_error_reports_unusable_reference_coverage_counts() -> None:
    dataset = _dataset(
        site_ids=["MAPK14;Y182;", "AKT1;T308;"],
        sample_names=["sample_a", "sample_b"],
    )
    references = _bundle(
        pd.DataFrame(
            {
                "kinase": ["MAP2K6", "AKT1"],
                "substrate_site": ["GSK3B;S9;", "STAT3;Y705;"],
            }
        )
    )
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
        ),
        activity_config=None,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflow().run(request)

    error = exc_info.value
    message = str(error)
    assert error.seam == "kinase.interpreter.reference_coverage"
    assert error.next_action is not None
    assert error.details["dataset_sites"] == 2
    assert error.details["reference_sites"] == 2
    assert error.details["overlap_sites"] == 0
    assert error.details["scoring_config_min_substrates"] == 2
    assert "seam=kinase.interpreter.reference_coverage" in message
    assert "dataset_sites=2" in message
    assert "reference_sites=2" in message
    assert "overlap_sites=0" in message
    assert "scoring_config_min_substrates=2" in message
    assert "next_action=" in message


def test_boundary_error_reports_empty_eligible_kinase_set_counts() -> None:
    dataset = _dataset(
        site_ids=["MAPK14;Y182;", "AKT1;T308;"],
        sample_names=["sample_a", "sample_b"],
    )
    references = _bundle(
        pd.DataFrame(
            {
                "kinase": ["MAP2K6", "AKT1"],
                "substrate_site": ["MAPK14;Y182;", "AKT1;T308;"],
            }
        )
    )
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=3,
            deterministic_max_selected_kinases=5,
            adaptive_ensemble_runs=5,
        ),
        activity_config=None,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflow().run(request)

    message = str(exc_info.value)
    assert "seam=kinase.interpreter.eligible_kinases" in message
    assert "reference_kinases=2" in message
    assert "kinases_with_overlap=2" in message
    assert "eligible_kinases=0" in message
    assert "max_quantified_sites_per_kinase=1" in message
    assert "scoring_config_min_substrates=2" in message
    assert "prediction_config_deterministic_max_selected_kinases=5" in message
    assert "prediction_config_adaptive_ensemble_runs=5" in message


def test_default_scoring_floor_rejects_single_substrate_kinase_profiles() -> None:
    dataset = _dataset(
        site_ids=["MAPK14;Y182;"],
        sample_names=["sample_a", "sample_b"],
    )
    references = _bundle(
        pd.DataFrame(
            {
                "kinase": ["MAP2K6"],
                "substrate_site": ["MAPK14;Y182;"],
            }
        )
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflow().run(
            KinaseWorkflowRequest(
                dataset=dataset,
                references=references,
            )
        )

    message = str(exc_info.value)
    assert "seam=kinase.interpreter.eligible_kinases" in message
    assert "scoring_config_min_substrates=2" in message
    assert "max_quantified_sites_per_kinase=1" in message


def test_eligibility_report_full_reference_overlap_counts() -> None:
    dataset = _dataset(
        site_ids=["MAPK14;Y182;", "GSK3B;S9;"],
        sample_names=["sample_a", "sample_b"],
    )
    references = _bundle(
        pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6", "AKT1", "AKT1"],
                "substrate_site": [
                    "MAPK14;Y182;",
                    "GSK3B;S9;",
                    "MAPK14;Y182;",
                    "GSK3B;S9;",
                ],
            }
        )
    )

    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=2,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=None,
        )
    )

    assert result.eligibility_report is not None
    assert result.eligibility_report.total_dataset_sites == 2
    assert result.eligibility_report.sequence_complete_sites == 2
    assert result.eligibility_report.localisation_eligible_sites is None
    assert result.eligibility_report.reference_overlap_sites == 2
    assert result.eligibility_report.excluded_no_reference_match == 0
    assert result.eligibility_report.excluded_low_localisation is None
    assert result.eligibility_report.eligible_kinases == 2
    assert result.eligibility_report.excluded_kinases_below_min_substrates == 0


def test_eligibility_report_partial_overlap_and_kinase_threshold_shortfall_counts() -> (
    None
):
    dataset = _dataset(
        site_ids=["MAPK14;Y182;", "GSK3B;S9;", "EXTRA;S1;"],
        sample_names=["sample_a", "sample_b"],
    )
    references = _bundle(
        pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6", "AKT1", "ERK1"],
                "substrate_site": [
                    "MAPK14;Y182;",
                    "GSK3B;S9;",
                    "MAPK14;Y182;",
                    "OFFSITE;T1;",
                ],
            }
        )
    )

    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=2,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=None,
        )
    )

    assert result.eligibility_report is not None
    assert result.eligibility_report.total_dataset_sites == 3
    assert result.eligibility_report.sequence_complete_sites == 3
    assert result.eligibility_report.reference_overlap_sites == 2
    assert result.eligibility_report.excluded_no_reference_match == 1
    assert result.eligibility_report.eligible_kinases == 1
    assert result.eligibility_report.excluded_kinases_below_min_substrates == 1


def test_eligibility_report_includes_localisation_counts_when_policy_available() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [2.0, 4.0]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "protein_id": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "S9"]
            ],
            "localisation_confidence": [0.98, 0.97],
        },
        index=phospho.index.copy(),
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=_bundle(
                pd.DataFrame(
                    {
                        "kinase": ["MAP2K6", "MAP2K6"],
                        "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
                    }
                )
            ),
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=2,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=None,
        )
    )

    assert result.eligibility_report is not None
    assert result.eligibility_report.localisation_eligible_sites == 2
    assert result.eligibility_report.excluded_low_localisation == 0


def test_boundary_error_reports_prediction_ensemble_collapse_counts() -> None:
    dataset = _dataset(
        site_ids=["MAPK14;Y182;", "GSK3B;S9;"],
        sample_names=["sample_a"],
    )
    references = _bundle(
        pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6"],
                "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
            }
        )
    )
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=3,
            adaptive_ensemble_runs=3,
        ),
        activity_config=None,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflow().run(request)

    message = str(exc_info.value)
    assert "seam=kinase.executor.prediction_ensemble" in message
    assert "eligible_kinases=1" in message
    assert "ranked_kinases=0" in message
    assert "prediction_config_deterministic_max_selected_kinases=3" in message
    assert "prediction_config_top_k=2" in message
    assert "dataset_samples=1" in message
    assert "dataset.phospho" in message
    assert "scoring_config.min_substrates" in message


def test_boundary_error_reports_activity_overlap_edge_case() -> None:
    dataset = _dataset(
        site_ids=["MAPK14;Y182;"],
        sample_names=["sample_a", "sample_b"],
    )
    references = _bundle(
        pd.DataFrame(
            {
                "kinase": ["MAP2K6"],
                "substrate_site": ["MAPK14;Y182;"],
            }
        )
    )
    request = _resolved_request(
        dataset=dataset,
        references=references,
        min_substrates=2,
        top_k=3,
        deterministic_max_selected_kinases=2,
        adaptive_ensemble_runs=2,
        threshold=0.7,
    )
    prediction_result = KinasePredictionResult(
        pred_mat=pd.DataFrame(
            {"MAP2K6": [0.8]},
            index=pd.Index(["OTHER;S1;"], name="site_id"),
        ),
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflowExecutor()._run_activity_stage(
            request=request,
            config=request.execution_config,
            prediction_result=prediction_result,
        )

    message = str(exc_info.value)
    assert "seam=kinase.activity.input_overlap" in message
    assert "overlap_sites=0" in message
    assert "pred_mat_sites=1" in message
    assert "phospho_sites=1" in message
    assert "min_overlap=1" in message
    assert "min_fraction=0.5" in message
    assert "prediction_result.pred_mat" in message
    assert "dataset.phospho" in message


def test_activity_inputs_reject_malformed_pred_mat_with_workflow_boundary_error() -> (
    None
):
    with pytest.raises(
        WorkflowBoundaryError, match="seam=kinase.activity.input_schema"
    ):
        KinaseActivityInputs(
            pred_mat=pd.DataFrame(
                {"MAP2K6": [1.5]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
            phospho_matrix=pd.DataFrame(
                {"sample_a": [1.0]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            ),
            threshold=0.5,
            min_substrates=1,
            top_n_substrates=1,
            overlap_summary=PredMatOverlapSummary(
                overlap_count=1,
                pred_mat_rows=1,
                phospho_rows=1,
            ),
        )


def test_boundary_error_reports_no_activity_candidates_after_filtering() -> None:
    dataset = _dataset(
        site_ids=["MAPK14;Y182;", "GSK3B;S9;"],
        sample_names=["sample_a", "sample_b"],
    )
    references = _bundle(
        pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6"],
                "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
            }
        )
    )
    request = _resolved_request(
        dataset=dataset,
        references=references,
        threshold=0.95,
        activity_min_substrates=3,
        activity_top_n_substrates=2,
    )
    prediction_result = KinasePredictionResult(
        pred_mat=pd.DataFrame(
            {"MAP2K6": [0.8, 0.7]},
            index=dataset.phospho.index.copy(),
        ),
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflowExecutor()._run_activity_stage(
            request=request,
            config=request.execution_config,
            prediction_result=prediction_result,
        )

    message = str(exc_info.value)
    assert "seam=kinase.activity.valid_candidates" in message
    assert "weighted_activity_kinases=0" in message
    assert "thresholded_mean_activity_kinases=0" in message
    assert "activity_config_threshold=0.95" in message
    assert "activity_config_min_substrates=3" in message
    assert "activity_config_top_n_substrates=2" in message


def test_activity_stage_returns_weighted_thresholded_mean_and_target_outputs() -> None:
    dataset = _dataset(
        site_ids=["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"],
        sample_names=["sample_a", "sample_b"],
    )
    references = _bundle(
        pd.DataFrame(
            {
                "kinase": ["MAP2K6", "AKT1"],
                "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
            }
        )
    )
    request = _resolved_request(
        dataset=dataset,
        references=references,
        threshold=0.3,
        activity_min_substrates=2,
        activity_top_n_substrates=2,
    )
    prediction_result = KinasePredictionResult(
        pred_mat=pd.DataFrame(
            {
                "MAP2K6": [0.9, 0.5, 0.0],
                "AKT1": [0.4, 0.2, 0.1],
            },
            index=dataset.phospho.index.copy(),
        ),
    )

    result = KinaseWorkflowExecutor()._run_activity_stage(
        request=request,
        config=request.execution_config,
        prediction_result=prediction_result,
    )
    assert result is not None
    assert result.weighted_activity.at["MAP2K6", "sample_a"] == pytest.approx(
        (1.0 * 0.9 + 2.0 * 0.5) / (0.9 + 0.5)
    )
    assert result.weighted_activity.at["AKT1", "sample_b"] == pytest.approx(
        (2.0 * 0.4 + 4.0 * 0.2) / (0.4 + 0.2)
    )
    assert list(result.thresholded_substrate_mean_activity.index) == ["MAP2K6"]
    assert result.thresholded_substrate_mean_activity.at[
        "MAP2K6", "sample_a"
    ] == pytest.approx(1.5)
    assert result.thresholded_substrate_counts.to_dict() == {"MAP2K6": 2}
    assert result.target_counts.to_dict() == {"MAP2K6": 2, "AKT1": 1}
    assert {"site_id", "site_key", "display_id", "kinase", "score"} <= set(
        result.target_table.columns
    )
    assert result.target_table.loc[:, "site_key"].notna().all()
    assert set(result.target_table.loc[:, "display_id"].astype(str)) == {
        "MAPK14;Y182;",
        "GSK3B;S9;",
    }
    assert int(result.target_table.shape[0]) == 3


def test_weighted_activity_is_stable_under_zero_padding_matrix_growth() -> None:
    dataset = _dataset(
        site_ids=[
            "MAPK14;Y182;",
            "GSK3B;S9;",
            "AKT1;T308;",
            "RPS6KB1;T412;",
            "EIF4EBP1;T70;",
        ],
        sample_names=["sample_a", "sample_b"],
    )
    references = _bundle(
        pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6"],
                "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
            }
        )
    )
    request = _resolved_request(dataset=dataset, references=references, threshold=0.5)
    compact_predictions = KinasePredictionResult(
        pred_mat=pd.DataFrame(
            {"MAP2K6": [0.8, 0.4]},
            index=pd.Index(
                site_key_index_from_display_ids(["MAPK14;Y182;", "GSK3B;S9;"]),
                name=dataset.phospho.index.name,
            ),
        ),
    )
    padded_predictions = KinasePredictionResult(
        pred_mat=pd.DataFrame(
            {"MAP2K6": [0.8, 0.4, 0.0, 0.0, 0.0]},
            index=dataset.phospho.index.copy(),
        ),
    )

    compact = KinaseWorkflowExecutor()._run_activity_stage(
        request=request,
        config=request.execution_config,
        prediction_result=compact_predictions,
    )
    padded = KinaseWorkflowExecutor()._run_activity_stage(
        request=request,
        config=request.execution_config,
        prediction_result=padded_predictions,
    )

    assert compact is not None
    assert padded is not None
    assert compact.weighted_activity.at["MAP2K6", "sample_a"] == pytest.approx(
        padded.weighted_activity.at["MAP2K6", "sample_a"]
    )
    assert compact.weighted_activity.at["MAP2K6", "sample_b"] == pytest.approx(
        padded.weighted_activity.at["MAP2K6", "sample_b"]
    )
