from __future__ import annotations

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyPhosphoDataset,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from phospy.errors import WorkflowBoundaryError
from phospy.prediction.models import KinasePredictionResult
from phospy.transformations.models import TransformationState
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseActivityExecutionConfig,
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.executor import KinaseWorkflowExecutor
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter


def _dataset(
    *,
    site_ids: list[str],
    sample_names: list[str],
) -> AnalysisReadyPhosphoDataset:
    phospho = pd.DataFrame(
        {
            sample: [
                float((index + 1) * (sample_position + 1))
                for index in range(len(site_ids))
            ]
            for sample_position, sample in enumerate(sample_names)
        },
        index=site_ids,
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": [site.split(";", 1)[0] for site in site_ids],
            "site": [site.split(";")[1] for site in site_ids],
            "site_sequence": ["A" * 31 for _ in site_ids],
        },
        index=site_ids,
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        transformation_state=TransformationState.raw(has_total_matrix=False),
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
            {"site_sequence": ["A" * 31 for _ in unique_sites]},
            index=unique_sites,
        ),
    )


def _resolved_request(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    references: ReferenceBundle,
    min_substrates: int = 2,
    top_k: int = 3,
    ensemble_size: int = 2,
    threshold: float = 0.7,
    activity_min_substrates: int = 2,
    activity_top_n_substrates: int = 3,
) -> ResolvedKinaseWorkflowRequest:
    scoring_site_index = dataset.phospho.index.intersection(
        references.site_sequences.index
    )
    execution_config = ResolvedKinaseExecutionConfig(
        scoring_min_substrates=int(min_substrates),
        include_diagnostic_scoring_tables=False,
        profile_missing_value_strategy="strict",
        prediction_top_k=int(top_k),
        prediction_ensemble_size=int(ensemble_size),
        prediction_mode="deterministic_ranking",
        prediction_adaptive_policy="stable",
        prediction_n_iterations=5,
        prediction_random_state=None,
        activity=ResolvedKinaseActivityExecutionConfig(
            threshold=float(threshold),
            min_substrates=int(activity_min_substrates),
            top_n_substrates=int(activity_top_n_substrates),
        ),
    )
    return ResolvedKinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        kinase_substrate_map=references.kinase_substrate_map,
        site_sequences=references.site_sequences,
        scoring_site_index=scoring_site_index,
        activity_phospho_matrix=dataset.phospho.loc[scoring_site_index, :].copy(
            deep=True
        ),
        uses_bundled_reference=False,
        execution_config=execution_config,
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
    assert interpreted.execution_config.prediction_ensemble_size == 10
    assert interpreted.execution_config.prediction_adaptive_policy == "stable"
    assert interpreted.execution_config.prediction_n_iterations == 5
    assert interpreted.execution_config.prediction_random_state is None
    assert interpreted.execution_config.activity is not None
    assert interpreted.execution_config.activity.threshold == pytest.approx(0.6)
    assert interpreted.execution_config.activity.min_substrates == 3
    assert interpreted.execution_config.activity.top_n_substrates == 20


def test_workflow_limits_scoring_to_sites_with_reference_sequences() -> None:
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
            prediction_config=KinasePredictionConfig(top_k=2, ensemble_size=2),
            activity_config=None,
        )
    )

    assert list(result.scoring_result.profile_scores.index) == [
        "MAPK14;Y182;",
        "GSK3B;S9;",
    ]
    assert list(result.prediction_result.pred_mat.index) == [
        "MAPK14;Y182;",
        "GSK3B;S9;",
    ]


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
        prediction_config=KinasePredictionConfig(top_k=2, ensemble_size=2),
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
        prediction_config=KinasePredictionConfig(top_k=3, ensemble_size=5),
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
    assert "prediction_config_ensemble_size=5" in message


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
        prediction_config=KinasePredictionConfig(top_k=2, ensemble_size=3),
        activity_config=None,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflow().run(request)

    message = str(exc_info.value)
    assert "seam=kinase.executor.prediction_ensemble" in message
    assert "eligible_kinases=1" in message
    assert "ranked_kinases=0" in message
    assert "prediction_config_ensemble_size=3" in message
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
        ensemble_size=2,
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
    assert "ksea_kinases=0" in message
    assert "activity_config_threshold=0.95" in message
    assert "activity_config_min_substrates=3" in message
    assert "activity_config_top_n_substrates=2" in message


def test_activity_stage_returns_weighted_ksea_and_target_outputs() -> None:
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
    assert list(result.ksea_scores.index) == ["MAP2K6"]
    assert result.ksea_scores.at["MAP2K6", "sample_a"] == pytest.approx(1.5)
    assert result.ksea_counts.to_dict() == {"MAP2K6": 2}
    assert result.target_counts.to_dict() == {"MAP2K6": 2, "AKT1": 1}
    assert set(result.target_table.columns) == {"site_id", "kinase", "score"}
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
                ["MAPK14;Y182;", "GSK3B;S9;"],
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
