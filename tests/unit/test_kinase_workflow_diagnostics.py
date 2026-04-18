from __future__ import annotations

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyPhosphoDataset,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from phospy.errors import WorkflowBoundaryError
from phospy.prediction.models import KinasePredictionResult
from phospy.workflows.kinase.contracts import ResolvedKinaseWorkflowRequest
from phospy.workflows.kinase.executor import KinaseWorkflowExecutor


def _dataset(
    *,
    site_ids: list[str],
    sample_names: list[str],
) -> AnalysisReadyPhosphoDataset:
    phospho = pd.DataFrame(
        {
            sample: [float(index + 1) for index in range(len(site_ids))]
            for sample in sample_names
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
    min_substrates: int = 1,
    top_k: int = 3,
    ensemble_size: int = 2,
    threshold: float = 0.7,
) -> ResolvedKinaseWorkflowRequest:
    return ResolvedKinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        kinase_substrate_map=references.kinase_substrate_map,
        site_sequences=references.site_sequences,
        scoring_site_index=dataset.phospho.index.copy(),
        scoring_config=KinaseScoringConfig(min_substrates=min_substrates),
        prediction_config=KinasePredictionConfig(
            top_k=top_k,
            ensemble_size=ensemble_size,
        ),
        activity_config=KinaseActivityConfig(enabled=True, threshold=threshold),
    )


def test_workflow_limits_scoring_to_sites_with_reference_sequences() -> None:
    dataset = _dataset(
        site_ids=["MAPK14;Y182;", "EXTRA;S1;"],
        sample_names=["sample_a", "sample_b"],
    )
    dataset.phospho.loc["MAPK14;Y182;", "sample_b"] = 2.0
    dataset.phospho.loc["EXTRA;S1;", "sample_b"] = 4.0
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6"],
                "substrate_site": ["MAPK14;Y182;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
    )
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(min_substrates=1),
            prediction_config=KinasePredictionConfig(top_k=2, ensemble_size=2),
            activity_config=None,
        )
    )

    assert list(result.scoring_result.profile_scores.index) == ["MAPK14;Y182;"]
    assert list(result.prediction_result.pred_mat.index) == ["MAPK14;Y182;"]


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
        scoring_config=KinaseScoringConfig(min_substrates=1),
        prediction_config=KinasePredictionConfig(top_k=2, ensemble_size=2),
        activity_config=None,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflow().run(request)

    message = str(exc_info.value)
    assert "seam=kinase.interpreter.reference_coverage" in message
    assert "dataset_sites=2" in message
    assert "reference_sites=2" in message
    assert "overlap_sites=0" in message
    assert "scoring_config_min_substrates=1" in message
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


def test_boundary_error_reports_prediction_ensemble_collapse_counts() -> None:
    dataset = _dataset(
        site_ids=["MAPK14;Y182;"],
        sample_names=["sample_a"],
    )
    references = _bundle(
        pd.DataFrame(
            {
                "kinase": ["MAP2K6"],
                "substrate_site": ["MAPK14;Y182;"],
            }
        )
    )
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(min_substrates=1),
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


def test_boundary_error_reports_activity_support_edge_case() -> None:
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
        min_substrates=1,
        top_k=3,
        ensemble_size=2,
        threshold=0.7,
    )
    prediction_result = KinasePredictionResult(
        pred_mat=pd.DataFrame(
            {"MAP2K6": [0.0]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflowExecutor()._run_activity_stage(
            request=request,
            prediction_result=prediction_result,
        )

    message = str(exc_info.value)
    assert "seam=kinase.executor.activity_support" in message
    assert "activity_kinases=1" in message
    assert "total_positive_predictions=0" in message
    assert "prediction_config_top_k=3" in message
    assert "scoring_config_min_substrates=1" in message
    assert "prediction_config.top_k" in message
    assert "scoring_config.min_substrates" in message
    assert "activity_config.threshold" not in message


def test_activity_score_uses_only_positive_supported_predictions() -> None:
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
    request = _resolved_request(dataset=dataset, references=references, threshold=0.6)
    prediction_result = KinasePredictionResult(
        pred_mat=pd.DataFrame(
            {
                "MAP2K6": [0.9, 0.0, 0.0],
                "AKT1": [0.4, 0.2, 0.0],
            },
            index=dataset.phospho.index.copy(),
        ),
    )

    result = KinaseWorkflowExecutor()._run_activity_stage(
        request=request,
        prediction_result=prediction_result,
    )
    assert result is not None
    scores = result.activity_scores
    assert scores.at["MAP2K6", "activity_score"] == pytest.approx(0.9)
    assert scores.at["MAP2K6", "n_predicted_sites"] == 1
    assert scores.at["AKT1", "activity_score"] == pytest.approx(0.3)
    assert scores.at["AKT1", "n_predicted_sites"] == 2


def test_activity_score_is_stable_under_zero_padding_matrix_growth() -> None:
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
        prediction_result=compact_predictions,
    )
    padded = KinaseWorkflowExecutor()._run_activity_stage(
        request=request,
        prediction_result=padded_predictions,
    )

    assert compact is not None
    assert padded is not None
    assert compact.activity_scores.at["MAP2K6", "activity_score"] == pytest.approx(
        padded.activity_scores.at["MAP2K6", "activity_score"]
    )
    assert compact.activity_scores.at["MAP2K6", "weighted_signal"] == pytest.approx(
        padded.activity_scores.at["MAP2K6", "weighted_signal"]
    )
    assert compact.activity_scores.at["MAP2K6", "n_predicted_sites"] == 2
    assert padded.activity_scores.at["MAP2K6", "n_predicted_sites"] == 2
