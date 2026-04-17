from __future__ import annotations

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyPhosphoDataset,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    Organism,
    ReferenceBundle,
    SimpleKinaseWorkflowRequest,
)
from phospy.errors import WorkflowBoundaryError
from phospy.prediction.models import KinasePredictionResult
from phospy.workflows.kinase.contracts import ResolvedKinaseWorkflowRequest
from phospy.workflows.kinase.executor import SimpleKinaseWorkflowExecutor


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
    request = SimpleKinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(min_substrates=1),
        prediction_config=KinasePredictionConfig(top_k=2, ensemble_size=2),
        activity_config=None,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflow().run(request)

    message = str(exc_info.value)
    assert "seam=simple_kinase.interpreter.reference_coverage" in message
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
    request = SimpleKinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(top_k=3, ensemble_size=5),
        activity_config=None,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflow().run(request)

    message = str(exc_info.value)
    assert "seam=simple_kinase.interpreter.eligible_kinases" in message
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
    request = SimpleKinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(min_substrates=1),
        prediction_config=KinasePredictionConfig(top_k=2, ensemble_size=3),
        activity_config=None,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflow().run(request)

    message = str(exc_info.value)
    assert "seam=simple_kinase.executor.prediction_ensemble" in message
    assert "eligible_kinases=1" in message
    assert "ranked_kinases=0" in message
    assert "prediction_config_ensemble_size=3" in message
    assert "prediction_config_top_k=2" in message
    assert "dataset_samples=1" in message


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
    request = ResolvedKinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        kinase_substrate_map=references.kinase_substrate_map,
        scoring_config=KinaseScoringConfig(min_substrates=1),
        prediction_config=KinasePredictionConfig(top_k=3, ensemble_size=2),
        activity_config=KinaseActivityConfig(enabled=True, threshold=0.7),
    )
    prediction_result = KinasePredictionResult(
        pred_mat=pd.DataFrame(
            {"MAP2K6": [0.0]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SimpleKinaseWorkflowExecutor()._run_activity_stage(
            request=request,
            prediction_result=prediction_result,
        )

    message = str(exc_info.value)
    assert "seam=simple_kinase.executor.activity_support" in message
    assert "activity_kinases=1" in message
    assert "total_positive_predictions=0" in message
    assert "prediction_config_top_k=3" in message
    assert "activity_config_threshold=0.7" in message
