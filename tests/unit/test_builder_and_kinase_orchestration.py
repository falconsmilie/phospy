from __future__ import annotations

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DatasetBuildRequest,
    KinaseActivityConfig,
    KinaseActivityResult,
    KinasePredictionConfig,
    KinasePredictionResult,
    KinaseScoringConfig,
    KinaseScoringResult,
    Organism,
    ReferenceBundle,
    ReferencePreset,
    SimpleKinaseWorkflow,
    SimpleKinaseWorkflowRequest,
    SimpleKinaseWorkflowResult,
)
from phospy.datasets.builders.contracts import InterpretedDatasetBuildRequest
from phospy.workflows.kinase.contracts import ResolvedKinaseWorkflowRequest


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=["MAPK14;Y182;"],
    )


def _site_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
        },
        index=["MAPK14;Y182;"],
    )


def _dataset() -> AnalysisReadyPhosphoDataset:
    return AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
    )


def _bundle() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["MAP2K6"], "substrate_site": ["MAPK14;Y182;"]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
    )


def test_request_config_and_result_models_construct() -> None:
    dataset = _dataset()
    build_request = DatasetBuildRequest(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
    )
    workflow_request = SimpleKinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(min_substrates=1),
        prediction_config=KinasePredictionConfig(top_k=10, ensemble_size=3),
        activity_config=KinaseActivityConfig(enabled=True, threshold=0.5),
    )
    scoring = KinaseScoringResult(profile_scores=pd.DataFrame({"profile_score": []}))
    prediction = KinasePredictionResult(pred_mat=pd.DataFrame())
    activity = KinaseActivityResult(
        activity_scores=pd.DataFrame({"activity_score": []})
    )
    workflow_result = SimpleKinaseWorkflowResult(
        dataset=dataset,
        references=_bundle(),
        scoring_result=scoring,
        prediction_result=prediction,
        activity_result=activity,
    )
    assert isinstance(build_request, DatasetBuildRequest)
    assert isinstance(workflow_request, SimpleKinaseWorkflowRequest)
    assert isinstance(workflow_result, SimpleKinaseWorkflowResult)


def test_builder_run_contract_builds_analysis_ready_dataset() -> None:
    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            organism=Organism.RAT,
        )
    )
    assert isinstance(built, AnalysisReadyPhosphoDataset)
    assert isinstance(built.transformation_state.label, str)


def test_builder_orchestration_uses_collaborators() -> None:
    dataset = _dataset()
    request = DatasetBuildRequest(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
    )
    calls: list[str] = []
    interpreted = InterpretedDatasetBuildRequest(
        phospho=request.phospho,
        site_metadata=request.site_metadata,
        sample_metadata=request.sample_metadata,
        total=request.total,
        organism=request.organism,
        transformation_state=request.transformation_state,
    )

    class ValidatorSpy:
        def run(self, req: DatasetBuildRequest) -> DatasetBuildRequest:
            calls.append("validator")
            return req

    class InterpreterSpy:
        def run(self, req: DatasetBuildRequest) -> InterpretedDatasetBuildRequest:
            calls.append("interpreter")
            return interpreted

    class ExecutorSpy:
        def run(
            self, req: InterpretedDatasetBuildRequest
        ) -> AnalysisReadyPhosphoDataset:
            calls.append("executor")
            return dataset

    builder = AnalysisReadyDatasetBuilder(
        validator=ValidatorSpy(),
        interpreter=InterpreterSpy(),
        executor=ExecutorSpy(),
    )
    built = builder.run(request)
    assert built is dataset
    assert calls == ["validator", "interpreter", "executor"]


def test_workflow_run_contract_returns_nested_results() -> None:
    result = SimpleKinaseWorkflow().run(
        SimpleKinaseWorkflowRequest(dataset=_dataset(), references=ReferencePreset.AUTO)
    )
    assert isinstance(result, SimpleKinaseWorkflowResult)
    assert isinstance(result.scoring_result, KinaseScoringResult)
    assert isinstance(result.prediction_result, KinasePredictionResult)
    assert not hasattr(result, "profile_scores")
    assert not hasattr(result, "combined_scores")
    assert not hasattr(result, "weights")
    assert not hasattr(result, "substrate_list")


def test_workflow_public_entrypoint_exercises_collaborators() -> None:
    calls: list[str] = []
    request = SimpleKinaseWorkflowRequest(
        dataset=_dataset(), references=ReferencePreset.AUTO
    )
    interpreted = ResolvedKinaseWorkflowRequest(
        dataset=request.dataset,
        references=_bundle(),
        scoring_config=request.scoring_config,
        prediction_config=request.prediction_config,
        activity_config=request.activity_config,
    )
    expected = SimpleKinaseWorkflowResult(
        dataset=request.dataset,
        references=interpreted.references,
        scoring_result=KinaseScoringResult(profile_scores=pd.DataFrame()),
        prediction_result=KinasePredictionResult(pred_mat=pd.DataFrame()),
        activity_result=None,
    )

    class ValidatorSpy:
        def run(
            self, workflow_request: SimpleKinaseWorkflowRequest
        ) -> SimpleKinaseWorkflowRequest:
            calls.append("validator")
            return workflow_request

    class InterpreterSpy:
        def run(
            self, workflow_request: SimpleKinaseWorkflowRequest
        ) -> ResolvedKinaseWorkflowRequest:
            calls.append("interpreter")
            return interpreted

    class ExecutorSpy:
        def run(
            self, resolved: ResolvedKinaseWorkflowRequest
        ) -> SimpleKinaseWorkflowResult:
            calls.append("executor")
            return expected

    workflow = SimpleKinaseWorkflow(
        validator=ValidatorSpy(),
        interpreter=InterpreterSpy(),
        executor=ExecutorSpy(),
    )
    observed = workflow.run(request)
    assert observed is expected
    assert calls == ["validator", "interpreter", "executor"]
