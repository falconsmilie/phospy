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
    KinaseWorkflow,
    KinaseWorkflowRequest,
    KinaseWorkflowResult,
    Organism,
    ReferenceBundle,
    ReferencePreset,
    SignalomeConfig,
    SignalomeWorkflow,
    SignalomeWorkflowRequest,
    SignalomeWorkflowResult,
)
from phospy.datasets.builders.contracts import InterpretedDatasetBuildRequest
from phospy.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
)
from phospy.transformations.models import TransformationState
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseActivityExecutionConfig,
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeExecutionConfig,
    ResolvedSignalomeWorkflowRequest,
)


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0, 0.8], "sample_b": [2.0, 1.2]},
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )


def _site_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
            ],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )


def _dataset() -> AnalysisReadyPhosphoDataset:
    return AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        transformation_state=TransformationState.raw(has_total_matrix=False),
    )


def _bundle() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6"],
                "substrate_site": ["MAPK14;Y182;", "GSK3B;S9;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                    "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
                ]
            },
            index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
        ),
    )


def _resolved_kinase_execution_config(
    request: KinaseWorkflowRequest,
) -> ResolvedKinaseExecutionConfig:
    activity = (
        None
        if request.activity_config is None or not request.activity_config.enabled
        else ResolvedKinaseActivityExecutionConfig(
            threshold=float(request.activity_config.threshold),
            min_substrates=int(request.activity_config.min_substrates),
            top_n_substrates=int(request.activity_config.top_n_substrates),
        )
    )
    return ResolvedKinaseExecutionConfig(
        scoring_min_substrates=int(request.scoring_config.min_substrates),
        include_diagnostic_scoring_tables=bool(
            request.scoring_config.include_diagnostic_scoring_tables
        ),
        profile_missing_value_strategy=request.scoring_config.profile_missing_value_strategy,
        prediction_top_k=int(request.prediction_config.top_k),
        prediction_ensemble_size=int(request.prediction_config.ensemble_size),
        prediction_mode=request.prediction_config.mode,
        prediction_adaptive_policy=request.prediction_config.adaptive_policy,
        prediction_n_iterations=int(request.prediction_config.n_iterations),
        prediction_random_state=request.prediction_config.random_state,
        activity=activity,
    )


def _resolved_signalome_execution_config(
    config: SignalomeConfig,
) -> ResolvedSignalomeExecutionConfig:
    return ResolvedSignalomeExecutionConfig(
        substrate_support_cutoff=float(config.substrate_support_cutoff),
        network_correlation_threshold=float(config.network_correlation_threshold),
        network_policy=config.network_policy,
        assignment_policy=config.assignment_policy,
        module_selection_primary_threshold=float(
            config.module_selection_primary_correlation_threshold
        ),
        module_selection_fallback_threshold=float(
            config.module_selection_fallback_correlation_threshold
        ),
        module_selection_max_clusters=int(config.module_selection_max_clusters),
        requested_module_count=(
            None if config.module_count is None else int(config.module_count)
        ),
    )


def test_request_config_and_result_models_construct() -> None:
    dataset = _dataset()
    build_request = DatasetBuildRequest(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
    )
    workflow_request = KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(top_k=10, ensemble_size=3),
        activity_config=KinaseActivityConfig(
            enabled=True,
            threshold=0.5,
            min_substrates=2,
            top_n_substrates=5,
        ),
    )
    scoring = KinaseScoringResult(profile_scores=pd.DataFrame({"profile_score": []}))
    prediction = KinasePredictionResult(pred_mat=pd.DataFrame())
    activity = KinaseActivityResult(
        weighted_activity=pd.DataFrame(),
        ksea_scores=pd.DataFrame(),
        ksea_counts=pd.Series(dtype="int64", name="n_substrates"),
        target_counts=pd.Series(dtype="int64", name="n_targets"),
        target_table=pd.DataFrame(columns=["site_id", "kinase", "score"]),
    )
    workflow_result = KinaseWorkflowResult(
        dataset=dataset,
        references=_bundle(),
        scoring_result=scoring,
        prediction_result=prediction,
        activity_result=activity,
    )
    assert isinstance(build_request, DatasetBuildRequest)
    assert isinstance(workflow_request, KinaseWorkflowRequest)
    assert isinstance(workflow_result, KinaseWorkflowResult)


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
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_bundle(),
            scoring_config=KinaseScoringConfig(min_substrates=2),
            activity_config=None,
        )
    )
    assert isinstance(result, KinaseWorkflowResult)
    assert isinstance(result.scoring_result, KinaseScoringResult)
    assert isinstance(result.prediction_result, KinasePredictionResult)
    assert not hasattr(result, "profile_scores")
    assert not hasattr(result, "combined_scores")
    assert not hasattr(result, "weights")
    assert not hasattr(result, "substrate_list")


def test_workflow_public_entrypoint_exercises_collaborators() -> None:
    calls: list[str] = []
    request = KinaseWorkflowRequest(dataset=_dataset(), references=ReferencePreset.AUTO)
    bundle = _bundle()
    interpreted = ResolvedKinaseWorkflowRequest(
        dataset=request.dataset,
        references=bundle,
        kinase_substrate_map=bundle.kinase_substrate_map,
        site_sequences=bundle.site_sequences,
        scoring_site_index=request.dataset.phospho.index.copy(),
        activity_phospho_matrix=request.dataset.phospho.copy(deep=True),
        uses_bundled_reference=True,
        execution_config=_resolved_kinase_execution_config(request),
    )
    expected = KinaseWorkflowResult(
        dataset=request.dataset,
        references=interpreted.references,
        scoring_result=KinaseScoringResult(profile_scores=pd.DataFrame()),
        prediction_result=KinasePredictionResult(pred_mat=pd.DataFrame()),
        activity_result=None,
    )

    class ValidatorSpy:
        def run(self, workflow_request: KinaseWorkflowRequest) -> KinaseWorkflowRequest:
            calls.append("validator")
            return workflow_request

    class InterpreterSpy:
        def run(
            self, workflow_request: KinaseWorkflowRequest
        ) -> ResolvedKinaseWorkflowRequest:
            calls.append("interpreter")
            return interpreted

    class ExecutorSpy:
        def run(self, resolved: ResolvedKinaseWorkflowRequest) -> KinaseWorkflowResult:
            calls.append("executor")
            return expected

    workflow = KinaseWorkflow(
        validator=ValidatorSpy(),
        interpreter=InterpreterSpy(),
        executor=ExecutorSpy(),
    )
    observed = workflow.run(request)
    assert observed is expected
    assert calls == ["validator", "interpreter", "executor"]


def test_signalome_workflow_public_entrypoint_exercises_collaborators() -> None:
    calls: list[str] = []
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_bundle(),
            scoring_config=KinaseScoringConfig(min_substrates=2),
            activity_config=None,
        )
    )
    request = SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )
    score_matrix = kinase_result.scoring_result.profile_scores
    interpreted = ResolvedSignalomeWorkflowRequest(
        dataset=kinase_result.dataset,
        kinase_result=kinase_result,
        execution_config=_resolved_signalome_execution_config(request.config),
        downstream_score_matrix=score_matrix,
        downstream_score_source="profile_scores",
        prediction_matrix=kinase_result.prediction_result.pred_mat,
        site_to_protein=pd.Series(
            ["MAPK14"],
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            name="protein_id",
            dtype=str,
        ),
    )
    expected = SignalomeWorkflowResult(
        dataset=kinase_result.dataset,
        kinase_result=kinase_result,
        module_assignments=SignalomeAssignments(
            table=pd.DataFrame(
                {"protein_id": ["MAPK14"], "module_id": [1]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            )
        ),
        signalome_modules=SignalomeModules(
            table=pd.DataFrame(
                {"MAP2K6": [100.0]},
                index=pd.Index([1], name="module_id"),
            )
        ),
        kinase_network=KinaseNetwork(
            edges=pd.DataFrame(
                columns=["source_kinase", "target_kinase", "correlation"]
            ),
            nodes=pd.DataFrame(
                {"degree": [0], "n_substrates": [1]},
                index=pd.Index(["MAP2K6"], name="kinase"),
            ),
        ),
        expanded_signalome=None,
    )

    class ValidatorSpy:
        def run(
            self, workflow_request: SignalomeWorkflowRequest
        ) -> SignalomeWorkflowRequest:
            calls.append("validator")
            return workflow_request

    class InterpreterSpy:
        def run(
            self, workflow_request: SignalomeWorkflowRequest
        ) -> ResolvedSignalomeWorkflowRequest:
            calls.append("interpreter")
            return interpreted

    class ExecutorSpy:
        def run(
            self, resolved: ResolvedSignalomeWorkflowRequest
        ) -> SignalomeWorkflowResult:
            calls.append("executor")
            return expected

    workflow = SignalomeWorkflow(
        validator=ValidatorSpy(),
        interpreter=InterpreterSpy(),
        executor=ExecutorSpy(),
    )
    observed = workflow.run(request)
    assert observed is expected
    assert calls == ["validator", "interpreter", "executor"]


def test_signalome_workflow_exposes_only_run_entrypoint() -> None:
    workflow = SignalomeWorkflow()
    assert callable(getattr(workflow, "run", None))
    assert not hasattr(workflow, "run_from_analysis_ready")
