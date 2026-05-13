from __future__ import annotations

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    KinaseWorkflow,
    SignalomeWorkflow,
)
from phospy.api import (
    DatasetBuildRequest,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    KinaseWorkflowResult,
    Organism,
    ReferenceBundle,
    ReferencePreset,
    SignalomeConfig,
    SignalomeWorkflowRequest,
    SignalomeWorkflowResult,
)
from phospy.api.results import (
    KinaseActivityResult,
    KinasePredictionResult,
    KinaseScoringResult,
)
from phospy.science.datasets.builders.contracts import InterpretedDatasetBuildRequest
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.science.prediction.policies import resolve_prediction_sampling_policy
from phospy.science.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
)
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseActivityExecutionConfig,
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeExecutionConfig,
    ResolvedSignalomeWorkflowRequest,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.signalome_config import build_signalome_config


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
            "localisation_confidence": [0.95, 0.9],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )


def _dataset() -> AnalysisReadyPhosphoDataset:
    return AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
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
            method=request.activity_config.method,
            threshold=float(request.activity_config.threshold),
            min_substrates=int(request.activity_config.min_substrates),
            top_n_substrates=int(request.activity_config.top_n_substrates),
            ksea_min_substrates=int(request.activity_config.ksea_min_substrates),
            ksea_evidence_threshold=float(
                request.activity_config.threshold
                if request.activity_config.ksea_evidence_threshold is None
                else request.activity_config.ksea_evidence_threshold
            ),
            ksea_p_value_method=request.activity_config.ksea_p_value_method,
            ksea_adjust_p_values=bool(request.activity_config.ksea_adjust_p_values),
        )
    )
    return ResolvedKinaseExecutionConfig(
        scoring_min_substrates=int(request.scoring_config.min_substrates),
        include_diagnostic_scoring_tables=bool(
            request.scoring_config.include_diagnostic_scoring_tables
        ),
        profile_missing_value_strategy=request.scoring_config.profile_missing_value_strategy,
        prediction_top_k=int(request.prediction_config.top_k),
        prediction_deterministic_max_selected_kinases=int(
            request.prediction_config.deterministic_max_selected_kinases
        ),
        prediction_adaptive_ensemble_runs=int(
            request.prediction_config.adaptive_ensemble_runs
        ),
        prediction_mode=request.prediction_config.mode,
        prediction_adaptive_policy=request.prediction_config.adaptive_policy,
        prediction_n_iterations=int(request.prediction_config.n_iterations),
        prediction_random_state=request.prediction_config.random_state,
        prediction_sampling_policy=resolve_prediction_sampling_policy(
            request.prediction_config.adaptive_policy
        ),
        activity=activity,
    )


def _resolved_signalome_execution_config(
    config: SignalomeConfig,
) -> ResolvedSignalomeExecutionConfig:
    return ResolvedSignalomeExecutionConfig(
        substrate_support_cutoff=float(config.scientific.substrate_support_cutoff),
        network_correlation_threshold=float(
            config.output.network_correlation_threshold
        ),
        network_policy=config.output.network_policy,
        assignment_policy=config.scientific.assignment_policy,
        score_preconditioning_policy=config.validation.score_preconditioning_policy,
        module_selection_primary_threshold=float(
            config.clustering.module_selection_primary_correlation_threshold
        ),
        module_selection_fallback_threshold=float(
            config.clustering.module_selection_fallback_correlation_threshold
        ),
        module_selection_max_clusters=int(
            config.clustering.module_selection_max_clusters
        ),
        candidate_scoring_policy=config.clustering.candidate_scoring_policy,
        max_exact_tree_sites=int(config.performance.max_exact_tree_sites),
        max_full_candidate_scoring_sites=int(
            config.performance.max_full_candidate_scoring_sites
        ),
        requested_module_count=(
            None
            if config.clustering.module_count is None
            else int(config.clustering.module_count)
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
        prediction_config=KinasePredictionConfig(
            top_k=10,
            deterministic_max_selected_kinases=3,
            adaptive_ensemble_runs=3,
        ),
        activity_config=KinaseActivityConfig(
            enabled=True,
            threshold=0.5,
            min_substrates=2,
            top_n_substrates=5,
        ),
    )
    scoring = KinaseScoringResult(
        profile_scores=pd.DataFrame(
            {"MAP2K6": [0.75]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        )
    )
    prediction = KinasePredictionResult(
        pred_mat=pd.DataFrame(
            {"MAP2K6": [0.75]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        )
    )
    activity = KinaseActivityResult(
        weighted_activity=pd.DataFrame(),
        thresholded_substrate_mean_activity=pd.DataFrame(),
        thresholded_substrate_counts=pd.Series(dtype="int64", name="n_substrates"),
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
    assert isinstance(built.intensity_scale_state.label, str)


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


def test_builder_orchestration_threads_preprocessing_config_to_executor() -> None:
    dataset = _dataset()
    request = DatasetBuildRequest(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        preprocessing_config=DatasetPreprocessingConfig(
            missing_data=DatasetMissingDataConfig(
                policy="impute_row_median",
                min_observed_values=2,
            ),
        ),
    )
    calls: list[str] = []

    class ValidatorSpy:
        def run(self, req: DatasetBuildRequest) -> DatasetBuildRequest:
            calls.append("validator")
            return req

    class InterpreterSpy:
        def run(self, req: DatasetBuildRequest) -> InterpretedDatasetBuildRequest:
            calls.append("interpreter")
            return InterpretedDatasetBuildRequest(
                phospho=req.phospho,
                site_metadata=req.site_metadata,
                sample_metadata=req.sample_metadata,
                total=req.total,
                organism=req.organism,
                preprocessing_plan=PreprocessingPlan.from_config(
                    req.preprocessing_config
                ),
            )

    class ExecutorSpy:
        def run(
            self, req: InterpretedDatasetBuildRequest
        ) -> AnalysisReadyPhosphoDataset:
            calls.append("executor")
            assert req.preprocessing_plan == PreprocessingPlan.from_config(
                request.preprocessing_config
            )
            return dataset

    built = AnalysisReadyDatasetBuilder(
        validator=ValidatorSpy(),
        interpreter=InterpreterSpy(),
        executor=ExecutorSpy(),
    ).run(request)
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
    assert not hasattr(result, "rank_weighted_fusion_scores")
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
        execution_config=_resolved_kinase_execution_config(request),
    )
    expected = KinaseWorkflowResult(
        dataset=request.dataset,
        references=interpreted.references,
        scoring_result=KinaseScoringResult(
            profile_scores=pd.DataFrame(
                {"MAP2K6": [0.75]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            )
        ),
        prediction_result=KinasePredictionResult(
            pred_mat=pd.DataFrame(
                {"MAP2K6": [0.75]},
                index=pd.Index(["MAPK14;Y182;"], name="site_id"),
            )
        ),
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
        config=build_signalome_config(substrate_support_cutoff=0.5),
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
            ["MAPK14", "GSK3B"],
            index=score_matrix.index.copy(),
            name="protein_id",
            dtype=str,
        ),
    )
    expected = SignalomeWorkflowResult(
        dataset=kinase_result.dataset,
        kinase_result=kinase_result,
        module_assignments=SignalomeAssignments(
            table=pd.DataFrame(
                {
                    "protein_id": ["MAPK14"],
                    "module_id": [1],
                    "top_kinase": ["MAP2K6"],
                    "top_score": [0.9],
                    "top_kinase_candidates": [("MAP2K6",)],
                    "top_kinase_weights": [(("MAP2K6", 1.0),)],
                    "top_kinase_tie_count": [1],
                    "top_kinase_is_ambiguous": [False],
                    "top_kinase_selection_policy": [
                        "max_score_then_lexicographic_tiebreak"
                    ],
                    "module_top_kinase": ["MAP2K6"],
                    "module_top_kinase_candidates": [("MAP2K6",)],
                    "module_top_kinase_tie_count": [1],
                    "module_top_kinase_is_ambiguous": [False],
                    "module_top_kinase_selection_policy": [
                        "max_score_then_lexicographic_tiebreak"
                    ],
                },
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
