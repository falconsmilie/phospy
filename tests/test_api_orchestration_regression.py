from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass
from types import SimpleNamespace

import pandas as pd
import pytest

import phospy
import phospy.api.signalome_workflows as signalome_workflows_module
import phospy.api.simple_workflows as simple_workflows_module
from phospy.api import (
    DatasetLoadOptions,
    KinaseActivityConfig,
    PredictionRunConfig,
    SignalomeRunConfig,
    SignalomeWorkflow,
    SimpleKinaseWorkflow,
)
from phospy.api.simple_workflow_composition import SimpleKinaseExecutionGraph
from phospy.api.workflow_results import SimpleKinaseWorkflowResult
from phospy.prediction import PredMatResult
from phospy.preprocessing import CorePreprocessingConfig
from phospy.references import ReferenceBundle


@dataclass
class _BuilderDouble:
    dataset: object
    calls: list[dict[str, object]]

    def build(self, **kwargs: object) -> object:
        self.calls.append(kwargs)
        return self.dataset


@dataclass
class _ProviderDouble:
    bundle: object
    calls: list[dict[str, object]]

    def resolve(self, **kwargs: object) -> object:
        self.calls.append(kwargs)
        return self.bundle


@dataclass
class _AnalyzerDouble:
    result: object
    calls: list[dict[str, object]]

    def run(self, **kwargs: object) -> object:
        self.calls.append(kwargs)
        return self.result


@dataclass
class _WorkflowExecutorDouble:
    validate_result: object
    execute_result: object
    validate_calls: list[dict[str, object]]
    execute_calls: list[object]

    def validate_request(self, **kwargs: object) -> object:
        self.validate_calls.append(kwargs)
        return self.validate_result

    def execute_validated_request(self, request: object) -> object:
        self.execute_calls.append(request)
        return self.execute_result


class _SignalomeRequestDouble:
    pass


class _SignalomeResultDouble:
    pass


class _ReferenceBundleDouble(ReferenceBundle):
    pass


@dataclass
class _PredictionResultDouble:
    pred_mat_result: object
    substrate_list: dict[str, list[str]] = field(default_factory=dict)
    close_calls: int = 0

    def close(self) -> None:
        self.close_calls += 1


def test_simple_workflow_result_contract_shape_is_stable() -> None:
    assert is_dataclass(SimpleKinaseWorkflowResult)
    assert [item.name for item in fields(SimpleKinaseWorkflowResult)] == [
        "analysis_ready_dataset",
        "reference_bundle",
        "scoring_result",
        "prediction_result",
        "kinase_activity_result",
    ]
    assert tuple(SimpleKinaseWorkflowResult.__slots__) == (
        "analysis_ready_dataset",
        "reference_bundle",
        "scoring_result",
        "prediction_result",
        "kinase_activity_result",
    )
    assert "pred_mat_result" not in SimpleKinaseWorkflowResult.__dataclass_fields__
    assert "profile_scores" not in SimpleKinaseWorkflowResult.__dataclass_fields__
    assert "combined_scores" not in SimpleKinaseWorkflowResult.__dataclass_fields__
    assert "weights" not in SimpleKinaseWorkflowResult.__dataclass_fields__
    assert "substrate_list" not in SimpleKinaseWorkflowResult.__dataclass_fields__


def make_small_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [3.0, 4.0]},
        index=["SITE_1", "SITE_2"],
    )


def test_simple_kinase_workflow_run_delegates_to_domain_services() -> None:
    phospho = pd.DataFrame({"uid": ["u1"], "gene_names": ["PRKACA"]})
    total = pd.DataFrame({"genes": ["PRKACA"]})
    phospho_matrix = make_small_matrix()
    site_sequences = pd.Series(
        {"SITE_1": "AAAA", "SITE_2": "BBBB"},
        name="site_sequence",
    )
    analysis_ready_dataset = SimpleNamespace(
        phospho_matrix=phospho_matrix,
        site_sequences=site_sequences,
    )
    builder_calls: list[dict[str, object]] = []
    provider_calls: list[dict[str, object]] = []
    validate_calls: list[dict[str, object]] = []
    execute_calls: list[object] = []
    analyzer_calls: list[dict[str, object]] = []
    reference_bundle = _ReferenceBundleDouble.__new__(_ReferenceBundleDouble)
    pred_mat_result = object()
    scoring_result = object()
    prediction_result = _PredictionResultDouble(pred_mat_result=pred_mat_result)
    validated_request = object()
    workflow_result = SimpleNamespace(
        scoring_result=scoring_result,
        prediction_result=prediction_result,
        profile_result=object(),
        motif_result=object(),
    )
    kinase_activity_result = object()
    workflow = SimpleKinaseWorkflow(
        analysis_ready_builder=_BuilderDouble(
            dataset=analysis_ready_dataset,
            calls=builder_calls,
        ),
        reference_provider=_ProviderDouble(
            bundle=reference_bundle,
            calls=provider_calls,
        ),
        workflow_executor=_WorkflowExecutorDouble(
            validate_result=validated_request,
            execute_result=workflow_result,
            validate_calls=validate_calls,
            execute_calls=execute_calls,
        ),
        activity_analyzer=_AnalyzerDouble(
            result=kinase_activity_result,
            calls=analyzer_calls,
        ),
    )

    result = workflow.run(
        phospho=phospho,
        total=total,
        species="human",
        reference="ochoa",
        dataset_options=DatasetLoadOptions(phospho_encoding="utf-8"),
        preprocessing_config=CorePreprocessingConfig(
            min_observed=5,
            max_unmatched_fraction=0.25,
        ),
        prediction_config=PredictionRunConfig(
            min_substrates=3,
            score_threshold=0.65,
        ),
        activity_config=KinaseActivityConfig(
            threshold=0.55,
            min_substrates=4,
            top_n_substrates=12,
        ),
    )

    assert isinstance(result, SimpleKinaseWorkflowResult)
    assert result.analysis_ready_dataset is analysis_ready_dataset
    assert result.reference_bundle is reference_bundle
    assert result.scoring_result is scoring_result
    assert result.prediction_result is prediction_result
    assert result.pred_mat_result is pred_mat_result
    assert result.kinase_activity_result is kinase_activity_result
    assert builder_calls == [
        {
            "phospho": phospho,
            "total": total,
            "phospho_encoding": "utf-8",
            "schema": builder_calls[0]["schema"],
            "comparisons": None,
            "preprocessing_config": CorePreprocessingConfig(
                min_observed=5,
                max_unmatched_fraction=0.25,
            ),
            "source": "simple kinase workflow",
            "phospho_only_source": "simple kinase workflow (phospho only)",
        }
    ]
    assert provider_calls == [{"species": "human", "reference": "ochoa"}]
    assert len(validate_calls) == 1
    assert validate_calls[0]["phospho_matrix"] is phospho_matrix
    assert validate_calls[0]["site_sequences"] is site_sequences
    assert validate_calls[0]["reference_bundle"] is reference_bundle
    assert validate_calls[0]["min_substrates"] == 3
    assert validate_calls[0]["score_threshold"] == 0.65
    assert execute_calls == [validated_request]
    assert analyzer_calls == [
        {
            "pred_mat": pred_mat_result,
            "phospho_matrix": phospho_matrix,
            "threshold": 0.55,
            "min_substrates": 4,
            "top_n_substrates": 12,
        }
    ]


def test_simple_workflow_result_pred_mat_result_is_delegated_only() -> None:
    first_pred_mat_result = PredMatResult(
        pd.DataFrame({"KINASE_A": [0.1]}, index=["SITE_1"])
    )
    second_pred_mat_result = PredMatResult(
        pd.DataFrame({"KINASE_A": [0.9]}, index=["SITE_1"])
    )
    profile_scores = pd.DataFrame({"KINASE_A": [0.1]}, index=["SITE_1"])
    combined_scores = pd.DataFrame({"KINASE_A": [0.2]}, index=["SITE_1"])
    weights = pd.DataFrame(
        {"motif_weight": [0.4], "profile_weight": [0.6]},
        index=["KINASE_A"],
    )
    substrate_list = {"KINASE_A": ["SITE_1"]}

    prediction_result = _PredictionResultDouble(
        pred_mat_result=first_pred_mat_result,
        substrate_list=substrate_list,
    )
    result = SimpleKinaseWorkflowResult(
        analysis_ready_dataset=SimpleNamespace(),
        reference_bundle=SimpleNamespace(),
        scoring_result=SimpleNamespace(
            profile_scores=profile_scores,
            combined_scores=combined_scores,
            weights=weights,
        ),
        prediction_result=prediction_result,
        kinase_activity_result=SimpleNamespace(),
    )

    assert isinstance(result.pred_mat_result, PredMatResult)
    assert result.pred_mat_result is first_pred_mat_result
    assert result.profile_scores is profile_scores
    assert result.combined_scores is combined_scores
    assert result.weights is weights
    assert result.substrate_list is substrate_list

    result.prediction_result.pred_mat_result = second_pred_mat_result

    assert result.pred_mat_result is second_pred_mat_result


def test_simple_workflow_result_close_delegates_to_prediction_result_close() -> None:
    prediction_result = _PredictionResultDouble(pred_mat_result=object())
    result = SimpleKinaseWorkflowResult(
        analysis_ready_dataset=SimpleNamespace(),
        reference_bundle=SimpleNamespace(),
        scoring_result=SimpleNamespace(
            profile_scores=pd.DataFrame(),
            combined_scores=None,
            weights=None,
        ),
        prediction_result=prediction_result,
        kinase_activity_result=SimpleNamespace(),
    )

    result.close()

    assert prediction_result.close_calls == 1


def test_simple_workflow_result_context_manager_closes_on_exit() -> None:
    prediction_result = _PredictionResultDouble(pred_mat_result=object())
    result = SimpleKinaseWorkflowResult(
        analysis_ready_dataset=SimpleNamespace(),
        reference_bundle=SimpleNamespace(),
        scoring_result=SimpleNamespace(
            profile_scores=pd.DataFrame(),
            combined_scores=None,
            weights=None,
        ),
        prediction_result=prediction_result,
        kinase_activity_result=SimpleNamespace(),
    )

    with result as managed:
        assert managed is result
        assert prediction_result.close_calls == 0

    assert prediction_result.close_calls == 1


def test_simple_kinase_workflow_run_delegates_to_execution_service() -> None:
    expected_result = object()
    calls: list[dict[str, object]] = []

    class _ExecutionServiceDouble:
        def run(self, **kwargs: object) -> object:
            calls.append(kwargs)
            return expected_result

    workflow = SimpleKinaseWorkflow(execution_service=_ExecutionServiceDouble())
    result = workflow.run(
        phospho=pd.DataFrame({"uid": ["u1"]}),
        species="rat",
    )

    assert result is expected_result
    assert len(calls) == 1
    assert calls[0]["species"] == "rat"
    assert calls[0]["reference"] == "auto"


def test_simple_kinase_workflow_uses_execution_graph_override() -> None:
    phospho = pd.DataFrame({"uid": ["u1"], "gene_names": ["PRKACA"]})
    phospho_matrix = make_small_matrix()
    site_sequences = pd.Series(
        {"SITE_1": "AAAA", "SITE_2": "BBBB"},
        name="site_sequence",
    )
    analysis_ready_dataset = SimpleNamespace(
        phospho_matrix=phospho_matrix,
        site_sequences=site_sequences,
    )
    builder_calls: list[dict[str, object]] = []
    provider_calls: list[dict[str, object]] = []
    validate_calls: list[dict[str, object]] = []
    execute_calls: list[object] = []
    analyzer_calls: list[dict[str, object]] = []
    reference_bundle = _ReferenceBundleDouble.__new__(_ReferenceBundleDouble)
    pred_mat_result = object()
    prediction_result = _PredictionResultDouble(pred_mat_result=pred_mat_result)
    validated_request = object()
    workflow_result = SimpleNamespace(
        scoring_result=object(),
        prediction_result=prediction_result,
        profile_result=object(),
        motif_result=object(),
    )
    workflow = SimpleKinaseWorkflow(
        execution_graph=SimpleKinaseExecutionGraph(
            analysis_ready_builder=_BuilderDouble(
                dataset=analysis_ready_dataset,
                calls=builder_calls,
            ),
            reference_provider=_ProviderDouble(
                bundle=reference_bundle,
                calls=provider_calls,
            ),
            workflow_executor=_WorkflowExecutorDouble(
                validate_result=validated_request,
                execute_result=workflow_result,
                validate_calls=validate_calls,
                execute_calls=execute_calls,
            ),
            activity_analyzer=_AnalyzerDouble(
                result=object(),
                calls=analyzer_calls,
            ),
        )
    )

    workflow.run(phospho=phospho, species="human")

    assert len(builder_calls) == 1
    assert provider_calls == [{"species": "human", "reference": "auto"}]
    assert validate_calls[0]["phospho_matrix"] is phospho_matrix
    assert execute_calls == [validated_request]
    assert analyzer_calls[0]["pred_mat"] is pred_mat_result


def test_simple_kinase_workflow_builds_default_graph_from_composition_layer(
    monkeypatch,
) -> None:
    graph_calls: list[dict[str, object]] = []
    fake_graph = SimpleKinaseExecutionGraph(
        analysis_ready_builder=_BuilderDouble(
            dataset=SimpleNamespace(
                phospho_matrix=make_small_matrix(),
                site_sequences=pd.Series({"SITE_1": "AAAA"}),
            ),
            calls=[],
        ),
        reference_provider=_ProviderDouble(
            bundle=_ReferenceBundleDouble.__new__(_ReferenceBundleDouble),
            calls=[],
        ),
        workflow_executor=_WorkflowExecutorDouble(
            validate_result=object(),
            execute_result=SimpleNamespace(
                scoring_result=object(),
                prediction_result=_PredictionResultDouble(pred_mat_result=object()),
                profile_result=object(),
                motif_result=object(),
            ),
            validate_calls=[],
            execute_calls=[],
        ),
        activity_analyzer=_AnalyzerDouble(result=object(), calls=[]),
    )

    def fake_create_default_simple_kinase_execution_graph(
        **kwargs: object,
    ) -> SimpleKinaseExecutionGraph:
        graph_calls.append(kwargs)
        return fake_graph

    monkeypatch.setattr(
        simple_workflows_module,
        "create_default_simple_kinase_execution_graph",
        fake_create_default_simple_kinase_execution_graph,
    )

    workflow = SimpleKinaseWorkflow(flank_size=9, kernel="linear")

    assert workflow._analysis_ready_builder is fake_graph.analysis_ready_builder
    assert workflow._reference_provider is fake_graph.reference_provider
    assert workflow._activity_analyzer is fake_graph.activity_analyzer
    assert workflow._workflow_executor is fake_graph.workflow_executor
    assert len(graph_calls) == 1
    assert graph_calls[0]["flank_size"] == 9
    assert graph_calls[0]["kernel"] == "linear"


def test_simple_kinase_workflow_rejects_conflicting_execution_inputs() -> None:
    with pytest.raises(ValueError, match="execution_service"):
        SimpleKinaseWorkflow(
            execution_service=SimpleNamespace(run=lambda **_: object()),
            execution_graph=SimpleKinaseExecutionGraph(
                analysis_ready_builder=_BuilderDouble(
                    dataset=SimpleNamespace(),
                    calls=[],
                ),
                reference_provider=_ProviderDouble(
                    bundle=SimpleNamespace(),
                    calls=[],
                ),
                workflow_executor=_WorkflowExecutorDouble(
                    validate_result=SimpleNamespace(),
                    execute_result=SimpleNamespace(),
                    validate_calls=[],
                    execute_calls=[],
                ),
                activity_analyzer=_AnalyzerDouble(
                    result=SimpleNamespace(),
                    calls=[],
                ),
            ),
        )


def test_signalome_workflow_run_delegates_to_validation_and_execution(
    monkeypatch,
) -> None:
    workflow = SignalomeWorkflow()
    request = _SignalomeRequestDouble()
    signalome_result = _SignalomeResultDouble()
    calls: list[tuple[str, object]] = []
    scoring_result = object()
    prediction_result = object()
    expression_matrix = make_small_matrix()

    def fake_validate_signalome_request(**kwargs: object) -> object:
        calls.append(("validate", kwargs))
        return request

    def fake_execute_signalome_inputs(received_request: object) -> object:
        calls.append(("execute", received_request))
        return signalome_result

    monkeypatch.setattr(
        signalome_workflows_module,
        "validate_signalome_request",
        fake_validate_signalome_request,
    )
    monkeypatch.setattr(
        signalome_workflows_module,
        "execute_signalome_inputs",
        fake_execute_signalome_inputs,
    )

    result = workflow.run(
        scoring_result=scoring_result,
        prediction_result=prediction_result,
        expression_matrix=expression_matrix,
        kinases_of_interest=["PRKACA", "BTK"],
        config=SignalomeRunConfig(
            kinase_network_threshold=0.8,
            kinase_network_policy="signed",
            assignment_policy="weighted_top",
            signalome_cutoff=0.4,
        ),
    )

    assert result is signalome_result
    assert calls[0][0] == "validate"
    assert calls[0][1]["scoring_result"] is scoring_result
    assert calls[0][1]["prediction_result"] is prediction_result
    assert calls[0][1]["expression_matrix"] is expression_matrix
    assert calls[0][1]["kinases_of_interest"] == ["PRKACA", "BTK"]
    assert calls[0][1]["kinase_network_threshold"] == 0.8
    assert calls[0][1]["kinase_network_policy"] == "signed"
    assert calls[0][1]["assignment_policy"] == "weighted_top"
    assert calls[0][1]["signalome_cutoff"] == 0.4
    assert calls[1] == ("execute", request)


def test_root_package_does_not_reexport_workflow_aliases() -> None:
    assert not hasattr(phospy, "KinaseWorkflow")
    assert not hasattr(phospy, "PredMatWorkflow")
    assert not hasattr(phospy, "SignalomeWorkflow")
    assert not hasattr(phospy, "SimpleKinaseWorkflow")
