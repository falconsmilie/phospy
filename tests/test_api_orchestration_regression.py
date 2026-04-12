from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import pandas as pd

import phospy
import phospy.api.workflows as api_workflows_module
from phospy.api import (
    KinaseWorkflow,
    PredMatWorkflow,
    SignalomeWorkflow,
    SimpleKinaseWorkflow,
)
from phospy.api.workflows import (
    KinaseWorkflowResult,
    PredMatWorkflowResult,
    SimpleKinaseWorkflowResult,
)
from phospy.internal.constants import DEFAULT_PHOSPHO_SENTINEL, DEFAULT_TOTAL_SENTINEL
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
class _PredMatWorkflowDouble:
    result: object
    calls: list[dict[str, object]]

    def run(self, **kwargs: object) -> object:
        self.calls.append(kwargs)
        return self.result


class _ClosablePredictionResult:
    def __init__(self, pred_mat_result: object) -> None:
        self.pred_mat_result = pred_mat_result
        self.closed = False

    def close(self) -> None:
        self.closed = True


class _SignalomeRequestDouble:
    pass


class _SignalomeResultDouble:
    pass


class _ReferenceBundleDouble(ReferenceBundle):
    pass


def make_small_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0, 2.0], "sample_b": [3.0, 4.0]},
        index=["SITE_1", "SITE_2"],
    )


def test_kinase_workflow_run_delegates_to_prediction_executor(
    monkeypatch,
) -> None:
    workflow = KinaseWorkflow()
    phospho_matrix = make_small_matrix()
    request = object()
    profile_result = object()
    motif_result = object()
    scoring_result = object()
    prediction_result = object()
    calls: list[tuple[str, object]] = []

    def fake_validate_request(**kwargs: object) -> object:
        calls.append(("validate", kwargs))
        return request

    def fake_execute_validated_request(received_request: object) -> object:
        calls.append(("execute", received_request))
        return SimpleNamespace(
            profile_result=profile_result,
            motif_result=motif_result,
            scoring_result=scoring_result,
            prediction_result=prediction_result,
        )

    monkeypatch.setattr(workflow._executor, "validate_request", fake_validate_request)
    monkeypatch.setattr(
        workflow._executor,
        "execute_validated_request",
        fake_execute_validated_request,
    )

    result = workflow.run(
        phospho_matrix=phospho_matrix,
        site_sequences={"SITE_1": "AAAA", "SITE_2": "BBBB"},
        min_substrates=2,
        ensemble_size=3,
        random_state=17,
    )

    assert isinstance(result, KinaseWorkflowResult)
    assert result.profile_result is profile_result
    assert result.motif_result is motif_result
    assert result.scoring_result is scoring_result
    assert result.prediction_result is prediction_result
    assert calls[0][0] == "validate"
    assert calls[0][1]["phospho_matrix"] is phospho_matrix
    assert calls[0][1]["min_substrates"] == 2
    assert calls[0][1]["ensemble_size"] == 3
    assert calls[0][1]["random_state"] == 17
    assert calls[1] == ("execute", request)


def test_pred_mat_workflow_run_delegates_to_prediction_executor(
    monkeypatch,
) -> None:
    workflow = PredMatWorkflow()
    phospho_matrix = make_small_matrix()
    request = object()
    scoring_result = object()
    pred_mat_result = object()
    prediction_result = _ClosablePredictionResult(pred_mat_result=pred_mat_result)
    calls: list[tuple[str, object]] = []

    def fake_validate_request(**kwargs: object) -> object:
        calls.append(("validate", kwargs))
        return request

    def fake_execute_validated_request(received_request: object) -> object:
        calls.append(("execute", received_request))
        return SimpleNamespace(
            scoring_result=scoring_result,
            prediction_result=prediction_result,
        )

    monkeypatch.setattr(workflow._executor, "validate_request", fake_validate_request)
    monkeypatch.setattr(
        workflow._executor,
        "execute_validated_request",
        fake_execute_validated_request,
    )

    result = workflow.run(
        phospho_matrix=phospho_matrix,
        site_sequences={"SITE_1": "AAAA", "SITE_2": "BBBB"},
        score_threshold=0.7,
        inclusion=11,
    )

    assert isinstance(result, PredMatWorkflowResult)
    assert result.scoring_result is scoring_result
    assert result.prediction_result is prediction_result
    assert result.pred_mat_result is pred_mat_result
    assert calls[0][0] == "validate"
    assert calls[0][1]["phospho_matrix"] is phospho_matrix
    assert calls[0][1]["score_threshold"] == 0.7
    assert calls[0][1]["inclusion"] == 11
    assert calls[1] == ("execute", request)


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
    pred_mat_calls: list[dict[str, object]] = []
    analyzer_calls: list[dict[str, object]] = []
    reference_bundle = _ReferenceBundleDouble.__new__(_ReferenceBundleDouble)
    pred_mat_result = object()
    workflow_result = SimpleNamespace(
        pred_mat_result=pred_mat_result,
        scoring_result=object(),
        prediction_result=object(),
        close=lambda: None,
    )
    kinase_activity_result = object()
    workflow = SimpleKinaseWorkflow(
        reference_provider=_ProviderDouble(
            bundle=reference_bundle,
            calls=provider_calls,
        ),
        activity_analyzer=_AnalyzerDouble(
            result=kinase_activity_result,
            calls=analyzer_calls,
        ),
        analysis_ready_builder=_BuilderDouble(
            dataset=analysis_ready_dataset,
            calls=builder_calls,
        ),
    )
    workflow.pred_mat_workflow = _PredMatWorkflowDouble(
        result=workflow_result,
        calls=pred_mat_calls,
    )

    result = workflow.run(
        phospho=phospho,
        total=total,
        species="human",
        reference="ochoa",
        phospho_encoding="utf-8",
        min_observed=5,
        max_unmatched_fraction=0.25,
        min_substrates=3,
        score_threshold=0.65,
        kinase_activity_threshold=0.55,
        kinase_activity_min_substrates=4,
        kinase_activity_top_n_substrates=12,
    )

    assert isinstance(result, SimpleKinaseWorkflowResult)
    assert result.analysis_ready_dataset is analysis_ready_dataset
    assert result.reference_bundle is reference_bundle
    assert result.workflow_result is workflow_result
    assert result.kinase_activity_result is kinase_activity_result
    assert builder_calls == [
        {
            "phospho": phospho,
            "total": total,
            "phospho_encoding": "utf-8",
            "schema": builder_calls[0]["schema"],
            "comparisons": None,
            "preprocessing_config": None,
            "localization_threshold": 0.75,
            "min_observed": 5,
            "max_unmatched_fraction": 0.25,
            "total_sentinel": DEFAULT_TOTAL_SENTINEL,
            "phospho_sentinel": DEFAULT_PHOSPHO_SENTINEL,
            "source": "simple kinase workflow",
            "phospho_only_source": "simple kinase workflow (phospho only)",
        }
    ]
    assert provider_calls == [{"species": "human", "reference": "ochoa"}]
    assert pred_mat_calls == [
        {
            "phospho_matrix": phospho_matrix,
            "site_sequences": site_sequences,
            "reference_bundle": reference_bundle,
            "min_substrates": 3,
            "min_motif_size": 1,
            "allow_profile_only_fallback": False,
            "ensemble_size": 10,
            "top": 50,
            "score_threshold": 0.65,
            "inclusion": 20,
            "n_iterations": 5,
            "random_state": None,
            "svm_mode": None,
        }
    ]
    assert analyzer_calls == [
        {
            "pred_mat": pred_mat_result,
            "phospho_matrix": phospho_matrix,
            "threshold": 0.55,
            "min_substrates": 4,
            "top_n_substrates": 12,
        }
    ]


def test_simple_kinase_workflow_uses_builder_for_phospho_only_mode() -> None:
    builder_calls: list[dict[str, object]] = []
    provider_calls: list[dict[str, object]] = []
    pred_mat_calls: list[dict[str, object]] = []
    analyzer_calls: list[dict[str, object]] = []
    analysis_ready_dataset = SimpleNamespace(
        phospho_matrix=make_small_matrix(),
        site_sequences=pd.Series({"SITE_1": "AAAA", "SITE_2": "BBBB"}),
    )
    reference_bundle = _ReferenceBundleDouble.__new__(_ReferenceBundleDouble)
    workflow = SimpleKinaseWorkflow(
        reference_provider=_ProviderDouble(reference_bundle, provider_calls),
        activity_analyzer=_AnalyzerDouble(object(), analyzer_calls),
        analysis_ready_builder=_BuilderDouble(analysis_ready_dataset, builder_calls),
    )
    workflow.pred_mat_workflow = _PredMatWorkflowDouble(
        SimpleNamespace(
            pred_mat_result=object(),
            scoring_result=object(),
            prediction_result=object(),
            close=lambda: None,
        ),
        pred_mat_calls,
    )

    workflow.run(phospho=pd.DataFrame({"uid": ["u1"]}), species="mouse")

    assert builder_calls[0]["total"] is None
    assert builder_calls[0]["source"] == "simple kinase workflow"
    assert (
        builder_calls[0]["phospho_only_source"]
        == "simple kinase workflow (phospho only)"
    )
    assert len(provider_calls) == 1
    assert len(pred_mat_calls) == 1
    assert len(analyzer_calls) == 1


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

    def fake_execute_validated_signalome_request(received_request: object) -> object:
        calls.append(("execute", received_request))
        return signalome_result

    monkeypatch.setattr(
        api_workflows_module,
        "validate_signalome_request",
        fake_validate_signalome_request,
    )
    monkeypatch.setattr(
        api_workflows_module,
        "execute_validated_signalome_request",
        fake_execute_validated_signalome_request,
    )

    result = workflow.run(
        scoring_result=scoring_result,
        prediction_result=prediction_result,
        expression_matrix=expression_matrix,
        kinases_of_interest=["PRKACA", "BTK"],
        kinase_network_threshold=0.8,
        signalome_cutoff=0.4,
    )

    assert result is signalome_result
    assert calls[0][0] == "validate"
    assert calls[0][1]["scoring_result"] is scoring_result
    assert calls[0][1]["prediction_result"] is prediction_result
    assert calls[0][1]["expression_matrix"] is expression_matrix
    assert calls[0][1]["kinases_of_interest"] == ["PRKACA", "BTK"]
    assert calls[0][1]["kinase_network_threshold"] == 0.8
    assert calls[0][1]["signalome_cutoff"] == 0.4
    assert calls[1] == ("execute", request)


def test_root_convenience_exports_remain_thin_aliases() -> None:
    assert phospy.KinaseWorkflow is api_workflows_module.KinaseWorkflow
    assert phospy.PredMatWorkflow is api_workflows_module.PredMatWorkflow
    assert phospy.SignalomeWorkflow is api_workflows_module.SignalomeWorkflow
    assert phospy.SimpleKinaseWorkflow is api_workflows_module.SimpleKinaseWorkflow
