from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import pandas as pd

import phospy
import phospy.api.signalome_workflows as signalome_workflows_module
from phospy.api import (
    DatasetLoadOptions,
    KinaseActivityConfig,
    PredictionRunConfig,
    SignalomeRunConfig,
    SignalomeWorkflow,
    SimpleKinaseWorkflow,
)
from phospy.api.workflow_results import SimpleKinaseWorkflowResult
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
class _PredMatWorkflowDouble:
    result: object
    calls: list[dict[str, object]]

    def run(self, **kwargs: object) -> object:
        self.calls.append(kwargs)
        return self.result


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
    scoring_result = object()
    prediction_result = object()
    workflow_result = SimpleNamespace(
        pred_mat_result=pred_mat_result,
        scoring_result=scoring_result,
        prediction_result=prediction_result,
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
    assert pred_mat_calls == [
        {
            "phospho_matrix": phospho_matrix,
            "site_sequences": site_sequences,
            "reference_bundle": reference_bundle,
            "prediction_config": PredictionRunConfig(
                min_substrates=3,
                score_threshold=0.65,
            ),
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


def test_simple_kinase_workflow_run_delegates_to_execution_service() -> None:
    workflow = SimpleKinaseWorkflow()
    expected_result = object()
    calls: list[dict[str, object]] = []

    class _ExecutionServiceDouble:
        def run(self, **kwargs: object) -> object:
            calls.append(kwargs)
            return expected_result

    workflow._execution_service = _ExecutionServiceDouble()
    result = workflow.run(
        phospho=pd.DataFrame({"uid": ["u1"]}),
        species="rat",
    )

    assert result is expected_result
    assert len(calls) == 1
    assert calls[0]["analysis_ready_builder"] is workflow.analysis_ready_builder
    assert calls[0]["reference_provider"] is workflow.reference_provider
    assert calls[0]["pred_mat_workflow"] is workflow.pred_mat_workflow
    assert calls[0]["activity_analyzer"] is workflow.activity_analyzer


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
