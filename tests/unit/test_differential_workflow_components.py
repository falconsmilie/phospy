from __future__ import annotations

import inspect
from typing import Any

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    Organism,
)
from phospy.api.results import DifferentialAnalysisResult
from phospy.errors import WorkflowBoundaryError, WorkflowValidationError
from phospy.workflows.differential.executor import DifferentialAnalysisExecutor
from phospy.workflows.differential.interpreter import DifferentialAnalysisInterpreter
from phospy.workflows.differential.models import (
    InterpretedDifferentialAnalysisRequest,
    ValidatedDifferentialAnalysisRequest,
)
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)


def _dataset() -> AnalysisReadyPhosphoDataset:
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0, 1.0],
            "A_2": [1.1, 2.1, 1.1],
            "B_1": [2.1, 2.0, 1.0],
            "B_2": [2.0, 2.2, 0.9],
        },
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1"],
            "site": ["Y182", "S9", "T308"],
            "site_sequence": ["A" * 31, "B" * 31, "C" * 31],
            "protein_id": ["MAPK14", "GSK3B", "AKT1"],
        },
        index=phospho.index.copy(),
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


def _request() -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=_dataset(),
        design=pd.DataFrame(
            {
                "A": [1.0, 1.0, 0.0, 0.0],
                "B": [0.0, 0.0, 1.0, 1.0],
            },
            index=pd.Index(["A_1", "A_2", "B_1", "B_2"], name="sample"),
        ),
        contrasts=pd.DataFrame(
            {"B_vs_A": [-1.0, 1.0]},
            index=pd.Index(["A", "B"], name="coefficient"),
        ),
    )


def _public_methods(cls: type[Any]) -> set[str]:
    return {
        name
        for name, value in inspect.getmembers(cls)
        if callable(value) and not name.startswith("_")
    }


def test_differential_workflow_calls_validator_interpreter_executor_in_order() -> None:
    events: list[str] = []
    validated = object()
    interpreted = object()
    expected_result = object()

    class _Validator:
        def run(self, request: object) -> object:
            events.append("validator")
            return validated

    class _Interpreter:
        def run(self, request: object) -> object:
            events.append("interpreter")
            assert request is validated
            return interpreted

    class _Executor:
        def run(self, request: object) -> object:
            events.append("executor")
            assert request is interpreted
            return expected_result

    result = DifferentialAnalysisWorkflow(
        validator=_Validator(),  # type: ignore[arg-type]
        interpreter=_Interpreter(),  # type: ignore[arg-type]
        executor=_Executor(),  # type: ignore[arg-type]
    ).run(_request())

    assert events == ["validator", "interpreter", "executor"]
    assert result is expected_result


def test_differential_workflow_dependency_injection_supports_real_stage_contracts() -> (
    None
):
    workflow = DifferentialAnalysisWorkflow(
        validator=DifferentialAnalysisValidator(),
        interpreter=DifferentialAnalysisInterpreter(),
        executor=DifferentialAnalysisExecutor(),
    )
    result = workflow.run(_request())
    assert isinstance(result, DifferentialAnalysisResult)


def test_differential_public_stages_expose_run_only() -> None:
    assert _public_methods(DifferentialAnalysisWorkflow) == {"run"}
    assert _public_methods(DifferentialAnalysisValidator) == {"run"}
    assert _public_methods(DifferentialAnalysisInterpreter) == {"run"}
    assert _public_methods(DifferentialAnalysisExecutor) == {"run"}
    assert not hasattr(DifferentialAnalysisWorkflow, "validate")
    assert not hasattr(DifferentialAnalysisWorkflow, "interpret")
    assert not hasattr(DifferentialAnalysisWorkflow, "execute")


def test_differential_validator_rejects_invalid_raw_request_type() -> None:
    validator = DifferentialAnalysisValidator()
    with pytest.raises(
        WorkflowValidationError,
        match="differential workflow input must be a DifferentialAnalysisRequest",
    ):
        validator.run(object())


def test_differential_validator_rejects_unknown_contrast_term_before_interpretation() -> (
    None
):
    request = _request()
    bad_contrasts = request.contrasts.to_dataframe().rename(index={"A": "A_bad"})
    bad_request = DifferentialAnalysisRequest(
        dataset=request.dataset,
        design=request.design,
        contrasts=bad_contrasts,
    )
    with pytest.raises(
        WorkflowValidationError,
        match="contrasts.index must match differential workflow request design.columns",
    ):
        DifferentialAnalysisValidator().run(bad_request)


def test_differential_interpreter_checks_sample_to_design_alignment() -> None:
    request = _request()
    validated = ValidatedDifferentialAnalysisRequest(
        dataset=request.dataset,
        design=request.design,
        contrasts=request.contrasts,
        empirical_bayes=request.empirical_bayes,
        multiple_testing=request.multiple_testing,
    )
    misaligned_design = validated.design.to_dataframe()
    misaligned_design.index = pd.Index(["x1", "x2", "x3", "x4"], name="sample")
    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.interpreter.sample_label_alignment",
    ):
        DifferentialAnalysisInterpreter().run(
            ValidatedDifferentialAnalysisRequest(
                dataset=validated.dataset,
                design=type(validated.design)(misaligned_design),
                contrasts=validated.contrasts,
                empirical_bayes=validated.empirical_bayes,
                multiple_testing=validated.multiple_testing,
            )
        )


def test_differential_executor_accepts_only_interpreted_requests() -> None:
    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.executor.interpreted_request_type",
    ):
        DifferentialAnalysisExecutor().run(object())  # type: ignore[arg-type]


def test_differential_invalid_contrast_fails_before_executor() -> None:
    calls = {"executor": 0}

    class _ExecutorSpy:
        def run(self, request: InterpretedDifferentialAnalysisRequest):
            calls["executor"] += 1
            raise AssertionError("executor should not be called")

    request = _request()
    bad_request = DifferentialAnalysisRequest(
        dataset=request.dataset,
        design=request.design,
        contrasts=request.contrasts.to_dataframe().rename(index={"A": "A_bad"}),
    )

    workflow = DifferentialAnalysisWorkflow(executor=_ExecutorSpy())  # type: ignore[arg-type]
    with pytest.raises(WorkflowValidationError):
        workflow.run(bad_request)
    assert calls["executor"] == 0


def test_differential_misaligned_design_fails_before_executor() -> None:
    calls = {"executor": 0}

    class _ExecutorSpy:
        def run(self, request: InterpretedDifferentialAnalysisRequest):
            calls["executor"] += 1
            raise AssertionError("executor should not be called")

    request = _request()
    design = request.design.to_dataframe()
    design.index = pd.Index(["x1", "x2", "x3", "x4"], name="sample")

    workflow = DifferentialAnalysisWorkflow(executor=_ExecutorSpy())  # type: ignore[arg-type]
    with pytest.raises(WorkflowBoundaryError):
        workflow.run(
            DifferentialAnalysisRequest(
                dataset=request.dataset,
                design=design,
                contrasts=request.contrasts,
            )
        )
    assert calls["executor"] == 0
