from __future__ import annotations

import inspect
from typing import Any

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    Contrast,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    MultipleTestingConfig,
    Organism,
    SampleDesignRecord,
    TechnicalReplicatePolicy,
)
from phospy.api.results import DifferentialAnalysisResult
from phospy.errors import WorkflowBoundaryError, WorkflowValidationError
from phospy.science.differential.models import EmpiricalBayesConfig
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
        design=ExperimentalDesign(
            samples=(
                SampleDesignRecord(
                    sample_id="A_1",
                    condition="A",
                    biological_replicate_id="A_r1",
                ),
                SampleDesignRecord(
                    sample_id="A_2",
                    condition="A",
                    biological_replicate_id="A_r2",
                ),
                SampleDesignRecord(
                    sample_id="B_1",
                    condition="B",
                    biological_replicate_id="B_r1",
                ),
                SampleDesignRecord(
                    sample_id="B_2",
                    condition="B",
                    biological_replicate_id="B_r2",
                ),
            )
        ),
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
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
    bad_contrasts = (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A_bad",
        ),
    )
    bad_request = DifferentialAnalysisRequest(
        dataset=request.dataset,
        design=request.design,
        contrasts=bad_contrasts,
    )
    with pytest.raises(
        WorkflowValidationError,
        match="unknown denominator condition",
    ):
        DifferentialAnalysisValidator().run(bad_request)


def test_differential_interpreter_checks_sample_to_design_alignment() -> None:
    request = _request()
    validated = DifferentialAnalysisValidator().run(request)
    misaligned_design = validated.design_matrix.to_dataframe()
    misaligned_design.index = pd.Index(["x1", "x2", "x3", "x4"], name="sample")
    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.interpreter.sample_label_alignment",
    ):
        DifferentialAnalysisInterpreter().run(
            ValidatedDifferentialAnalysisRequest(
                dataset=validated.dataset,
                design=validated.design,
                contrasts=validated.contrasts,
                analysis_sample_ids=validated.analysis_sample_ids,
                design_matrix=type(validated.design_matrix)(misaligned_design),
                contrast_matrix=validated.contrast_matrix,
                config=validated.config,
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
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A_bad",
            ),
        ),
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
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="x1", condition="A"),
            SampleDesignRecord(sample_id="x2", condition="A"),
            SampleDesignRecord(sample_id="x3", condition="B"),
            SampleDesignRecord(sample_id="x4", condition="B"),
        )
    )

    workflow = DifferentialAnalysisWorkflow(executor=_ExecutorSpy())  # type: ignore[arg-type]
    with pytest.raises(WorkflowValidationError):
        workflow.run(
            DifferentialAnalysisRequest(
                dataset=request.dataset,
                design=design,
                contrasts=request.contrasts,
            )
        )
    assert calls["executor"] == 0


def test_differential_validator_rejects_non_differential_config_type() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="differential workflow request config must be DifferentialAnalysisConfig",
    ):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=_request().design,
                contrasts=_request().contrasts,
                config=object(),  # type: ignore[arg-type]
            )
        )


def test_differential_validator_rejects_raw_string_technical_replicate_policy() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="technical_replicate_policy must be TechnicalReplicatePolicy",
    ):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=_request().design,
                contrasts=_request().contrasts,
                config=DifferentialAnalysisConfig(
                    technical_replicate_policy="mean",  # type: ignore[arg-type]
                ),
            )
        )


def test_differential_validator_rejects_non_bool_allow_design_subset() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="allow_design_subset must be a bool",
    ):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=_request().design,
                contrasts=_request().contrasts,
                config=DifferentialAnalysisConfig(
                    allow_design_subset=1,  # type: ignore[arg-type]
                ),
            )
        )


@pytest.mark.parametrize("value", [0, -1, 1.5, True])
def test_differential_validator_rejects_invalid_minimum_condition_replicates(
    value: object,
) -> None:
    pattern = (
        "minimum_condition_replicates must be >= 1"
        if isinstance(value, int) and not isinstance(value, bool)
        else "minimum_condition_replicates must be an int"
    )
    with pytest.raises(WorkflowValidationError, match=pattern):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=_request().design,
                contrasts=_request().contrasts,
                config=DifferentialAnalysisConfig(
                    minimum_condition_replicates=value,  # type: ignore[arg-type]
                ),
            )
        )


def test_differential_validator_rejects_non_empirical_bayes_config() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="empirical_bayes must be EmpiricalBayesConfig",
    ):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=_request().design,
                contrasts=_request().contrasts,
                config=DifferentialAnalysisConfig(
                    empirical_bayes=object(),  # type: ignore[arg-type]
                ),
            )
        )


def test_differential_validator_rejects_non_multiple_testing_config() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="multiple_testing must be MultipleTestingConfig",
    ):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=_request().design,
                contrasts=_request().contrasts,
                config=DifferentialAnalysisConfig(
                    multiple_testing=object(),  # type: ignore[arg-type]
                ),
            )
        )


def test_differential_validator_passes_config_values_to_collaborators() -> None:
    calls: dict[str, object] = {}
    base_request = _request()

    class _FakeTechnicalReplicateResolver:
        def run(self, *, dataset, design, technical_replicate_policy):
            calls["technical_replicate_policy"] = technical_replicate_policy
            from phospy.workflows.differential.replicates import (
                TechnicalReplicateResolution,
            )

            return TechnicalReplicateResolution(
                dataset=dataset,
                design=design,
                workflow_provenance=None,
            )

    class _FakeDesignValidator:
        def run(
            self,
            *,
            dataset,
            design,
            contrasts,
            allow_design_subset,
            minimum_condition_replicates,
        ):
            calls["allow_design_subset"] = allow_design_subset
            calls["minimum_condition_replicates"] = minimum_condition_replicates
            from phospy.validation.workflows.differential import (
                ValidatedExperimentalDesignContract,
            )

            return ValidatedExperimentalDesignContract(
                design=design,
                contrasts=contrasts,
                analysis_sample_ids=tuple(
                    record.sample_id for record in design.samples
                ),
                condition_labels=("A", "B"),
                design_frame=pd.DataFrame(
                    {
                        "A": [1.0, 1.0, 0.0, 0.0],
                        "B": [0.0, 0.0, 1.0, 1.0],
                    },
                    index=pd.Index(["A_1", "A_2", "B_1", "B_2"], name="sample"),
                ),
                contrast_frame=pd.DataFrame(
                    {"B_vs_A": [-1.0, 1.0]},
                    index=pd.Index(["A", "B"], name="coefficient"),
                ),
            )

    request = DifferentialAnalysisRequest(
        dataset=base_request.dataset,
        design=base_request.design,
        contrasts=base_request.contrasts,
        config=DifferentialAnalysisConfig(
            technical_replicate_policy=TechnicalReplicatePolicy.MEAN,
            allow_design_subset=True,
            minimum_condition_replicates=1,
            empirical_bayes=EmpiricalBayesConfig(method="standard"),
            multiple_testing=MultipleTestingConfig(),
        ),
    )
    DifferentialAnalysisValidator(
        technical_replicate_resolver=_FakeTechnicalReplicateResolver(),  # type: ignore[arg-type]
        design_validator=_FakeDesignValidator(),  # type: ignore[arg-type]
    ).run(request)

    assert calls["technical_replicate_policy"] is TechnicalReplicatePolicy.MEAN
    assert calls["allow_design_subset"] is True
    assert calls["minimum_condition_replicates"] == 1
