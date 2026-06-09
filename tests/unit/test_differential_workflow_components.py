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
from phospy.science.datasets.models import DatasetPreprocessingReport
from phospy.science.differential.models import EmpiricalBayesConfig
from phospy.science.transformations.models import (
    IntensityScaleState,
    MatrixIntensityScaleState,
)
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
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns


def _dataset() -> AnalysisReadyPhosphoDataset:
    genes = ["MAPK14", "GSK3B", "AKT1"]
    sites = ["Y182", "S9", "T308"]
    site_index = protein_site_key_index(protein_identifiers=genes, sites=sites)
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0, 1.0],
            "A_2": [1.1, 2.1, 1.1],
            "B_1": [2.1, 2.0, 1.0],
            "B_2": [2.0, 2.2, 0.9],
        },
        index=site_index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": ["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"],
            **site_key_context_columns(site_index),
            "gene_symbol": genes,
            "site": sites,
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15) for site in sites
            ],
            "protein_id": genes,
        },
        index=phospho.index.copy(),
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
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


def test_differential_workflow_uses_explicit_design_not_sample_metadata_conditions() -> (
    None
):
    base_request = _request()
    base_dataset = base_request.dataset
    dataset_with_passive_metadata = AnalysisReadyPhosphoDataset(
        phospho=base_dataset.phospho,
        site_metadata=base_dataset.site_metadata,
        sample_metadata=pd.DataFrame(
            {
                "condition": ["metadata_only"] * 4,
                "batch": ["batch_1", "batch_1", "batch_2", "batch_2"],
            },
            index=base_dataset.phospho.columns.copy(),
        ),
        organism=base_dataset.organism,
        intensity_scale_state=base_dataset.intensity_scale_state,
        processing_state=base_dataset.processing_state,
    )

    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset_with_passive_metadata,
            design=base_request.design,
            contrasts=base_request.contrasts,
        )
    )

    assert isinstance(result, DifferentialAnalysisResult)


def test_differential_result_references_input_dataset_preprocessing_report() -> None:
    base_dataset = _dataset()
    preprocessing_report = DatasetPreprocessingReport.from_rows()
    dataset_with_report = AnalysisReadyPhosphoDataset(
        phospho=base_dataset.phospho,
        site_metadata=base_dataset.site_metadata,
        sample_metadata=base_dataset.sample_metadata,
        total=base_dataset.total,
        comparisons=base_dataset.comparisons,
        organism=base_dataset.organism,
        intensity_scale_state=base_dataset.intensity_scale_state,
        processing_state=base_dataset.processing_state,
        preprocessing_report=preprocessing_report,
        provenance=base_dataset.provenance,
    )
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset_with_report,
            design=_request().design,
            contrasts=_request().contrasts,
        )
    )
    assert result.input_dataset_preprocessing_report is preprocessing_report


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


def test_differential_validator_rejects_established_linear_scale() -> None:
    valid_dataset = _dataset()
    linear_dataset = AnalysisReadyPhosphoDataset(
        phospho=valid_dataset.phospho,
        site_metadata=valid_dataset.site_metadata,
        organism=valid_dataset.organism,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="requires established log2-scale phospho intensities",
    ):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=linear_dataset,
                design=_request().design,
                contrasts=_request().contrasts,
            )
        )


def test_differential_validator_rejects_declared_but_unestablished_log2_scale() -> None:
    dataset = _dataset()
    object.__setattr__(
        dataset,
        "intensity_scale_state",
        IntensityScaleState(
            phospho=MatrixIntensityScaleState.log2(established_by="test.declaration"),
            total=None,
        ),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="requires established log2-scale phospho intensities",
    ):
        DifferentialAnalysisValidator().run(
            DifferentialAnalysisRequest(
                dataset=dataset,
                design=_request().design,
                contrasts=_request().contrasts,
            )
        )


def test_differential_invalid_scale_fails_before_executor() -> None:
    calls = {"executor": 0}

    class _ExecutorSpy:
        def run(self, request: InterpretedDifferentialAnalysisRequest):
            calls["executor"] += 1
            raise AssertionError("executor should not be called")

    valid_dataset = _dataset()
    linear_dataset = AnalysisReadyPhosphoDataset(
        phospho=valid_dataset.phospho,
        site_metadata=valid_dataset.site_metadata,
        organism=valid_dataset.organism,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )
    with pytest.raises(WorkflowValidationError):
        DifferentialAnalysisWorkflow(executor=_ExecutorSpy()).run(  # type: ignore[arg-type]
            DifferentialAnalysisRequest(
                dataset=linear_dataset,
                design=_request().design,
                contrasts=_request().contrasts,
            )
        )
    assert calls["executor"] == 0


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

    class _FakeTechnicalReplicatePlanner:
        def run(self, *, dataset, design, technical_replicate_policy):
            calls["technical_replicate_policy"] = technical_replicate_policy
            from phospy.workflows.differential.replicates import (
                TechnicalReplicateAggregationPlan,
            )

            return TechnicalReplicateAggregationPlan(
                technical_replicate_policy=technical_replicate_policy,
                groups=(),
                aggregate_phospho=False,
                aggregate_total_protein=False,
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
        technical_replicate_planner=_FakeTechnicalReplicatePlanner(),  # type: ignore[arg-type]
        design_validator=_FakeDesignValidator(),  # type: ignore[arg-type]
    ).run(request)

    assert calls["technical_replicate_policy"] is TechnicalReplicatePolicy.MEAN
    assert calls["allow_design_subset"] is True
    assert calls["minimum_condition_replicates"] == 1
