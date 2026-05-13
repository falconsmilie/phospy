from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    Contrast,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
    TechnicalReplicatePolicy,
)
from phospy.errors import WorkflowValidationError
from phospy.science.differential.executor import (
    DifferentialAnalysisExecutor as DifferentialComputationExecutor,
)
from phospy.science.differential.models import (
    DifferentialAnalysisRequest as ComputationRequest,
)
from phospy.workflows.differential.executor import DifferentialAnalysisExecutor
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)


def _dataset_with_technical_replicates() -> AnalysisReadyPhosphoDataset:
    phospho = pd.DataFrame(
        {
            "A1_T1": [1.0, 10.0],
            "A1_T2": [3.0, 8.0],
            "A2_T1": [2.0, 2.0],
            "A2_T2": [4.0, 4.0],
            "B1_T1": [5.0, 20.0],
            "B1_T2": [7.0, 18.0],
            "B2_T1": [6.0, 6.0],
            "B2_T2": [8.0, 8.0],
        },
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "S9"]
            ],
            "protein_id": ["MAPK14", "GSK3B"],
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


def _repeated_design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A1_T1",
                condition="A",
                biological_replicate_id="A1",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="A1_T2",
                condition="A",
                biological_replicate_id="A1",
                technical_replicate_id="T2",
            ),
            SampleDesignRecord(
                sample_id="A2_T1",
                condition="A",
                biological_replicate_id="A2",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="A2_T2",
                condition="A",
                biological_replicate_id="A2",
                technical_replicate_id="T2",
            ),
            SampleDesignRecord(
                sample_id="B1_T1",
                condition="B",
                biological_replicate_id="B1",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="B1_T2",
                condition="B",
                biological_replicate_id="B1",
                technical_replicate_id="T2",
            ),
            SampleDesignRecord(
                sample_id="B2_T1",
                condition="B",
                biological_replicate_id="B2",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="B2_T2",
                condition="B",
                biological_replicate_id="B2",
                technical_replicate_id="T2",
            ),
        )
    )


def _independent_design() -> ExperimentalDesign:
    return ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A1_T1",
                condition="A",
                biological_replicate_id="A1",
            ),
            SampleDesignRecord(
                sample_id="A1_T2",
                condition="A",
                biological_replicate_id="A2",
            ),
            SampleDesignRecord(
                sample_id="B1_T1",
                condition="B",
                biological_replicate_id="B1",
            ),
            SampleDesignRecord(
                sample_id="B1_T2",
                condition="B",
                biological_replicate_id="B2",
            ),
        )
    )


def _contrasts() -> tuple[Contrast, ...]:
    return (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )


def _request(
    *,
    dataset: AnalysisReadyPhosphoDataset | None = None,
    design: ExperimentalDesign | None = None,
    technical_replicate_policy: TechnicalReplicatePolicy = (
        TechnicalReplicatePolicy.REJECT
    ),
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=_dataset_with_technical_replicates() if dataset is None else dataset,
        design=_repeated_design() if design is None else design,
        contrasts=_contrasts(),
        config=DifferentialAnalysisConfig(
            technical_replicate_policy=technical_replicate_policy
        ),
    )


def test_independent_biological_replicates_pass_unchanged() -> None:
    dataset = _dataset_with_technical_replicates()
    request = DifferentialAnalysisRequest(
        dataset=dataset,
        design=_independent_design(),
        contrasts=_contrasts(),
        config=DifferentialAnalysisConfig(allow_design_subset=True),
    )
    validated = DifferentialAnalysisValidator().run(request)
    assert validated.analysis_sample_ids == ("A1_T1", "A1_T2", "B1_T1", "B1_T2")
    assert validated.workflow_provenance is None


def test_repeated_biological_replicates_fail_with_default_reject_policy() -> None:
    with pytest.raises(
        WorkflowValidationError,
        match="Technical replicates require explicit aggregation",
    ):
        DifferentialAnalysisValidator().run(_request())


def test_repeated_biological_replicates_aggregate_with_mean() -> None:
    validated = DifferentialAnalysisValidator().run(
        _request(technical_replicate_policy=TechnicalReplicatePolicy.MEAN)
    )
    aggregated = validated.dataset.phospho
    assert aggregated.columns.tolist() == ["A1", "A2", "B1", "B2"]
    assert aggregated.loc["MAPK14;Y182;", "A1"] == pytest.approx(2.0)
    assert aggregated.loc["MAPK14;Y182;", "A2"] == pytest.approx(3.0)
    assert aggregated.loc["MAPK14;Y182;", "B1"] == pytest.approx(6.0)
    assert aggregated.loc["MAPK14;Y182;", "B2"] == pytest.approx(7.0)


def test_repeated_biological_replicates_aggregate_with_median() -> None:
    validated = DifferentialAnalysisValidator().run(
        _request(technical_replicate_policy=TechnicalReplicatePolicy.MEDIAN)
    )
    aggregated = validated.dataset.phospho
    assert aggregated.columns.tolist() == ["A1", "A2", "B1", "B2"]
    assert aggregated.loc["GSK3B;S9;", "A1"] == pytest.approx(9.0)
    assert aggregated.loc["GSK3B;S9;", "B1"] == pytest.approx(19.0)


def test_aggregation_groups_by_condition_plus_biological_replicate_id() -> None:
    phospho = pd.DataFrame(
        {
            "A_R1_T1": [1.0, 2.0],
            "A_R1_T2": [3.0, 4.0],
            "B_R1_T1": [10.0, 20.0],
            "B_R1_T2": [30.0, 40.0],
        },
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "S9"]
            ],
            "protein_id": ["MAPK14", "GSK3B"],
        },
        index=phospho.index.copy(),
    )
    dataset = AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(
                sample_id="A_R1_T1",
                condition="A",
                biological_replicate_id="R1",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="A_R1_T2",
                condition="A",
                biological_replicate_id="R1",
                technical_replicate_id="T2",
            ),
            SampleDesignRecord(
                sample_id="B_R1_T1",
                condition="B",
                biological_replicate_id="R1",
                technical_replicate_id="T1",
            ),
            SampleDesignRecord(
                sample_id="B_R1_T2",
                condition="B",
                biological_replicate_id="R1",
                technical_replicate_id="T2",
            ),
        )
    )
    validated = DifferentialAnalysisValidator().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=design,
            contrasts=_contrasts(),
            config=DifferentialAnalysisConfig(
                technical_replicate_policy=TechnicalReplicatePolicy.MEAN,
                minimum_condition_replicates=1,
            ),
        )
    )
    assert validated.dataset.phospho.columns.tolist() == ["A__R1", "B__R1"]
    sample_ids = tuple(record.sample_id for record in validated.design.samples)
    assert sample_ids == ("A__R1", "B__R1")


def test_original_dataset_is_not_mutated_by_aggregation() -> None:
    dataset = _dataset_with_technical_replicates()
    before = dataset.phospho.copy(deep=True)
    DifferentialAnalysisValidator().run(
        _request(
            dataset=dataset,
            technical_replicate_policy=TechnicalReplicatePolicy.MEAN,
        )
    )
    pdt.assert_frame_equal(before, dataset.phospho)


def test_design_after_aggregation_has_one_row_per_biological_replicate() -> None:
    validated = DifferentialAnalysisValidator().run(
        _request(technical_replicate_policy=TechnicalReplicatePolicy.MEAN)
    )
    assert tuple(record.sample_id for record in validated.design.samples) == (
        "A1",
        "A2",
        "B1",
        "B2",
    )
    assert all(
        record.technical_replicate_id is None for record in validated.design.samples
    )


def test_executor_receives_aggregated_matrix() -> None:
    observed_columns: list[str] = []

    class _ComputationExecutorSpy:
        def run(self, request: ComputationRequest):
            observed_columns.extend(request.matrix.columns.astype(str).tolist())
            return DifferentialComputationExecutor().run(request)

    workflow = DifferentialAnalysisWorkflow(
        executor=DifferentialAnalysisExecutor(
            computation_executor=_ComputationExecutorSpy()  # type: ignore[arg-type]
        )
    )
    workflow.run(_request(technical_replicate_policy=TechnicalReplicatePolicy.MEAN))
    assert observed_columns == ["A1", "A2", "B1", "B2"]


def test_provenance_records_technical_replicate_lineage() -> None:
    result = DifferentialAnalysisWorkflow().run(
        _request(technical_replicate_policy=TechnicalReplicatePolicy.MEAN)
    )
    assert result.workflow_provenance is not None
    assert result.workflow_provenance["technical_replicate_policy"] == "mean"
    groups = result.workflow_provenance["groups"]
    assert isinstance(groups, list)
    a1_group = next(
        group
        for group in groups
        if group["condition"] == "A" and group["biological_replicate_id"] == "A1"
    )
    assert a1_group["output_sample_id"] == "A1"
    assert a1_group["input_sample_ids"] == ["A1_T1", "A1_T2"]
    assert a1_group["technical_replicate_ids"] == ["T1", "T2"]
    assert a1_group["n_technical_replicates"] == 2


def test_invalid_technical_replicate_policy_fails() -> None:
    request = DifferentialAnalysisRequest(
        dataset=_dataset_with_technical_replicates(),
        design=_repeated_design(),
        contrasts=_contrasts(),
        config=DifferentialAnalysisConfig(
            technical_replicate_policy="invalid"  # type: ignore[arg-type]
        ),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="technical_replicate_policy must be TechnicalReplicatePolicy",
    ):
        DifferentialAnalysisValidator().run(request)
