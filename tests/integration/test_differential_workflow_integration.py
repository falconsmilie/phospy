from __future__ import annotations

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DifferentialAnalysisWorkflow,
)
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.errors import WorkflowValidationError
from phospy.science.transformations.models import (
    IntensityScaleState,
    MatrixIntensityScaleState,
    QuantitativeMeaning,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)


def _dataset() -> AnalysisReadyPhosphoDataset:
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0, 1.0, 1.2],
            "A_2": [1.1, 2.1, 1.1, 1.1],
            "B_1": [2.1, 2.0, 1.0, 1.3],
            "B_2": [2.0, 2.2, 0.9, 1.2],
        },
        index=pd.Index(
            ["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;", "RPS6KB1;T389;"],
            name="site_id",
        ),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1", "RPS6KB1"],
            "site": ["Y182", "S9", "T308", "T389"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "S9", "T308", "T389"]
            ],
            "localisation_confidence": [0.95] * phospho.shape[0],
            "protein_id": ["MAPK14", "GSK3B", "AKT1", "RPS6KB1"],
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


def _design() -> ExperimentalDesign:
    return ExperimentalDesign(
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
    )


def _contrasts() -> tuple[Contrast, ...]:
    return (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
    )


def test_differential_workflow_runs_on_analysis_ready_dataset() -> None:
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=_dataset(),
            design=_design(),
            contrasts=_contrasts(),
        )
    )

    table = result.table_for("B_vs_A")
    assert list(table.columns) == ["logFC", "t", "P.Value", "adj.P.Val"]
    assert table.index.tolist() == _dataset().phospho.index.tolist()
    assert (table.loc[:, "adj.P.Val"] >= 0.0).all()
    assert (table.loc[:, "adj.P.Val"] <= 1.0).all()


def test_differential_workflow_invalid_contrast_fails_before_execution() -> None:
    with pytest.raises(WorkflowValidationError):
        DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=_design(),
                contrasts=(
                    Contrast(
                        name="B_vs_A",
                        numerator_condition="B",
                        denominator_condition="A_bad",
                    ),
                ),
            )
        )


def test_differential_workflow_misaligned_design_fails_before_execution() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="X1", condition="A"),
            SampleDesignRecord(sample_id="X2", condition="A"),
            SampleDesignRecord(sample_id="X3", condition="B"),
            SampleDesignRecord(sample_id="X4", condition="B"),
        )
    )
    with pytest.raises(WorkflowValidationError):
        DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=design,
                contrasts=_contrasts(),
            )
        )


def test_differential_workflow_runs_on_builder_log2_dataset() -> None:
    phospho = pd.DataFrame(
        {
            "A_1": [10.0, 8.0, 6.0, 4.0],
            "A_2": [11.0, 8.5, 5.5, 4.5],
            "B_1": [15.0, 9.0, 6.5, 4.8],
            "B_2": [16.0, 9.2, 6.6, 5.0],
        },
        index=pd.Index(
            ["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;", "RPS6KB1;T389;"],
            name="site_id",
        ),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1", "RPS6KB1"],
            "site": ["Y182", "S9", "T308", "T389"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "S9", "T308", "T389"]
            ],
            "localisation_confidence": [0.95] * phospho.shape[0],
            "protein_id": ["MAPK14", "GSK3B", "AKT1", "RPS6KB1"],
        },
        index=phospho.index.copy(),
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                )
            ),
        )
    )

    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
            dataset=dataset,
            design=_design(),
            contrasts=_contrasts(),
        )
    )
    assert "B_vs_A" in result.contrast_tables


def test_differential_workflow_rejects_linear_scale_before_execution() -> None:
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
        DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=linear_dataset,
                design=_design(),
                contrasts=_contrasts(),
            )
        )


def test_differential_workflow_rejects_unestablished_log2_scale_before_execution() -> (
    None
):
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
        DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=dataset,
                design=_design(),
                contrasts=_contrasts(),
            )
        )


def test_differential_workflow_rejects_unknown_scale_before_execution() -> None:
    dataset = _dataset()
    unknown_state = IntensityScaleState.raw(
        has_total_matrix=False
    ).with_quantitative_meaning(QuantitativeMeaning.UNKNOWN)
    object.__setattr__(dataset, "intensity_scale_state", unknown_state)

    with pytest.raises(
        WorkflowValidationError,
        match="requires established log2-scale phospho intensities",
    ):
        DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=dataset,
                design=_design(),
                contrasts=_contrasts(),
            )
        )
