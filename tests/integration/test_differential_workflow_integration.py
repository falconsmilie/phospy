from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.errors import WorkflowValidationError
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
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
            "site_sequence": ["A" * 31, "B" * 31, "C" * 31, "D" * 31],
            "protein_id": ["MAPK14", "GSK3B", "AKT1", "RPS6KB1"],
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
