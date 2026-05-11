from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    Organism,
)
from phospy.errors import WorkflowBoundaryError, WorkflowValidationError
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


def _design() -> pd.DataFrame:
    return pd.DataFrame(
        {"A": [1.0, 1.0, 0.0, 0.0], "B": [0.0, 0.0, 1.0, 1.0]},
        index=pd.Index(["A_1", "A_2", "B_1", "B_2"], name="sample"),
    )


def _contrasts() -> pd.DataFrame:
    return pd.DataFrame(
        {"B_vs_A": [-1.0, 1.0]},
        index=pd.Index(["A", "B"], name="coefficient"),
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
                contrasts=_contrasts().rename(index={"A": "A_bad"}),
            )
        )


def test_differential_workflow_misaligned_design_fails_before_execution() -> None:
    design = _design().copy()
    design.index = pd.Index(["X1", "X2", "X3", "X4"], name="sample")
    with pytest.raises(WorkflowBoundaryError):
        DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=_dataset(),
                design=design,
                contrasts=_contrasts(),
            )
        )
