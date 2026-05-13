from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy.api import (
    Contrast,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.differential.executor import DifferentialAnalysisExecutor
from phospy.differential.models import DifferentialAnalysisRequest as CoreDiffRequest
from phospy.errors import (
    PhosPyInputError,
    WorkflowBoundaryError,
    WorkflowValidationError,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)

NEGATIVE_FIXTURE_DIR = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "rewrite_parity"
    / "differential_contract_negative_cases"
)


def test_negative_case_fixtures_are_source_labelled() -> None:
    provenance = NEGATIVE_FIXTURE_DIR / "PROVENANCE.md"
    assert provenance.is_file()
    text = provenance.read_text(encoding="utf-8")
    assert "ADR-0019" in text
    assert "rank-deficient" in text.lower()
    assert "missing values" in text.lower()


def _dataset(
    *,
    samples: tuple[str, ...] = ("A_1", "A_2", "B_1", "B_2"),
) -> object:
    from phospy import AnalysisReadyPhosphoDataset

    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0, 0.9],
            "A_2": [1.2, 2.1, 1.1],
            "B_1": [2.0, 1.8, 0.8],
            "B_2": [2.2, 2.0, 1.0],
        },
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"], name="site_id"),
    ).loc[:, list(samples)]
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


def _workflow_request(
    *,
    design: ExperimentalDesign,
    samples: tuple[str, ...] = ("A_1", "A_2", "B_1", "B_2"),
    minimum_condition_replicates: int = 2,
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=_dataset(samples=samples),
        design=design,
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
        ),
        config=DifferentialAnalysisConfig(
            minimum_condition_replicates=minimum_condition_replicates
        ),
    )


def test_workflow_rejects_batch_modelling_in_current_release() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="B_1", condition="B", batch="batch_2"),
            SampleDesignRecord(sample_id="B_2", condition="B", batch="batch_2"),
        )
    )
    with pytest.raises(WorkflowValidationError, match="unsupported design features"):
        DifferentialAnalysisWorkflow().run(_workflow_request(design=design))


def test_workflow_rejects_block_or_paired_modelling_in_current_release() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", block="pair_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", block="pair_2"),
            SampleDesignRecord(sample_id="B_1", condition="B", block="pair_1"),
            SampleDesignRecord(sample_id="B_2", condition="B", block="pair_2"),
        )
    )
    with pytest.raises(
        WorkflowValidationError, match="blocking/paired differential modelling"
    ):
        DifferentialAnalysisWorkflow().run(_workflow_request(design=design))


def test_workflow_rejects_non_positive_residual_dof_for_small_n_design() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B"),
        )
    )
    with pytest.raises(
        WorkflowBoundaryError,
        match="differential.interpreter.residual_dof",
    ):
        DifferentialAnalysisWorkflow().run(
            _workflow_request(
                design=design,
                samples=("A_1", "B_1"),
                minimum_condition_replicates=1,
            )
        )


def test_core_executor_rejects_rank_deficient_design_matrix() -> None:
    matrix = pd.read_csv(NEGATIVE_FIXTURE_DIR / "rank_deficient_matrix.csv").set_index(
        "site_id"
    )
    design = pd.read_csv(NEGATIVE_FIXTURE_DIR / "rank_deficient_design.csv").set_index(
        "sample"
    )
    contrasts = pd.read_csv(
        NEGATIVE_FIXTURE_DIR / "rank_deficient_contrasts.csv"
    ).set_index("coefficient")
    request = CoreDiffRequest(matrix=matrix, design=design, contrasts=contrasts)

    with pytest.raises(PhosPyInputError, match="full column rank"):
        DifferentialAnalysisExecutor().run(request)
