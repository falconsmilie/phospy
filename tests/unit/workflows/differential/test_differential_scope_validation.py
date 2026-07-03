from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    BatchCovariate,
    Contrast,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.errors import WorkflowBoundaryError, WorkflowValidationError
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT,
)
from phospy.science.transformations.models import (
    IntensityScaleState,
    MatrixIntensityScaleState,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.processing_state import (
    imputed_processing_state as valid_imputed_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns
from tests.support.unsafe_dataset_states import (
    unsafe_replace_dataset_intensity_scale_state,
)

ROOT = Path(__file__).resolve().parents[4]


class _ExecutorSpy:
    def __init__(self) -> None:
        self.calls = 0

    def run(self, request: object) -> object:
        del request
        self.calls += 1
        raise AssertionError("executor should not be called")


def _dataset(
    *,
    samples: tuple[str, ...] = ("A_1", "A_2", "B_1", "B_2"),
    matrix: pd.DataFrame | None = None,
) -> AnalysisReadyPhosphoDataset:
    genes = ["MAPK14", "GSK3B", "AKT1"]
    sites = ["Y182", "S9", "T308"]
    site_index = protein_site_key_index(protein_identifiers=genes, sites=sites)
    phospho = (
        pd.DataFrame(
            {
                "A_1": [1.0, 2.0, 1.0],
                "A_2": [1.2, 2.1, 1.1],
                "B_1": [2.0, 1.8, 0.8],
                "B_2": [2.2, 2.0, 1.0],
            },
            index=site_index.copy(),
        )
        if matrix is None
        else matrix
    ).loc[:, list(samples)]
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
        index=site_index.copy(),
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


def _condition_design(
    samples: tuple[str, ...] = ("A_1", "A_2", "B_1", "B_2"),
) -> ExperimentalDesign:
    records: list[SampleDesignRecord] = []
    for sample_id in samples:
        condition = sample_id.split("_", maxsplit=1)[0]
        records.append(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=condition,
                biological_replicate_id=f"{condition}_{sample_id}",
            )
        )
    return ExperimentalDesign(samples=tuple(records))


def _contrast() -> tuple[Contrast, ...]:
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
    contrasts: tuple[Contrast, ...] | None = None,
    minimum_condition_replicates: int = 2,
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=_dataset() if dataset is None else dataset,
        design=_condition_design() if design is None else design,
        contrasts=_contrast() if contrasts is None else contrasts,
        config=DifferentialAnalysisConfig(
            minimum_condition_replicates=minimum_condition_replicates
        ),
    )


def test_rank_deficient_design_fails_before_execution() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="B_1", condition="B", batch="batch_2"),
            SampleDesignRecord(sample_id="B_2", condition="B", batch="batch_2"),
        ),
        fixed_effects=(BatchCovariate(),),
    )
    executor = _ExecutorSpy()

    with pytest.raises(WorkflowValidationError, match="rank deficient"):
        DifferentialAnalysisWorkflow(executor=executor).run(_request(design=design))
    assert executor.calls == 0


def test_unknown_contrast_condition_fails_before_execution() -> None:
    executor = _ExecutorSpy()

    with pytest.raises(WorkflowValidationError, match="unknown denominator condition"):
        DifferentialAnalysisWorkflow(executor=executor).run(
            _request(
                contrasts=(
                    Contrast(
                        name="B_vs_missing",
                        numerator_condition="B",
                        denominator_condition="missing",
                    ),
                )
            )
        )
    assert executor.calls == 0


def test_insufficient_residual_degrees_of_freedom_fails_before_execution() -> None:
    executor = _ExecutorSpy()

    with pytest.raises(WorkflowBoundaryError, match="residual_dof"):
        DifferentialAnalysisWorkflow(executor=executor).run(
            _request(
                dataset=_dataset(samples=("A_1", "B_1")),
                design=_condition_design(samples=("A_1", "B_1")),
                minimum_condition_replicates=1,
            )
        )
    assert executor.calls == 0


def test_imputed_data_are_rejected_under_default_policy() -> None:
    dataset = _dataset()
    imputed_dataset = AnalysisReadyPhosphoDataset(
        phospho=dataset.phospho,
        site_metadata=dataset.site_metadata,
        organism=dataset.organism,
        intensity_scale_state=dataset.intensity_scale_state,
        processing_state=valid_imputed_processing_state(dataset.processing_state),
    )

    with pytest.raises(WorkflowValidationError, match="imputed cells"):
        DifferentialAnalysisWorkflow().run(_request(dataset=imputed_dataset))


def test_established_intensity_scale_is_required() -> None:
    dataset = _dataset()
    unsafe_replace_dataset_intensity_scale_state(
        dataset,
        IntensityScaleState(
            phospho=MatrixIntensityScaleState.log2(established_by="test.declaration"),
            total=None,
        ),
    )

    with pytest.raises(WorkflowValidationError, match="established log2-scale"):
        DifferentialAnalysisWorkflow().run(_request(dataset=dataset))


def test_all_constant_site_intensities_are_withheld_per_feature() -> None:
    dataset = _dataset()
    matrix = dataset.phospho
    matrix.iloc[0, :] = 5.0
    constant_dataset = _dataset(matrix=matrix)

    result = DifferentialAnalysisWorkflow().run(_request(dataset=constant_dataset))
    table = result.table_for("B_vs_A")

    assert table.iloc[0][DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT
    )
    assert table.iloc[0][["logFC", "t", "P.Value", "adj.P.Val"]].isna().all()
    assert table.iloc[1:][DIFFERENTIAL_RESULT_STATUS_COLUMN].tolist() == [
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        DIFFERENTIAL_RESULT_STATUS_TESTED,
    ]


def test_differential_docs_do_not_claim_broad_phosr_or_limma_parity() -> None:
    docs = {
        "parity": ROOT / "docs" / "parity.md",
        "coverage": ROOT / "docs" / "scientific-coverage.md",
        "api": ROOT / "docs" / "api" / "differential-analysis.md",
    }
    text = "\n".join(path.read_text(encoding="utf-8") for path in docs.values())
    normalized = " ".join(text.lower().split())

    assert "current differential analysis is not full phosr or limma parity" in (
        normalized
    )
    assert "supported designs are limited to tested design and contrast envelopes" in (
        normalized
    )
    assert "upstream-imputed datasets are rejected by default" in normalized
    assert "fixed-effect covariates are not full batch correction" in normalized
    forbidden_claims = (
        "differential analysis is full phosr parity",
        "differential analysis is full limma parity",
        "limma-equivalent differential analysis",
    )
    for claim in forbidden_claims:
        assert claim not in normalized
