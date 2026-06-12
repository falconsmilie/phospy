from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy.api import (
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    BatchCovariate,
    Contrast,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.errors import (
    PhosPyInputError,
    WorkflowBoundaryError,
    WorkflowValidationError,
)
from phospy.science.differential.executor import DifferentialAnalysisExecutor
from phospy.science.differential.models import (
    DifferentialAnalysisRequest as CoreDiffRequest,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

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

    genes = ["MAPK14", "GSK3B", "AKT1"]
    sites = ["Y182", "S9", "T308"]
    site_index = protein_site_key_index(protein_identifiers=genes, sites=sites)
    phospho = pd.DataFrame(
        {
            "A_1": [1.0, 2.0, 0.9],
            "A_2": [1.2, 2.1, 1.1],
            "B_1": [2.0, 1.8, 0.8],
            "B_2": [2.2, 2.0, 1.0],
        },
        index=site_index,
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


def _workflow_request(
    *,
    design: ExperimentalDesign,
    samples: tuple[str, ...] = ("A_1", "A_2", "B_1", "B_2"),
    minimum_condition_replicates: int = 2,
    paired_design_policy: str = "reject",
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
            minimum_condition_replicates=minimum_condition_replicates,
            paired_design_policy=paired_design_policy,
        ),
    )


def test_workflow_differential_validation_rejects_confounded_batch_design() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="B_1", condition="B", batch="batch_2"),
            SampleDesignRecord(sample_id="B_2", condition="B", batch="batch_2"),
        ),
        fixed_effects=(BatchCovariate(),),
    )
    with pytest.raises(
        WorkflowValidationError,
        match="rank deficient.*confounded",
    ):
        DifferentialAnalysisWorkflow().run(_workflow_request(design=design))


def test_workflow_differential_validation_accepts_balanced_batch_fixed_effect() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", batch="batch_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", batch="batch_2"),
            SampleDesignRecord(sample_id="B_1", condition="B", batch="batch_1"),
            SampleDesignRecord(sample_id="B_2", condition="B", batch="batch_2"),
        ),
        fixed_effects=(BatchCovariate(),),
    )

    result = DifferentialAnalysisWorkflow().run(_workflow_request(design=design))

    assert result.policy_provenance is not None
    assert result.policy_provenance.design.coefficient_labels == (
        "A",
        "B",
        "batch[batch_2]",
    )


def test_differential_block_default_policy_rejects_explicit_block_ids() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", block_id="pair_2"),
            SampleDesignRecord(sample_id="B_1", condition="B", block_id="pair_1"),
            SampleDesignRecord(sample_id="B_2", condition="B", block_id="pair_2"),
        )
    )
    with pytest.raises(
        WorkflowValidationError,
        match="paired_design_policy='fixed_block'",
    ):
        DifferentialAnalysisWorkflow().run(_workflow_request(design=design))


def test_differential_block_default_policy_allows_unblocked_designs() -> None:
    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A"),
            SampleDesignRecord(sample_id="A_2", condition="A"),
            SampleDesignRecord(sample_id="B_1", condition="B"),
            SampleDesignRecord(sample_id="B_2", condition="B"),
        )
    )

    result = DifferentialAnalysisWorkflow().run(_workflow_request(design=design))

    assert result.table_for("B_vs_A").shape[0] == 3


def test_differential_block_fixed_block_invalid_design_skips_executor() -> None:
    calls = {"executor": 0}

    class _ExecutorSpy:
        def run(self, request: object) -> object:
            calls["executor"] += 1
            raise AssertionError("executor should not be called")

    design = ExperimentalDesign(
        samples=(
            SampleDesignRecord(sample_id="A_1", condition="A", block_id="pair_1"),
            SampleDesignRecord(sample_id="A_2", condition="A", block_id="pair_2"),
            SampleDesignRecord(sample_id="B_1", condition="B", block_id="pair_1"),
            SampleDesignRecord(sample_id="B_2", condition="B", block_id="pair_3"),
        )
    )
    with pytest.raises(
        WorkflowValidationError,
        match="fixed_block.*at least 2 samples.*incomplete blocks",
    ):
        DifferentialAnalysisWorkflow(
            executor=_ExecutorSpy(),  # type: ignore[arg-type]
        ).run(
            _workflow_request(
                design=design,
                paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
            )
        )
    assert calls["executor"] == 0


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
