from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    EmpiricalBayesConfig,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.errors import DatasetValidationError
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)

FIXTURE_DIR = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "rewrite_parity"
    / "differential_limma_envelope"
)
NEGATIVE_FIXTURE_DIR = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "rewrite_parity"
    / "differential_contract_negative_cases"
)


def _canonical_site_id(raw_site_id: str, *, ordinal: int) -> str:
    tokens = [token.strip() for token in raw_site_id.split(";") if token.strip()]
    if len(tokens) >= 2:
        return f"{tokens[0]};{tokens[1]};"
    return f"SITE{ordinal};S{ordinal};"


def _build_dataset(matrix: pd.DataFrame):
    from phospy import AnalysisReadyPhosphoDataset

    parsed = [site_id.split(";") for site_id in matrix.index.astype(str).tolist()]
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": [parts[0] for parts in parsed],
            "site": [parts[1] for parts in parsed],
            "site_sequence": ["A" * 31] * matrix.shape[0],
            "protein_id": [parts[0] for parts in parsed],
        },
        index=matrix.index.copy(),
    )
    return AnalysisReadyPhosphoDataset(
        phospho=matrix,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _load_matrix() -> pd.DataFrame:
    frame = pd.read_csv(FIXTURE_DIR / "matrix.csv")
    matrix = frame.set_index("site_id")
    raw_site_ids = matrix.index.astype(str).tolist()
    canonical_ids = [
        _canonical_site_id(site_id, ordinal=idx)
        for idx, site_id in enumerate(raw_site_ids, start=1)
    ]
    matrix.index = pd.Index(canonical_ids, name=matrix.index.name)
    return matrix


def _request_for_reverse_contrasts(
    matrix: pd.DataFrame,
) -> DifferentialAnalysisRequest:
    design = ExperimentalDesign(
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
    contrasts = (
        Contrast(
            name="B_vs_A",
            numerator_condition="B",
            denominator_condition="A",
        ),
        Contrast(
            name="A_vs_B",
            numerator_condition="A",
            denominator_condition="B",
        ),
    )
    return DifferentialAnalysisRequest(
        dataset=_build_dataset(matrix),
        design=design,
        contrasts=contrasts,
        empirical_bayes=EmpiricalBayesConfig(method="standard"),
    )


def test_result_tables_follow_public_differential_contract() -> None:
    result = DifferentialAnalysisWorkflow().run(
        _request_for_reverse_contrasts(_load_matrix())
    )

    assert list(result.contrast_tables) == ["B_vs_A", "A_vs_B"]
    for contrast_name in ("B_vs_A", "A_vs_B"):
        table = result.table_for(contrast_name)
        assert list(table.columns) == ["logFC", "t", "P.Value", "adj.P.Val"]
        assert np.isfinite(table.loc[:, "logFC"]).all()
        assert np.isfinite(table.loc[:, "t"]).all()
        assert (table.loc[:, "P.Value"] >= 0.0).all()
        assert (table.loc[:, "P.Value"] <= 1.0).all()
        assert (table.loc[:, "adj.P.Val"] >= 0.0).all()
        assert (table.loc[:, "adj.P.Val"] <= 1.0).all()


def test_reverse_contrasts_are_directionally_consistent() -> None:
    result = DifferentialAnalysisWorkflow().run(
        _request_for_reverse_contrasts(_load_matrix())
    )
    b_vs_a = result.table_for("B_vs_A")
    a_vs_b = result.table_for("A_vs_B")

    np.testing.assert_allclose(
        b_vs_a.loc[:, "logFC"].to_numpy(dtype=float),
        -a_vs_b.loc[:, "logFC"].to_numpy(dtype=float),
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        b_vs_a.loc[:, "t"].to_numpy(dtype=float),
        -a_vs_b.loc[:, "t"].to_numpy(dtype=float),
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        b_vs_a.loc[:, "P.Value"].to_numpy(dtype=float),
        a_vs_b.loc[:, "P.Value"].to_numpy(dtype=float),
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        b_vs_a.loc[:, "adj.P.Val"].to_numpy(dtype=float),
        a_vs_b.loc[:, "adj.P.Val"].to_numpy(dtype=float),
        rtol=1e-12,
        atol=1e-12,
    )


def test_missing_values_are_rejected_before_differential_execution() -> None:
    matrix = pd.read_csv(NEGATIVE_FIXTURE_DIR / "matrix_with_missing.csv").set_index(
        "site_id"
    )
    raw_site_ids = matrix.index.astype(str).tolist()
    canonical_ids = [
        _canonical_site_id(site_id, ordinal=idx)
        for idx, site_id in enumerate(raw_site_ids, start=1)
    ]
    matrix.index = pd.Index(canonical_ids, name=matrix.index.name)

    with pytest.raises(DatasetValidationError, match="dataset.phospho"):
        _build_dataset(matrix)
