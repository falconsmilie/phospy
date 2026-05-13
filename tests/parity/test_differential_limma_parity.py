from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api import (
    Contrast,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    EmpiricalBayesConfig,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)

pytestmark = [pytest.mark.parity]

DIFF_PARITY_DIR = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "rewrite_parity"
    / "differential_limma_envelope"
)
PARITY_RTOL = 1e-6
PARITY_ATOL = 1e-8


def _canonical_site_id(raw_site_id: str, *, ordinal: int) -> str:
    tokens = [token.strip() for token in raw_site_id.split(";") if token.strip()]
    if len(tokens) >= 2:
        return f"{tokens[0]};{tokens[1]};"
    return f"SITE{ordinal};S{ordinal};"


def _load_matrix() -> pd.DataFrame:
    frame = pd.read_csv(DIFF_PARITY_DIR / "matrix.csv")
    matrix = frame.set_index("site_id")
    raw_site_ids = matrix.index.astype(str).tolist()
    canonical_ids = [
        _canonical_site_id(site_id, ordinal=idx)
        for idx, site_id in enumerate(raw_site_ids, start=1)
    ]
    matrix.index = pd.Index(canonical_ids, name=matrix.index.name)
    return matrix


def _load_design() -> pd.DataFrame:
    frame = pd.read_csv(DIFF_PARITY_DIR / "design.csv")
    return frame.set_index("sample")


def _load_contrasts_matrix() -> pd.DataFrame:
    frame = pd.read_csv(DIFF_PARITY_DIR / "contrasts.csv")
    return frame.set_index("coefficient")


def _load_expected(contrast_name: str) -> pd.DataFrame:
    frame = pd.read_csv(DIFF_PARITY_DIR / f"limma_{contrast_name}.csv")
    expected = frame.set_index("site_id").loc[:, ["logFC", "t", "P.Value", "adj.P.Val"]]
    raw_site_ids = expected.index.astype(str).tolist()
    canonical_ids = [
        _canonical_site_id(site_id, ordinal=idx)
        for idx, site_id in enumerate(raw_site_ids, start=1)
    ]
    expected.index = pd.Index(canonical_ids, name=expected.index.name)
    return expected


def _dataset_from_matrix(matrix: pd.DataFrame) -> AnalysisReadyPhosphoDataset:
    parsed = [site_id.split(";") for site_id in matrix.index.astype(str).tolist()]
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": [parts[0] for parts in parsed],
            "site": [parts[1] for parts in parsed],
            "site_sequence": ["A" * 31] * matrix.shape[0],
            "localisation_confidence": [0.95] * matrix.shape[0],
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


def _design_from_matrix(design: pd.DataFrame) -> ExperimentalDesign:
    records: list[SampleDesignRecord] = []
    replicate_counts: dict[str, int] = {}
    for sample_id, row in design.iterrows():
        active_terms = [term for term, value in row.items() if float(value) == 1.0]
        if len(active_terms) != 1:
            raise AssertionError(
                "parity fixture design must be one-hot encoded per sample"
            )
        condition = str(active_terms[0])
        replicate_counts.setdefault(condition, 0)
        replicate_counts[condition] += 1
        records.append(
            SampleDesignRecord(
                sample_id=str(sample_id),
                condition=condition,
                biological_replicate_id=f"{condition}_r{replicate_counts[condition]}",
            )
        )
    return ExperimentalDesign(samples=tuple(records))


def _contrasts_from_matrix(contrasts: pd.DataFrame) -> tuple[Contrast, ...]:
    typed: list[Contrast] = []
    for name in contrasts.columns:
        vector = contrasts.loc[:, name]
        numerator_terms = [
            str(term) for term, value in vector.items() if float(value) == 1.0
        ]
        denominator_terms = [
            str(term) for term, value in vector.items() if float(value) == -1.0
        ]
        if len(numerator_terms) != 1 or len(denominator_terms) != 1:
            raise AssertionError(
                "parity fixture contrasts must have one +1 and one -1 per column"
            )
        typed.append(
            Contrast(
                name=str(name),
                numerator_condition=numerator_terms[0],
                denominator_condition=denominator_terms[0],
            )
        )
    return tuple(typed)


def _run_fixture_workflow() -> object:
    matrix = _load_matrix()
    design_matrix = _load_design()
    contrast_matrix = _load_contrasts_matrix()
    request = DifferentialAnalysisRequest(
        dataset=_dataset_from_matrix(matrix),
        design=_design_from_matrix(design_matrix),
        contrasts=_contrasts_from_matrix(contrast_matrix),
        config=DifferentialAnalysisConfig(
            empirical_bayes=EmpiricalBayesConfig(method="standard")
        ),
    )
    return DifferentialAnalysisWorkflow().run(request)


def test_differential_limma_envelope_fixture_is_source_labelled() -> None:
    required = (
        "matrix.csv",
        "design.csv",
        "contrasts.csv",
        "limma_B_vs_A.csv",
        "limma_A_vs_B.csv",
        "PROVENANCE.md",
    )
    for name in required:
        path = DIFF_PARITY_DIR / name
        assert path.is_file(), f"missing differential parity fixture: {name}"
        assert path.stat().st_size > 0, f"empty differential parity fixture: {name}"

    provenance = (DIFF_PARITY_DIR / "PROVENANCE.md").read_text(encoding="utf-8")
    assert "Generated with R version" in provenance
    assert "limma version:" in provenance
    assert "Design: ~0 + group with groups A/B and 2 replicates each" in provenance
    assert "Contrasts (column order): B_vs_A then A_vs_B" in provenance


def test_two_condition_small_n_parity_matches_limma_with_explicit_tolerance() -> None:
    result = _run_fixture_workflow()

    assert result.residual_degrees_of_freedom == pytest.approx(2.0)
    assert list(result.contrast_tables) == ["B_vs_A", "A_vs_B"]

    for contrast_name in ("B_vs_A", "A_vs_B"):
        observed = result.table_for(contrast_name)
        expected = _load_expected(contrast_name)
        pdt.assert_frame_equal(
            observed,
            expected,
            check_exact=False,
            rtol=PARITY_RTOL,
            atol=PARITY_ATOL,
        )


def test_reverse_contrast_preserves_sign_convention_and_order() -> None:
    result = _run_fixture_workflow()
    b_vs_a = result.table_for("B_vs_A")
    a_vs_b = result.table_for("A_vs_B")

    np.testing.assert_allclose(
        b_vs_a.loc[:, "logFC"].to_numpy(dtype=float),
        -a_vs_b.loc[:, "logFC"].to_numpy(dtype=float),
        rtol=PARITY_RTOL,
        atol=PARITY_ATOL,
    )
    np.testing.assert_allclose(
        b_vs_a.loc[:, "t"].to_numpy(dtype=float),
        -a_vs_b.loc[:, "t"].to_numpy(dtype=float),
        rtol=PARITY_RTOL,
        atol=PARITY_ATOL,
    )
    np.testing.assert_allclose(
        b_vs_a.loc[:, "P.Value"].to_numpy(dtype=float),
        a_vs_b.loc[:, "P.Value"].to_numpy(dtype=float),
        rtol=PARITY_RTOL,
        atol=PARITY_ATOL,
    )
    np.testing.assert_allclose(
        b_vs_a.loc[:, "adj.P.Val"].to_numpy(dtype=float),
        a_vs_b.loc[:, "adj.P.Val"].to_numpy(dtype=float),
        rtol=PARITY_RTOL,
        atol=PARITY_ATOL,
    )


def test_zero_variance_feature_matches_limma_and_remains_finite() -> None:
    observed = _run_fixture_workflow().table_for("B_vs_A")
    expected = _load_expected("B_vs_A")
    site_id = "SITE5;S5;"

    pdt.assert_series_equal(
        observed.loc[site_id],
        expected.loc[site_id],
        check_exact=False,
        rtol=PARITY_RTOL,
        atol=PARITY_ATOL,
    )
    assert observed.at[site_id, "logFC"] == pytest.approx(0.0)
    assert observed.at[site_id, "t"] == pytest.approx(0.0)
    assert observed.at[site_id, "P.Value"] == pytest.approx(1.0)
    assert observed.at[site_id, "adj.P.Val"] == pytest.approx(1.0)
