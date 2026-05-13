from __future__ import annotations

from pathlib import Path

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
    / "differential_r_reference"
)
PARITY_RTOL = 1e-6
PARITY_ATOL = 1e-8


def _canonical_site_id(raw_site_id: str, *, ordinal: int) -> str:
    tokens = [token.strip() for token in raw_site_id.split(";") if token.strip()]
    if len(tokens) >= 2:
        return f"{tokens[0]};{tokens[1]};"
    return f"SITE_{ordinal};S{ordinal};"


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


def _load_expected(contrast_name: str, *, matrix_site_index: pd.Index) -> pd.DataFrame:
    frame = pd.read_csv(DIFF_PARITY_DIR / f"limma_{contrast_name}.csv")
    frame = frame.set_index("site_id")
    raw_site_ids = frame.index.astype(str).tolist()
    canonical_ids = [
        _canonical_site_id(site_id, ordinal=idx)
        for idx, site_id in enumerate(raw_site_ids, start=1)
    ]
    expected = frame.loc[:, ["logFC", "t", "P.Value", "adj.P.Val"]].copy(deep=True)
    expected.index = pd.Index(canonical_ids, name=expected.index.name)
    expected = expected.reindex(matrix_site_index)
    return expected


def _dataset_from_matrix(matrix: pd.DataFrame) -> AnalysisReadyPhosphoDataset:
    site_ids = matrix.index.astype(str).tolist()
    parsed = [
        _canonical_site_id(site_id, ordinal=idx).split(";")
        for idx, site_id in enumerate(site_ids, start=1)
    ]
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": [parts[0] for parts in parsed],
            "site": [parts[1] for parts in parsed],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in [parts[1] for parts in parsed]
            ],
            "localisation_confidence": [0.95] * len(parsed),
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


def test_differential_parity_fixtures_are_present_and_readable() -> None:
    required = (
        "matrix.csv",
        "design.csv",
        "contrasts.csv",
        "limma_B_vs_A.csv",
        "limma_C_vs_A.csv",
        "PROVENANCE.md",
    )
    for name in required:
        path = DIFF_PARITY_DIR / name
        assert path.is_file(), f"missing differential parity fixture: {name}"
        assert path.stat().st_size > 0, f"empty differential parity fixture: {name}"

    provenance = (DIFF_PARITY_DIR / "PROVENANCE.md").read_text(encoding="utf-8")
    assert "limma version:" in provenance
    assert "Design: ~0 + group" in provenance


def test_differential_analysis_matches_limma_parity_within_tolerance() -> None:
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
    result = DifferentialAnalysisWorkflow().run(request)

    for contrast_name in ("B_vs_A", "C_vs_A"):
        observed = result.table_for(contrast_name).sort_index()
        expected = _load_expected(
            contrast_name,
            matrix_site_index=matrix.index,
        ).sort_index()
        pdt.assert_frame_equal(
            observed,
            expected,
            check_exact=False,
            rtol=PARITY_RTOL,
            atol=PARITY_ATOL,
        )
