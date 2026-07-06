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
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
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


def _load_expected(contrast_name: str, *, result_site_index: pd.Index) -> pd.DataFrame:
    frame = pd.read_csv(DIFF_PARITY_DIR / f"limma_{contrast_name}.csv")
    expected = frame.set_index("site_id").loc[:, ["logFC", "t", "P.Value", "adj.P.Val"]]
    raw_site_ids = expected.index.astype(str).tolist()
    canonical_ids = [
        _canonical_site_id(site_id, ordinal=idx)
        for idx, site_id in enumerate(raw_site_ids, start=1)
    ]
    expected.index = site_key_index_from_display_ids(
        canonical_ids,
        protein_namespace="gene_symbol",
    )
    return expected.reindex(result_site_index)


def _dataset_from_matrix(matrix: pd.DataFrame) -> AnalysisReadyPhosphoDataset:
    display_ids = matrix.index.astype(str).tolist()
    site_index = site_key_index_from_display_ids(
        display_ids,
        protein_namespace="gene_symbol",
    )
    parsed = [site_id.split(";") for site_id in display_ids]
    phospho = matrix.copy(deep=True)
    phospho.index = site_index
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": [parts[0] for parts in parsed],
            "site": [parts[1] for parts in parsed],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in [parts[1] for parts in parsed]
            ],
            "localisation_confidence": [0.95] * matrix.shape[0],
            "protein_id": [parts[0] for parts in parsed],
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


def _all_constant_row_mask(matrix: pd.DataFrame) -> pd.Series:
    values = matrix.to_numpy(dtype=float)
    mask = (values == values[:, [0]]).all(axis=1)
    return pd.Series(mask, index=matrix.index.copy())


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


def test_differential_parity_fixture_records_withheld_feature_count() -> None:
    matrix = _load_matrix()
    withheld_mask = _all_constant_row_mask(matrix)
    provenance = (DIFF_PARITY_DIR / "PROVENANCE.md").read_text(encoding="utf-8")

    assert int(withheld_mask.sum()) == 1
    assert matrix.index[withheld_mask].tolist() == ["SITE5;S5;"]
    assert "PhosPy withheld feature count: 1 all-constant row (SITE_05)" in provenance


def test_full_limma_envelope_fixture_withholds_all_constant_row() -> None:
    result = _run_fixture_workflow()
    table = result.table_for("B_vs_A")
    all_constant_rows = (
        table[DIFFERENTIAL_RESULT_STATUS_COLUMN]
        == DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT
    )

    assert all_constant_rows.any()
    assert (
        table.loc[all_constant_rows, ["logFC", "t", "P.Value", "adj.P.Val"]]
        .isna()
        .all()
        .all()
    )


def test_differential_parity_excludes_withheld_rows_from_comparison() -> None:
    result = _run_fixture_workflow()
    table = result.table_for("B_vs_A")
    tested_rows = (
        table[DIFFERENTIAL_RESULT_STATUS_COLUMN] == DIFFERENTIAL_RESULT_STATUS_TESTED
    )
    withheld_rows = (
        table[DIFFERENTIAL_RESULT_STATUS_COLUMN]
        == DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT
    )
    observed = table.loc[tested_rows, ["logFC"]]
    expected = _load_expected(
        "B_vs_A",
        result_site_index=table.index,
    ).loc[tested_rows, ["logFC"]]

    assert int(withheld_rows.sum()) == 1
    assert not observed.index.intersection(table.index[withheld_rows]).any()
    pdt.assert_frame_equal(
        observed,
        expected,
        check_exact=False,
        rtol=PARITY_RTOL,
        atol=PARITY_ATOL,
    )


def test_two_condition_small_n_coefficients_match_limma_with_explicit_tolerance() -> (
    None
):
    result = _run_fixture_workflow()

    assert result.residual_degrees_of_freedom == pytest.approx(2.0)
    assert list(result.contrast_tables) == ["B_vs_A", "A_vs_B"]

    for contrast_name in ("B_vs_A", "A_vs_B"):
        table = result.table_for(contrast_name)
        tested_rows = (
            table[DIFFERENTIAL_RESULT_STATUS_COLUMN]
            == DIFFERENTIAL_RESULT_STATUS_TESTED
        )
        observed = table.loc[tested_rows, ["logFC"]]
        expected = _load_expected(
            contrast_name,
            result_site_index=table.index,
        ).loc[tested_rows, ["logFC"]]
        pdt.assert_frame_equal(
            observed,
            expected,
            check_exact=False,
            rtol=PARITY_RTOL,
            atol=PARITY_ATOL,
        )
        assert np.isfinite(
            table.loc[tested_rows, ["t", "P.Value", "adj.P.Val"]].to_numpy()
        ).all()


def test_reverse_contrast_preserves_sign_convention_and_order() -> None:
    result = _run_fixture_workflow()
    b_vs_a_table = result.table_for("B_vs_A")
    a_vs_b_table = result.table_for("A_vs_B")
    tested_rows = (
        b_vs_a_table[DIFFERENTIAL_RESULT_STATUS_COLUMN]
        == DIFFERENTIAL_RESULT_STATUS_TESTED
    )
    assert tested_rows.equals(
        a_vs_b_table[DIFFERENTIAL_RESULT_STATUS_COLUMN]
        == DIFFERENTIAL_RESULT_STATUS_TESTED
    )
    b_vs_a = b_vs_a_table.loc[tested_rows, ["logFC", "t", "P.Value"]]
    a_vs_b = a_vs_b_table.loc[tested_rows, ["logFC", "t", "P.Value"]]

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
