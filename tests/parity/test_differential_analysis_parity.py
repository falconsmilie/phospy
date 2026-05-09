from __future__ import annotations

from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.api import (
    ContrastMatrix,
    DesignMatrix,
    DifferentialAnalysis,
    DifferentialAnalysisRequest,
    EmpiricalBayesConfig,
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


def _load_matrix() -> pd.DataFrame:
    frame = pd.read_csv(DIFF_PARITY_DIR / "matrix.csv")
    return frame.set_index("site_id")


def _load_design() -> pd.DataFrame:
    frame = pd.read_csv(DIFF_PARITY_DIR / "design.csv")
    return frame.set_index("sample")


def _load_contrasts() -> pd.DataFrame:
    frame = pd.read_csv(DIFF_PARITY_DIR / "contrasts.csv")
    return frame.set_index("coefficient")


def _load_expected(contrast_name: str) -> pd.DataFrame:
    frame = pd.read_csv(DIFF_PARITY_DIR / f"limma_{contrast_name}.csv")
    frame = frame.set_index("site_id")
    return frame.loc[:, ["logFC", "t", "P.Value", "adj.P.Val"]]


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
    request = DifferentialAnalysisRequest(
        matrix=_load_matrix(),
        design=DesignMatrix(_load_design()),
        contrasts=ContrastMatrix(_load_contrasts()),
        empirical_bayes=EmpiricalBayesConfig(method="standard"),
    )
    result = DifferentialAnalysis().run(request)

    for contrast_name in ("B_vs_A", "C_vs_A"):
        observed = result.table_for(contrast_name).sort_index()
        expected = _load_expected(contrast_name).sort_index()
        pdt.assert_frame_equal(
            observed,
            expected,
            check_exact=False,
            rtol=PARITY_RTOL,
            atol=PARITY_ATOL,
        )
