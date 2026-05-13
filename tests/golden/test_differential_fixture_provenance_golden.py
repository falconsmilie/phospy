from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
PARITY_FIXTURE_ROOT = ROOT / "tests" / "fixtures" / "rewrite_parity"

pytestmark = [pytest.mark.golden, pytest.mark.reproducibility, pytest.mark.release_gate]


def _require_common_provenance_fields(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    assert "Generated with R version" in text
    assert "limma version:" in text
    assert "Design:" in text
    assert "Contrasts" in text


def _require_expected_schema(path: Path) -> None:
    frame = pd.read_csv(path)
    assert list(frame.columns) == ["site_id", "logFC", "t", "P.Value", "adj.P.Val"]


def test_differential_limma_fixture_provenance_is_source_labelled() -> None:
    base_fixture = PARITY_FIXTURE_ROOT / "differential_r_reference" / "PROVENANCE.md"
    envelope_fixture = (
        PARITY_FIXTURE_ROOT / "differential_limma_envelope" / "PROVENANCE.md"
    )

    assert base_fixture.is_file()
    assert envelope_fixture.is_file()
    _require_common_provenance_fields(base_fixture)
    _require_common_provenance_fields(envelope_fixture)


def test_differential_limma_expected_tables_keep_stable_schema() -> None:
    _require_expected_schema(
        PARITY_FIXTURE_ROOT / "differential_r_reference" / "limma_B_vs_A.csv"
    )
    _require_expected_schema(
        PARITY_FIXTURE_ROOT / "differential_r_reference" / "limma_C_vs_A.csv"
    )
    _require_expected_schema(
        PARITY_FIXTURE_ROOT / "differential_limma_envelope" / "limma_B_vs_A.csv"
    )
    _require_expected_schema(
        PARITY_FIXTURE_ROOT / "differential_limma_envelope" / "limma_A_vs_B.csv"
    )
