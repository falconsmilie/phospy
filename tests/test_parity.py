from __future__ import annotations

from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from phosrpy import KinaseActivityAnalyzer, PhosphoDataset

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_DATA = ROOT / "examples" / "data"
R_FIXTURES = ROOT / "tests" / "fixtures" / "r_reference"

CORE_FIXTURE_FILES = [
    "df_total_unique.csv",
    "df_total_filtered.csv",
    "df_phospho_filtered.csv",
    "df_phospho_corrected.csv",
    "phosr_input.csv",
    "mat_phospho_corrected.csv",
    "site_sequences.csv",
]

KINASE_FIXTURE_FILES = [
    "predMat.csv",
    "kinase_activity_matrix.csv",
    "ksea_scores.csv",
    "ksea_counts.csv",
    "kinase_target_counts.csv",
]


def _require_fixture_files(names: list[str]) -> None:
    missing = [name for name in names if not (R_FIXTURES / name).exists()]
    if missing:
        pytest.skip(
            "R reference fixtures are not present yet. "
            "Generate them with `Rscript scripts/generate_r_fixtures.R` before running parity tests. "
            f"Missing: {', '.join(missing)}"
        )


def _read_table(name: str) -> pd.DataFrame:
    return pd.read_csv(R_FIXTURES / name)


def _read_indexed_table(name: str) -> pd.DataFrame:
    return pd.read_csv(R_FIXTURES / name, index_col=0)


def _read_sequences(name: str = "site_sequences.csv") -> pd.Series:
    frame = pd.read_csv(R_FIXTURES / name)
    if {"site_id", "centralized_sequence"} <= set(frame.columns):
        series = frame.set_index("site_id")["centralized_sequence"]
    else:
        series = frame.set_index(frame.columns[0]).iloc[:, 0]
        series.index.name = "site_id"
        series.name = "centralized_sequence"
    return series.sort_index()


def _sort_table(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    present = [col for col in columns if col in df.columns]
    return df.sort_values(present).reset_index(drop=True) if present else df.reset_index(drop=True)


def _normalize_numeric_frame(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in out.columns:
        out[col] = pd.to_numeric(out[col], errors="ignore")
    return out


def _assert_frame_close(actual: pd.DataFrame, expected: pd.DataFrame, *, atol: float = 1e-8, rtol: float = 1e-6) -> None:
    actual = _normalize_numeric_frame(actual)
    expected = _normalize_numeric_frame(expected)
    pdt.assert_frame_equal(actual, expected, check_dtype=False, atol=atol, rtol=rtol)


def _assert_series_close(actual: pd.Series, expected: pd.Series, *, atol: float = 1e-8, rtol: float = 1e-6) -> None:
    actual = pd.to_numeric(actual, errors="ignore")
    expected = pd.to_numeric(expected, errors="ignore")
    pdt.assert_series_equal(actual, expected, check_dtype=False, atol=atol, rtol=rtol)


@pytest.mark.parity
def test_core_outputs_match_r_reference() -> None:
    _require_fixture_files(CORE_FIXTURE_FILES)

    dataset = PhosphoDataset.from_files(
        total_path=EXAMPLE_DATA / "total.tsv",
        phospho_path=EXAMPLE_DATA / "phospho.tsv",
    )
    result = dataset.process_core()

    actual_total_unique = _sort_table(result.total_unique, ["genes"])
    expected_total_unique = _sort_table(_read_table("df_total_unique.csv"), ["genes"])
    _assert_frame_close(actual_total_unique, expected_total_unique)

    actual_total_filtered = _sort_table(result.total_filtered, ["genes"])
    expected_total_filtered = _sort_table(_read_table("df_total_filtered.csv"), ["genes"])
    _assert_frame_close(actual_total_filtered, expected_total_filtered)

    actual_phospho_filtered = _sort_table(result.phospho_filtered, ["gene_p_site", "uid"])
    expected_phospho_filtered = _sort_table(_read_table("df_phospho_filtered.csv"), ["gene_p_site", "uid"])
    _assert_frame_close(actual_phospho_filtered, expected_phospho_filtered)

    actual_phospho_corrected = _sort_table(result.phospho_corrected, ["gene_p_site", "uid"])
    expected_phospho_corrected = _sort_table(_read_table("df_phospho_corrected.csv"), ["gene_p_site", "uid"])
    _assert_frame_close(actual_phospho_corrected, expected_phospho_corrected)

    actual_phosr_input = _sort_table(result.site_matrix.phosr_input, ["site_id", "uid"])
    expected_phosr_input = _sort_table(_read_table("phosr_input.csv"), ["site_id", "uid"])
    _assert_frame_close(actual_phosr_input, expected_phosr_input)

    actual_matrix = result.site_matrix.matrix.sort_index().sort_index(axis=1)
    expected_matrix = _read_indexed_table("mat_phospho_corrected.csv").sort_index().sort_index(axis=1)
    _assert_frame_close(actual_matrix, expected_matrix)

    actual_sequences = result.site_matrix.sequences.sort_index()
    expected_sequences = _read_sequences().sort_index()
    _assert_series_close(actual_sequences, expected_sequences)


@pytest.mark.parity
def test_kinase_outputs_match_r_reference() -> None:
    _require_fixture_files(CORE_FIXTURE_FILES + KINASE_FIXTURE_FILES)

    dataset = PhosphoDataset.from_files(
        total_path=EXAMPLE_DATA / "total.tsv",
        phospho_path=EXAMPLE_DATA / "phospho.tsv",
    )
    core = dataset.process_core()

    pred_mat = _read_indexed_table("predMat.csv")
    analyzer = KinaseActivityAnalyzer(pred_mat=pred_mat)
    result = analyzer.analyze(core.site_matrix.matrix)

    actual_weighted = result.weighted_activity.sort_index().sort_index(axis=1)
    expected_weighted = _read_indexed_table("kinase_activity_matrix.csv").sort_index().sort_index(axis=1)
    _assert_frame_close(actual_weighted, expected_weighted)

    actual_ksea = result.ksea_scores.sort_index().sort_index(axis=1)
    expected_ksea = _read_indexed_table("ksea_scores.csv").sort_index().sort_index(axis=1)
    _assert_frame_close(actual_ksea, expected_ksea)

    actual_ksea_counts = result.ksea_counts.sort_index()
    expected_ksea_counts = _read_indexed_table("ksea_counts.csv").iloc[:, 0].sort_index()
    expected_ksea_counts.name = actual_ksea_counts.name
    _assert_series_close(actual_ksea_counts, expected_ksea_counts)

    actual_target_counts = result.target_counts.sort_index()
    expected_target_counts = _read_indexed_table("kinase_target_counts.csv").iloc[:, 0].sort_index()
    expected_target_counts.name = actual_target_counts.name
    _assert_series_close(actual_target_counts, expected_target_counts)
