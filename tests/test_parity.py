from __future__ import annotations

from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import KinaseActivityAnalyzer, PhosphoDataset
from phospy.internal.constants import (
    CENTRALIZED_SEQUENCE_COLUMN,
    CORE_OUTPUT_ARTIFACT_BASENAMES,
    CORE_SITE_SEQUENCES_BASENAME,
    KINASE_ACTIVITY_MATRIX_FILENAME,
    KINASE_OUTPUT_FILENAMES,
    KINASE_TARGET_COUNTS_FILENAME,
    KSEA_COUNTS_FILENAME,
    KSEA_SCORES_FILENAME,
    SITE_MATRIX_ID_COLUMN,
)
from phospy.preprocessing import CorePreprocessingConfig

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_DATA = ROOT / "examples" / "data"
R_FIXTURES = ROOT / "tests" / "fixtures" / "r_reference"
R_FIXTURES_L6 = ROOT / "tests" / "fixtures" / "r_reference_l6"

CORE_FIXTURE_FILES = [f"{basename}.csv" for basename in CORE_OUTPUT_ARTIFACT_BASENAMES]

KINASE_FIXTURE_FILES = ["predMat.csv", *KINASE_OUTPUT_FILENAMES[:-1]]


EXAMPLE_COMPARISONS = [
    ("group1", "group4"),
    ("group2", "group5"),
    ("group3", "group6"),
    ("group1", "group2"),
    ("group1", "group3"),
    ("group2", "group3"),
    ("group4", "group5"),
    ("group4", "group6"),
    ("group5", "group6"),
]

PHOSPHO_FIXTURE_ENCODING = "utf-16le"

L6_FIXTURE_FILES = [
    "l6_phospho_matrix.csv",
    "predMat.csv",
    "kinase_activity_matrix.csv",
    "ksea_scores.csv",
    "ksea_counts.csv",
    "kinase_target_counts.csv",
]


def _require_fixture_files(names: list[str], fixture_dir: Path = R_FIXTURES) -> None:
    missing = [name for name in names if not (fixture_dir / name).exists()]
    if missing:
        pytest.skip(
            "R reference fixtures are not present yet. "
            "Generate the relevant fixtures under `tests/fixtures/` before running parity tests. "
            f"Missing: {', '.join(missing)}"
        )


def _read_table(name: str, fixture_dir: Path = R_FIXTURES) -> pd.DataFrame:
    return pd.read_csv(fixture_dir / name)


def _read_indexed_table(name: str, fixture_dir: Path = R_FIXTURES) -> pd.DataFrame:
    return pd.read_csv(fixture_dir / name, index_col=0)


def _read_sequences(
    name: str = f"{CORE_SITE_SEQUENCES_BASENAME}.csv", fixture_dir: Path = R_FIXTURES
) -> pd.Series:
    frame = pd.read_csv(fixture_dir / name)
    if {SITE_MATRIX_ID_COLUMN, CENTRALIZED_SEQUENCE_COLUMN} <= set(frame.columns):
        series = frame.set_index(SITE_MATRIX_ID_COLUMN)[CENTRALIZED_SEQUENCE_COLUMN]
    else:
        series = frame.set_index(frame.columns[0]).iloc[:, 0]
        series.index.name = SITE_MATRIX_ID_COLUMN
        series.name = CENTRALIZED_SEQUENCE_COLUMN
    return series.sort_index()


def _sort_table(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    present = [col for col in columns if col in df.columns]
    return (
        df.sort_values(present).reset_index(drop=True)
        if present
        else df.reset_index(drop=True)
    )


def _normalize_numeric_frame(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    for col in out.columns:
        try:
            out[col] = pd.to_numeric(out[col], errors="raise")
        except (ValueError, TypeError):
            # Leave genuinely non-numeric columns, like gene names, unchanged.
            pass

    return out


def _assert_frame_close(
    actual: pd.DataFrame,
    expected: pd.DataFrame,
    *,
    atol: float = 1e-8,
    rtol: float = 1e-6,
) -> None:
    actual = _normalize_numeric_frame(actual)
    expected = _normalize_numeric_frame(expected)

    if actual.empty and expected.empty:
        if list(actual.columns) != list(expected.columns):
            assert len(actual.columns) == len(expected.columns)
        return

    pdt.assert_frame_equal(
        actual,
        expected,
        check_dtype=False,
        check_index_type=False,
        check_column_type=False,
        check_names=False,
        atol=atol,
        rtol=rtol,
    )


def _normalize_series(series: pd.Series) -> pd.Series:
    try:
        return pd.to_numeric(series, errors="raise")
    except (ValueError, TypeError):
        return series


def _assert_series_close(
    actual: pd.Series, expected: pd.Series, *, atol: float = 1e-8, rtol: float = 1e-6
) -> None:
    actual = _normalize_series(actual)
    expected = _normalize_series(expected)
    pdt.assert_series_equal(
        actual,
        expected,
        check_dtype=False,
        check_index_type=False,
        check_names=False,
        atol=atol,
        rtol=rtol,
    )


@pytest.mark.parity
def test_core_outputs_match_r_reference() -> None:
    _require_fixture_files(CORE_FIXTURE_FILES)

    dataset = PhosphoDataset.from_files(
        total_path=EXAMPLE_DATA / "total.tsv",
        phospho_path=EXAMPLE_DATA / "phospho.tsv",
        comparisons=EXAMPLE_COMPARISONS,
        phospho_encoding=PHOSPHO_FIXTURE_ENCODING,
    )
    result = dataset.preprocessing.run(config=CorePreprocessingConfig())

    actual_total_unique = _sort_table(result.total_unique, ["genes"])
    expected_total_unique = _sort_table(_read_table("df_total_unique.csv"), ["genes"])
    _assert_frame_close(actual_total_unique, expected_total_unique)

    actual_total_filtered = _sort_table(result.total_filtered, ["genes"])
    expected_total_filtered = _sort_table(
        _read_table("df_total_filtered.csv"), ["genes"]
    )
    _assert_frame_close(actual_total_filtered, expected_total_filtered)

    actual_phospho_filtered = _sort_table(
        result.phospho_filtered, ["gene_p_site", "uid"]
    )
    expected_phospho_filtered = _sort_table(
        _read_table("df_phospho_filtered.csv"), ["gene_p_site", "uid"]
    )
    _assert_frame_close(actual_phospho_filtered, expected_phospho_filtered)

    actual_phospho_corrected = _sort_table(
        result.phospho_corrected, ["gene_p_site", "uid"]
    )
    expected_phospho_corrected = _sort_table(
        _read_table("df_phospho_corrected.csv"), ["gene_p_site", "uid"]
    )
    _assert_frame_close(actual_phospho_corrected, expected_phospho_corrected)

    actual_phosr_input = _sort_table(result.site_matrix.phosr_input, ["site_id", "uid"])
    expected_phosr_input = _sort_table(
        _read_table("phosr_input.csv"), ["site_id", "uid"]
    )
    _assert_frame_close(actual_phosr_input, expected_phosr_input)

    actual_matrix = result.site_matrix.matrix.sort_index().sort_index(axis=1)
    expected_matrix = (
        _read_indexed_table("mat_phospho_corrected.csv").sort_index().sort_index(axis=1)
    )
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
        comparisons=EXAMPLE_COMPARISONS,
        phospho_encoding=PHOSPHO_FIXTURE_ENCODING,
    )
    core = dataset.preprocessing.run(config=CorePreprocessingConfig())

    pred_mat = _read_indexed_table("predMat.csv")
    result = KinaseActivityAnalyzer().run(
        pred_mat=pred_mat, phospho_matrix=core.site_matrix.matrix
    )

    actual_weighted = result.weighted_activity.sort_index().sort_index(axis=1)
    expected_weighted = (
        _read_indexed_table(KINASE_ACTIVITY_MATRIX_FILENAME)
        .sort_index()
        .sort_index(axis=1)
    )
    _assert_frame_close(actual_weighted, expected_weighted)

    actual_ksea = result.ksea_scores.sort_index().sort_index(axis=1)
    expected_ksea = (
        _read_indexed_table(KSEA_SCORES_FILENAME).sort_index().sort_index(axis=1)
    )
    _assert_frame_close(actual_ksea, expected_ksea)

    actual_ksea_counts = result.ksea_counts.sort_index()
    expected_ksea_counts_frame = _read_indexed_table(KSEA_COUNTS_FILENAME)
    expected_ksea_counts = (
        expected_ksea_counts_frame.iloc[:, 0].sort_index()
        if not expected_ksea_counts_frame.empty
        else pd.Series(dtype=int)
    )
    expected_ksea_counts.name = actual_ksea_counts.name
    _assert_series_close(actual_ksea_counts, expected_ksea_counts)

    actual_target_counts = result.target_counts.sort_index()
    expected_target_counts = (
        _read_indexed_table(KINASE_TARGET_COUNTS_FILENAME).iloc[:, 0].sort_index()
    )
    expected_target_counts.name = actual_target_counts.name
    _assert_series_close(actual_target_counts, expected_target_counts)


@pytest.mark.parity
def test_l6_kinase_outputs_match_r_reference() -> None:
    _require_fixture_files(L6_FIXTURE_FILES, fixture_dir=R_FIXTURES_L6)

    phospho_matrix = (
        _read_indexed_table("l6_phospho_matrix.csv", fixture_dir=R_FIXTURES_L6)
        .sort_index()
        .sort_index(axis=1)
    )
    pred_mat = _read_indexed_table("predMat.csv", fixture_dir=R_FIXTURES_L6)

    result = KinaseActivityAnalyzer().run(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=0.6,
        min_substrates=3,
        top_n_substrates=20,
    )

    actual_weighted = result.weighted_activity.sort_index().sort_index(axis=1)
    expected_weighted = (
        _read_indexed_table(KINASE_ACTIVITY_MATRIX_FILENAME, fixture_dir=R_FIXTURES_L6)
        .sort_index()
        .sort_index(axis=1)
    )
    _assert_frame_close(actual_weighted, expected_weighted)

    actual_ksea = result.ksea_scores.sort_index().sort_index(axis=1)
    expected_ksea = (
        _read_indexed_table(KSEA_SCORES_FILENAME, fixture_dir=R_FIXTURES_L6)
        .sort_index()
        .sort_index(axis=1)
    )
    _assert_frame_close(actual_ksea, expected_ksea)

    actual_ksea_counts = result.ksea_counts.sort_index()
    expected_ksea_counts_frame = _read_indexed_table(
        KSEA_COUNTS_FILENAME, fixture_dir=R_FIXTURES_L6
    )
    expected_ksea_counts = (
        expected_ksea_counts_frame.iloc[:, 0].sort_index()
        if not expected_ksea_counts_frame.empty
        else pd.Series(dtype=int)
    )
    expected_ksea_counts.name = actual_ksea_counts.name
    _assert_series_close(actual_ksea_counts, expected_ksea_counts)

    actual_target_counts = result.target_counts.sort_index()
    expected_target_counts = (
        _read_indexed_table(KINASE_TARGET_COUNTS_FILENAME, fixture_dir=R_FIXTURES_L6)
        .iloc[:, 0]
        .sort_index()
    )
    expected_target_counts.name = actual_target_counts.name
    _assert_series_close(actual_target_counts, expected_target_counts)
