from __future__ import annotations

from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import KinaseActivityAnalyzer, PhosphoDataset
from phospy.prediction import KinasePredictor, build_candidate_substrate_list
from phospy.profiles import build_kinase_substrate_profiles
from phospy.scoring import KinaseScorer, combine_profile_and_motif_scores

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_DATA = ROOT / "examples" / "data"
R_FIXTURES = ROOT / "tests" / "fixtures" / "r_reference"
R_FIXTURES_L6 = ROOT / "tests" / "fixtures" / "r_reference_l6"

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

L6_FIXTURE_FILES = [
    "l6_phospho_matrix.csv",
    "predMat.csv",
    "kinase_activity_matrix.csv",
    "ksea_scores.csv",
    "ksea_counts.csv",
    "kinase_target_counts.csv",
]

L6_NATIVE_FIXTURE_FILES = [
    "native_substrate_map.csv",
    "native_profile_matrix.csv",
    "native_profile_scores.csv",
    "native_motif_scores.csv",
    "native_motif_sizes.csv",
    "native_combined_scores.csv",
    "native_combined_weights.csv",
    "native_candidate_substrates.csv",
    "native_prediction_top30.csv",
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
    name: str = "site_sequences.csv", fixture_dir: Path = R_FIXTURES
) -> pd.Series:
    frame = pd.read_csv(fixture_dir / name)
    if {"site_id", "centralized_sequence"} <= set(frame.columns):
        series = frame.set_index("site_id")["centralized_sequence"]
    else:
        series = frame.set_index(frame.columns[0]).iloc[:, 0]
        series.index.name = "site_id"
        series.name = "centralized_sequence"
    return series.sort_index()


def _read_named_series(
    name: str,
    *,
    index_col: str,
    value_col: str,
    fixture_dir: Path = R_FIXTURES,
) -> pd.Series:
    frame = pd.read_csv(fixture_dir / name)
    series = frame.set_index(index_col)[value_col]
    series.index.name = index_col
    series.name = value_col
    return series


def _read_grouped_mapping(
    name: str,
    *,
    key_col: str = "kinase",
    value_col: str = "site_id",
    fixture_dir: Path = R_FIXTURES,
) -> dict[str, list[str]]:
    frame = pd.read_csv(fixture_dir / name)
    if frame.empty:
        return {}
    return {
        str(key): [str(value) for value in group[value_col].tolist()]
        for key, group in frame.groupby(key_col, sort=False)
    }


def _read_ranked_sites_by_kinase(
    name: str, *, fixture_dir: Path = R_FIXTURES
) -> dict[str, list[str]]:
    frame = pd.read_csv(fixture_dir / name)
    if frame.empty:
        return {}
    frame = frame.sort_values(["kinase", "rank"]).reset_index(drop=True)
    return {
        str(key): [str(value) for value in group["site_id"].tolist()]
        for key, group in frame.groupby("kinase", sort=False)
    }


def _top_n_overlap(expected: list[str], actual: list[str], n: int) -> float:
    expected_top = expected[:n]
    actual_top = actual[:n]
    if not expected_top:
        return 0.0
    return len(set(expected_top) & set(actual_top)) / float(len(expected_top))


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
    )
    result = dataset.process_core()

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
    )
    core = dataset.process_core()

    pred_mat = _read_indexed_table("predMat.csv")
    analyzer = KinaseActivityAnalyzer(pred_mat=pred_mat)
    result = analyzer.analyze(core.site_matrix.matrix)

    actual_weighted = result.weighted_activity.sort_index().sort_index(axis=1)
    expected_weighted = (
        _read_indexed_table("kinase_activity_matrix.csv")
        .sort_index()
        .sort_index(axis=1)
    )
    _assert_frame_close(actual_weighted, expected_weighted)

    actual_ksea = result.ksea_scores.sort_index().sort_index(axis=1)
    expected_ksea = (
        _read_indexed_table("ksea_scores.csv").sort_index().sort_index(axis=1)
    )
    _assert_frame_close(actual_ksea, expected_ksea)

    actual_ksea_counts = result.ksea_counts.sort_index()
    expected_ksea_counts_frame = _read_indexed_table("ksea_counts.csv")
    expected_ksea_counts = (
        expected_ksea_counts_frame.iloc[:, 0].sort_index()
        if not expected_ksea_counts_frame.empty
        else pd.Series(dtype=int)
    )
    expected_ksea_counts.name = actual_ksea_counts.name
    _assert_series_close(actual_ksea_counts, expected_ksea_counts)

    actual_target_counts = result.target_counts.sort_index()
    expected_target_counts = (
        _read_indexed_table("kinase_target_counts.csv").iloc[:, 0].sort_index()
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

    analyzer = KinaseActivityAnalyzer(pred_mat=pred_mat)
    result = analyzer.analyze(
        phospho_matrix, threshold=0.6, min_substrates=3, top_n_substrates=20
    )

    actual_weighted = result.weighted_activity.sort_index().sort_index(axis=1)
    expected_weighted = (
        _read_indexed_table("kinase_activity_matrix.csv", fixture_dir=R_FIXTURES_L6)
        .sort_index()
        .sort_index(axis=1)
    )
    _assert_frame_close(actual_weighted, expected_weighted)

    actual_ksea = result.ksea_scores.sort_index().sort_index(axis=1)
    expected_ksea = (
        _read_indexed_table("ksea_scores.csv", fixture_dir=R_FIXTURES_L6)
        .sort_index()
        .sort_index(axis=1)
    )
    _assert_frame_close(actual_ksea, expected_ksea)

    actual_ksea_counts = result.ksea_counts.sort_index()
    expected_ksea_counts_frame = _read_indexed_table(
        "ksea_counts.csv", fixture_dir=R_FIXTURES_L6
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
        _read_indexed_table("kinase_target_counts.csv", fixture_dir=R_FIXTURES_L6)
        .iloc[:, 0]
        .sort_index()
    )
    expected_target_counts.name = actual_target_counts.name
    _assert_series_close(actual_target_counts, expected_target_counts)


@pytest.mark.parity
def test_l6_native_profile_construction_matches_r_reference() -> None:
    _require_fixture_files(
        L6_FIXTURE_FILES + L6_NATIVE_FIXTURE_FILES,
        fixture_dir=R_FIXTURES_L6,
    )

    phospho_matrix = (
        _read_indexed_table("l6_phospho_matrix.csv", fixture_dir=R_FIXTURES_L6)
        .sort_index()
        .sort_index(axis=1)
    )
    substrate_map = _read_grouped_mapping(
        "native_substrate_map.csv",
        fixture_dir=R_FIXTURES_L6,
    )

    result = build_kinase_substrate_profiles(
        substrate_map=substrate_map,
        phospho_matrix=phospho_matrix,
        min_substrates=1,
    )

    actual_profile = result.profile_matrix.sort_index().sort_index(axis=1)
    expected_profile = (
        _read_indexed_table("native_profile_matrix.csv", fixture_dir=R_FIXTURES_L6)
        .sort_index()
        .sort_index(axis=1)
    )
    _assert_frame_close(actual_profile, expected_profile)

    actual_quantified = {
        kinase: sites for kinase, sites in result.quantified_substrates.items()
    }
    assert actual_quantified == substrate_map


@pytest.mark.parity
def test_l6_native_profile_scores_match_r_reference() -> None:
    _require_fixture_files(
        L6_FIXTURE_FILES + L6_NATIVE_FIXTURE_FILES,
        fixture_dir=R_FIXTURES_L6,
    )

    phospho_matrix = (
        _read_indexed_table("l6_phospho_matrix.csv", fixture_dir=R_FIXTURES_L6)
        .sort_index()
        .sort_index(axis=1)
    )
    expected_profile = (
        _read_indexed_table("native_profile_matrix.csv", fixture_dir=R_FIXTURES_L6)
        .sort_index()
        .sort_index(axis=1)
    )

    scorer = KinaseScorer(expected_profile)
    actual_profile_scores = (
        scorer.score_phosphosite_profiles(phospho_matrix)
        .sort_index()
        .sort_index(axis=1)
    )
    expected_profile_scores = (
        _read_indexed_table("native_profile_scores.csv", fixture_dir=R_FIXTURES_L6)
        .sort_index()
        .sort_index(axis=1)
    )
    _assert_frame_close(actual_profile_scores, expected_profile_scores, atol=1e-7)


@pytest.mark.parity
def test_l6_native_combined_scores_match_r_reference() -> None:
    _require_fixture_files(
        L6_FIXTURE_FILES + L6_NATIVE_FIXTURE_FILES,
        fixture_dir=R_FIXTURES_L6,
    )

    phospho_matrix = (
        _read_indexed_table("l6_phospho_matrix.csv", fixture_dir=R_FIXTURES_L6)
        .sort_index()
        .sort_index(axis=1)
    )
    substrate_map = _read_grouped_mapping(
        "native_substrate_map.csv",
        fixture_dir=R_FIXTURES_L6,
    )
    profile_result = build_kinase_substrate_profiles(
        substrate_map=substrate_map,
        phospho_matrix=phospho_matrix,
        min_substrates=1,
    )
    profile_scores = KinaseScorer(
        profile_result.profile_matrix
    ).score_phosphosite_profiles(phospho_matrix)
    motif_scores = _read_indexed_table(
        "native_motif_scores.csv", fixture_dir=R_FIXTURES_L6
    )
    motif_sizes = _read_named_series(
        "native_motif_sizes.csv",
        index_col="kinase",
        value_col="motif_size",
        fixture_dir=R_FIXTURES_L6,
    )

    actual_combined, actual_weights = combine_profile_and_motif_scores(
        motif_scores=motif_scores,
        profile_scores=profile_scores,
        motif_sizes=motif_sizes,
        profile_sizes=profile_result.substrate_counts.astype(float),
    )

    expected_combined = (
        _read_indexed_table("native_combined_scores.csv", fixture_dir=R_FIXTURES_L6)
        .sort_index()
        .sort_index(axis=1)
    )
    expected_weights = (
        _read_table("native_combined_weights.csv", fixture_dir=R_FIXTURES_L6)
        .set_index("kinase")
        .sort_index()
    )

    _assert_frame_close(
        actual_combined.sort_index().sort_index(axis=1),
        expected_combined,
        atol=1e-7,
    )
    _assert_frame_close(actual_weights.sort_index(), expected_weights, atol=1e-7)


@pytest.mark.parity
def test_l6_native_candidate_substrates_match_r_reference() -> None:
    _require_fixture_files(
        L6_FIXTURE_FILES + L6_NATIVE_FIXTURE_FILES,
        fixture_dir=R_FIXTURES_L6,
    )

    combined_scores = _read_indexed_table(
        "native_combined_scores.csv",
        fixture_dir=R_FIXTURES_L6,
    )
    actual = build_candidate_substrate_list(
        combined_scores,
        top=30,
        score_threshold=0.6,
        inclusion=5,
    )
    expected = _read_grouped_mapping(
        "native_candidate_substrates.csv",
        fixture_dir=R_FIXTURES_L6,
    )
    assert actual == expected


@pytest.mark.parity
def test_l6_native_prediction_rankings_agree_with_r_reference() -> None:
    _require_fixture_files(
        L6_FIXTURE_FILES + L6_NATIVE_FIXTURE_FILES,
        fixture_dir=R_FIXTURES_L6,
    )

    combined_scores = _read_indexed_table(
        "native_combined_scores.csv",
        fixture_dir=R_FIXTURES_L6,
    )
    expected_pred = _read_indexed_table("predMat.csv", fixture_dir=R_FIXTURES_L6)
    expected_top30 = _read_ranked_sites_by_kinase(
        "native_prediction_top30.csv",
        fixture_dir=R_FIXTURES_L6,
    )

    result = KinasePredictor().predict(
        combined_scores=combined_scores,
        ensemble_size=10,
        top=30,
        score_threshold=0.6,
        inclusion=5,
        n_iterations=5,
        random_state=1,
    )
    actual_pred = result.pred_matrix

    expected_kinases = sorted(expected_top30)
    assert expected_kinases
    assert set(expected_kinases) <= set(actual_pred.columns)

    common_kinases = [
        kinase for kinase in expected_kinases if kinase in expected_pred.columns
    ]
    assert len(common_kinases) == len(expected_kinases)

    top10_overlaps: list[float] = []
    top20_overlaps: list[float] = []
    top30_overlaps: list[float] = []
    rank_correlations: list[float] = []

    for kinase in common_kinases:
        actual_ranked_sites = (
            actual_pred.loc[:, kinase].sort_values(ascending=False).index.tolist()
        )
        expected_ranked_sites = expected_top30[kinase]

        top10_overlaps.append(
            _top_n_overlap(expected_ranked_sites, actual_ranked_sites, 10)
        )
        top20_overlaps.append(
            _top_n_overlap(expected_ranked_sites, actual_ranked_sites, 20)
        )
        top30_overlaps.append(
            _top_n_overlap(expected_ranked_sites, actual_ranked_sites, 30)
        )

        expected_ranks = expected_pred.loc[:, kinase].rank(
            ascending=False, method="average"
        )
        actual_ranks = actual_pred.loc[:, kinase].rank(
            ascending=False, method="average"
        )
        rank_correlations.append(expected_ranks.corr(actual_ranks, method="spearman"))

    assert pd.Series(rank_correlations).mean() >= 0.96
    assert pd.Series(top20_overlaps).mean() >= 0.85
    assert pd.Series(top30_overlaps).mean() >= 0.88
    assert sum(overlap >= 0.7 for overlap in top10_overlaps) >= 20
