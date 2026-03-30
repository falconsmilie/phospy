from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import KinaseActivityAnalyzer, PhosphoDataset
from phospy.prediction import (
    KinasePredictor,
    PredictionSamplingTrace,
    build_candidate_substrate_list,
    prediction_debug_trace_tables,
)
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


def _env_flag(name: str) -> bool:
    value = os.getenv(name, "")
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _show_parity_metrics() -> bool:
    return _env_flag("PHOSPY_SHOW_PARITY")


def _show_profile_construction_metrics() -> bool:
    return _show_parity_metrics() and _env_flag("PHOSPY_SHOW_PROFILE_CONSTRUCTION")


def _show_prediction_mode_comparison_metrics() -> bool:
    return _show_parity_metrics() and _env_flag(
        "PHOSPY_SHOW_PREDICTION_MODE_COMPARISON"
    )


def _show_replayed_prediction_mode_comparison_metrics() -> bool:
    return _show_parity_metrics() and _env_flag(
        "PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON"
    )


def _mean_abs_diff(actual: pd.DataFrame, expected: pd.DataFrame) -> float:
    aligned_actual = actual.sort_index().sort_index(axis=1)
    aligned_expected = expected.sort_index().sort_index(axis=1)
    diff = (aligned_actual - aligned_expected).abs().to_numpy().ravel()
    if diff.size == 0:
        return 0.0
    return float(pd.Series(diff).mean())


def _max_abs_diff(actual: pd.DataFrame, expected: pd.DataFrame) -> float:
    aligned_actual = actual.sort_index().sort_index(axis=1)
    aligned_expected = expected.sort_index().sort_index(axis=1)
    diff = (aligned_actual - aligned_expected).abs().to_numpy().ravel()
    if diff.size == 0:
        return 0.0
    return float(pd.Series(diff).max())


def _mean_column_correlation(
    actual: pd.DataFrame, expected: pd.DataFrame, *, method: str = "pearson"
) -> float:
    common_index = actual.index.intersection(expected.index)
    common_columns = actual.columns.intersection(expected.columns)
    corrs: list[float] = []
    for column in common_columns:
        actual_col = actual.loc[common_index, column]
        expected_col = expected.loc[common_index, column]
        corr = actual_col.corr(expected_col, method=method)
        if pd.notna(corr):
            corrs.append(float(corr))
    if not corrs:
        return float("nan")
    return sum(corrs) / len(corrs)


def _maybe_print_metrics(
    title: str, lines: list[str], *, optional: bool = False
) -> None:
    if optional:
        if not _show_profile_construction_metrics():
            return
    elif not _show_parity_metrics():
        return

    print(f"\n{title}:")
    for line in lines:
        print(f"  {line}")


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

    _maybe_print_metrics(
        "Optional profile-construction parity metrics",
        [
            f"kinases compared: {len(actual_profile.index.intersection(expected_profile.index))}",
            f"profile matrix shape: {actual_profile.shape} vs {expected_profile.shape}",
            f"mean per-kinase Pearson correlation: {_mean_column_correlation(actual_profile.T, expected_profile.T, method='pearson') * 100:.2f}%",
            f"mean absolute difference: {_mean_abs_diff(actual_profile, expected_profile):.6g}",
            f"max absolute difference: {_max_abs_diff(actual_profile, expected_profile):.6g}",
        ],
        optional=True,
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

    _maybe_print_metrics(
        "Profile-scoring parity metrics",
        [
            f"sites compared: {len(actual_profile_scores.index)}",
            f"kinases compared: {len(actual_profile_scores.columns)}",
            f"mean per-kinase Pearson correlation: {_mean_column_correlation(actual_profile_scores, expected_profile_scores, method='pearson') * 100:.2f}%",
            f"mean per-kinase Spearman correlation: {_mean_column_correlation(actual_profile_scores, expected_profile_scores, method='spearman') * 100:.2f}%",
            f"mean absolute difference: {_mean_abs_diff(actual_profile_scores, expected_profile_scores):.6g}",
            f"max absolute difference: {_max_abs_diff(actual_profile_scores, expected_profile_scores):.6g}",
        ],
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

    actual_combined = actual_combined.sort_index().sort_index(axis=1)

    _maybe_print_metrics(
        "Combined-score parity metrics",
        [
            f"sites compared: {len(actual_combined.index)}",
            f"kinases compared: {len(actual_combined.columns)}",
            f"mean per-kinase Pearson correlation: {_mean_column_correlation(actual_combined, expected_combined, method='pearson') * 100:.2f}%",
            f"mean per-kinase Spearman correlation: {_mean_column_correlation(actual_combined, expected_combined, method='spearman') * 100:.2f}%",
            f"mean absolute difference: {_mean_abs_diff(actual_combined, expected_combined):.6g}",
            f"max absolute difference: {_max_abs_diff(actual_combined, expected_combined):.6g}",
            f"mean weight absolute difference: {_mean_abs_diff(actual_weights.sort_index(), expected_weights):.6g}",
        ],
    )

    _assert_frame_close(
        actual_combined,
        expected_combined,
        atol=1e-7,
    )
    _assert_frame_close(actual_weights.sort_index(), expected_weights, atol=1e-7)


def _debug_trace_to_tables(result) -> dict[str, pd.DataFrame]:
    return prediction_debug_trace_tables(result)


def _read_prediction_trace_table(name: str) -> pd.DataFrame:
    return pd.read_csv(R_FIXTURES_L6 / "prediction_trace" / name)


def _replayed_prediction_trace_metrics(
    *, svm_mode: str
) -> dict[str, float | int | str]:
    combined_scores = _read_indexed_table(
        "native_combined_scores.csv",
        fixture_dir=R_FIXTURES_L6,
    )
    current_candidates = build_candidate_substrate_list(
        combined_scores,
        top=30,
        score_threshold=0.6,
        inclusion=5,
    )
    expected_candidates = _read_prediction_trace_table("trace_candidates.csv")
    initial_trace = _read_prediction_trace_table("trace_initial_negatives.csv")
    sample_trace = _read_prediction_trace_table("trace_iteration_samples.csv")

    trace_selected_counts = (
        expected_candidates.loc[expected_candidates["selected_candidate"]]
        .groupby("kinase")
        .size()
        .to_dict()
    )

    eligible_kinases: list[str] = []
    skipped_kinases: list[str] = []
    for kinase in sorted(trace_selected_counts):
        current_count = len(current_candidates.get(kinase, []))
        trace_count = int(trace_selected_counts[kinase])
        initial_counts = (
            initial_trace.loc[initial_trace["kinase"] == kinase]
            .groupby("ensemble")
            .size()
        )
        sample_counts = (
            sample_trace.loc[sample_trace["kinase"] == kinase]
            .groupby(["ensemble", "iteration", "class_label"])
            .size()
        )
        initial_ok = not initial_counts.empty and bool(
            (initial_counts == trace_count).all()
        )
        sample_ok = not sample_counts.empty and bool(
            (sample_counts == trace_count).all()
        )
        if current_count == trace_count and initial_ok and sample_ok:
            eligible_kinases.append(str(kinase))
        else:
            skipped_kinases.append(str(kinase))

    if not eligible_kinases:
        pytest.skip(
            "No replay-aligned trace kinases are available. Regenerate the R "
            "prediction trace fixtures so candidate counts and sampling rows "
            "match the current native candidate set."
        )

    sampling_trace = PredictionSamplingTrace.from_trace_directory(
        R_FIXTURES_L6 / "prediction_trace"
    ).subset_kinases(eligible_kinases)

    with tempfile.TemporaryDirectory(
        prefix="phospy_trace_metrics_"
    ) as trace_output_dir:
        result = KinasePredictor(svm_mode=svm_mode).predict(
            combined_scores=combined_scores,
            ensemble_size=10,
            top=30,
            score_threshold=0.6,
            inclusion=5,
            n_iterations=5,
            random_state=1,
            capture_debug_trace=True,
            debug_kinases=eligible_kinases,
            debug_top_n=10,
            sampling_trace=sampling_trace,
            trace_level="full",
            trace_sink=trace_output_dir,
        )
        actual_tables = _debug_trace_to_tables(result)
    expected_initial = _normalize_numeric_frame(
        _sort_table(
            _read_prediction_trace_table("trace_initial_negatives.csv").loc[
                lambda df: df["kinase"].astype(str).isin(eligible_kinases)
            ],
            ["kinase", "ensemble", "draw", "site"],
        )
    )
    actual_initial = _normalize_numeric_frame(
        _sort_table(
            actual_tables["trace_initial_negatives"],
            ["kinase", "ensemble", "draw", "site"],
        )
    )
    expected_samples = _normalize_numeric_frame(
        _sort_table(
            _read_prediction_trace_table("trace_iteration_samples.csv").loc[
                lambda df: df["kinase"].astype(str).isin(eligible_kinases)
            ],
            ["kinase", "ensemble", "iteration", "class_label", "draw", "site"],
        )
    )
    actual_samples = _normalize_numeric_frame(
        _sort_table(
            actual_tables["trace_iteration_samples"],
            ["kinase", "ensemble", "iteration", "class_label", "draw", "site"],
        )
    )

    initial_exact_matches = int(actual_initial.eq(expected_initial).all(axis=1).sum())
    initial_total_rows = int(len(expected_initial))
    sample_exact_matches = int(actual_samples.eq(expected_samples).all(axis=1).sum())
    sample_total_rows = int(len(expected_samples))

    expected_prob = _normalize_numeric_frame(
        _sort_table(
            _read_prediction_trace_table("trace_iteration_probabilities.csv").loc[
                lambda df: df["kinase"].astype(str).isin(eligible_kinases)
            ],
            ["kinase", "ensemble", "iteration", "site"],
        )
    )
    actual_prob = _normalize_numeric_frame(
        _sort_table(
            actual_tables["trace_iteration_probabilities"],
            ["kinase", "ensemble", "iteration", "site"],
        )
    )
    merged_prob = actual_prob.merge(
        expected_prob,
        on=["kinase", "ensemble", "iteration", "site", "label"],
        suffixes=("_py", "_r"),
        validate="one_to_one",
    )

    expected_decision = _normalize_numeric_frame(
        _sort_table(
            _read_prediction_trace_table("trace_iteration_decision_values.csv").loc[
                lambda df: df["kinase"].astype(str).isin(eligible_kinases)
            ],
            ["kinase", "ensemble", "iteration", "site"],
        )
    )
    actual_decision = _normalize_numeric_frame(
        _sort_table(
            actual_tables["trace_iteration_decision_values"],
            ["kinase", "ensemble", "iteration", "site"],
        )
    )
    merged_decision = actual_decision.merge(
        expected_decision,
        on=["kinase", "ensemble", "iteration", "site", "label"],
        suffixes=("_py", "_r"),
        validate="one_to_one",
    )

    expected_top = _normalize_numeric_frame(
        _sort_table(
            _read_prediction_trace_table("trace_final_ensemble_top.csv").loc[
                lambda df: df["kinase"].astype(str).isin(eligible_kinases)
            ],
            ["kinase", "ensemble", "rank"],
        )
    )
    actual_top = _normalize_numeric_frame(
        _sort_table(
            actual_tables["trace_final_ensemble_top"],
            ["kinase", "ensemble", "rank"],
        )
    )
    merged_top = actual_top.merge(
        expected_top,
        on=["kinase", "ensemble", "rank"],
        suffixes=("_py", "_r"),
        validate="one_to_one",
    )

    prob_class1_corr = merged_prob["prob_class_1_py"].corr(
        merged_prob["prob_class_1_r"], method="pearson"
    )
    prob_class2_corr = merged_prob["prob_class_2_py"].corr(
        merged_prob["prob_class_2_r"], method="pearson"
    )
    final_top_site_matches = int((merged_top["site_py"] == merged_top["site_r"]).sum())
    final_top_total = int(len(merged_top))

    return {
        "n_kinases": len(eligible_kinases),
        "trace_kinases": ", ".join(eligible_kinases),
        "skipped_trace_kinases": ", ".join(skipped_kinases),
        "initial_exact_matches": initial_exact_matches,
        "initial_total_rows": initial_total_rows,
        "sample_exact_matches": sample_exact_matches,
        "sample_total_rows": sample_total_rows,
        "iteration_prob_class1_corr": float(prob_class1_corr),
        "iteration_prob_class2_corr": float(prob_class2_corr),
        "iteration_decision_class1_corr": float(
            merged_decision["decision_value_class_1_py"].corr(
                merged_decision["decision_value_class_1_r"], method="pearson"
            )
        ),
        "iteration_decision_mae": float(
            (
                merged_decision["decision_value_class_1_py"]
                - merged_decision["decision_value_class_1_r"]
            )
            .abs()
            .mean()
        ),
        "iteration_prob_mae": float(
            (merged_prob["prob_class_1_py"] - merged_prob["prob_class_1_r"])
            .abs()
            .mean()
        ),
        "iteration_prob_max_abs_diff": float(
            (merged_prob["prob_class_1_py"] - merged_prob["prob_class_1_r"]).abs().max()
        ),
        "final_top_site_matches": final_top_site_matches,
        "final_top_total": final_top_total,
        "final_top_prob_mae": float(
            (merged_top["prob_class_1_py"] - merged_top["prob_class_1_r"]).abs().mean()
        ),
    }


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


def _prediction_parity_metrics(*, svm_mode: str) -> dict[str, float | int]:
    combined_scores = _read_indexed_table(
        "native_combined_scores.csv",
        fixture_dir=R_FIXTURES_L6,
    )
    expected_pred = _read_indexed_table("predMat.csv", fixture_dir=R_FIXTURES_L6)
    expected_top30 = _read_ranked_sites_by_kinase(
        "native_prediction_top30.csv",
        fixture_dir=R_FIXTURES_L6,
    )

    result = KinasePredictor(svm_mode=svm_mode).predict(
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

    mean_spearman = float(pd.Series(rank_correlations).mean())
    mean_top10_overlap = float(pd.Series(top10_overlaps).mean())
    mean_top20_overlap = float(pd.Series(top20_overlaps).mean())
    mean_top30_overlap = float(pd.Series(top30_overlaps).mean())
    n_good_top10 = int(sum(overlap >= 0.7 for overlap in top10_overlaps))
    n_kinases = len(common_kinases)

    return {
        "mean_spearman": mean_spearman,
        "mean_top10_overlap": mean_top10_overlap,
        "mean_top20_overlap": mean_top20_overlap,
        "mean_top30_overlap": mean_top30_overlap,
        "n_good_top10": n_good_top10,
        "n_kinases": n_kinases,
    }


@pytest.mark.parity
def test_l6_native_prediction_rankings_agree_with_r_reference() -> None:
    _require_fixture_files(
        L6_FIXTURE_FILES + L6_NATIVE_FIXTURE_FILES,
        fixture_dir=R_FIXTURES_L6,
    )

    metrics = _prediction_parity_metrics(svm_mode="default")

    _maybe_print_metrics(
        "Prediction parity metrics",
        [
            "svm_mode: default",
            f"kinases compared: {metrics['n_kinases']}",
            f"mean Spearman rank agreement: {metrics['mean_spearman'] * 100:.2f}%",
            f"mean top-10 overlap: {metrics['mean_top10_overlap'] * 100:.2f}%",
            f"mean top-20 overlap: {metrics['mean_top20_overlap'] * 100:.2f}%",
            f"mean top-30 overlap: {metrics['mean_top30_overlap'] * 100:.2f}%",
            f"kinases with top-10 overlap >= 70%: {metrics['n_good_top10']}/{metrics['n_kinases']}",
        ],
    )

    assert metrics["mean_spearman"] >= 0.96
    assert metrics["mean_top20_overlap"] >= 0.85
    assert metrics["mean_top30_overlap"] >= 0.88
    assert metrics["n_good_top10"] >= 20


@pytest.mark.parity
def test_l6_native_prediction_mode_comparison_metrics() -> None:
    _require_fixture_files(
        L6_FIXTURE_FILES + L6_NATIVE_FIXTURE_FILES,
        fixture_dir=R_FIXTURES_L6,
    )

    default_metrics = _prediction_parity_metrics(svm_mode="default")
    r_parity_metrics = _prediction_parity_metrics(svm_mode="r_parity")

    if _show_prediction_mode_comparison_metrics():
        _maybe_print_metrics(
            "Prediction parity mode comparison",
            [
                f"default mean Spearman rank agreement: {default_metrics['mean_spearman'] * 100:.2f}%",
                f"r_parity mean Spearman rank agreement: {r_parity_metrics['mean_spearman'] * 100:.2f}%",
                f"default mean top-10 overlap: {default_metrics['mean_top10_overlap'] * 100:.2f}%",
                f"r_parity mean top-10 overlap: {r_parity_metrics['mean_top10_overlap'] * 100:.2f}%",
                f"default mean top-20 overlap: {default_metrics['mean_top20_overlap'] * 100:.2f}%",
                f"r_parity mean top-20 overlap: {r_parity_metrics['mean_top20_overlap'] * 100:.2f}%",
                f"default mean top-30 overlap: {default_metrics['mean_top30_overlap'] * 100:.2f}%",
                f"r_parity mean top-30 overlap: {r_parity_metrics['mean_top30_overlap'] * 100:.2f}%",
            ],
        )

    assert default_metrics["n_kinases"] == r_parity_metrics["n_kinases"]


@pytest.mark.parity
def test_l6_replayed_prediction_trace_matches_r_sampling_path() -> None:
    _require_fixture_files(
        L6_FIXTURE_FILES + L6_NATIVE_FIXTURE_FILES,
        fixture_dir=R_FIXTURES_L6,
    )

    metrics = _replayed_prediction_trace_metrics(svm_mode="r_parity")

    _maybe_print_metrics(
        "Replayed prediction trace parity metrics",
        [
            "svm_mode: r_parity",
            f"kinases compared: {metrics['n_kinases']}",
            f"trace kinases used: {metrics['trace_kinases']}",
            f"trace kinases skipped: {metrics['skipped_trace_kinases'] or 'none'}",
            f"initial negative exact matches: {metrics['initial_exact_matches']}/{metrics['initial_total_rows']}",
            f"iteration sample exact matches: {metrics['sample_exact_matches']}/{metrics['sample_total_rows']}",
            f"iteration prob class-1 Pearson correlation: {metrics['iteration_prob_class1_corr'] * 100:.3f}%",
            f"iteration prob class-2 Pearson correlation: {metrics['iteration_prob_class2_corr'] * 100:.3f}%",
            f"iteration decision class-1 Pearson correlation: {metrics['iteration_decision_class1_corr'] * 100:.3f}%",
            f"iteration decision class-1 mean absolute difference: {metrics['iteration_decision_mae']:.6g}",
            f"iteration prob class-1 mean absolute difference: {metrics['iteration_prob_mae']:.6g}",
            f"final top-site matches: {metrics['final_top_site_matches']}/{metrics['final_top_total']}",
            f"final top class-1 mean absolute difference: {metrics['final_top_prob_mae']:.6g}",
        ],
    )

    assert metrics["initial_exact_matches"] == metrics["initial_total_rows"]
    assert metrics["sample_exact_matches"] == metrics["sample_total_rows"]
    assert metrics["iteration_decision_class1_corr"] >= 0.999999
    assert metrics["iteration_decision_mae"] <= 1e-12
    assert metrics["iteration_prob_class1_corr"] >= 0.998
    assert metrics["iteration_prob_mae"] <= 0.015
    assert metrics["final_top_site_matches"] >= 95


@pytest.mark.parity
def test_l6_replayed_prediction_mode_comparison_metrics() -> None:
    _require_fixture_files(
        L6_FIXTURE_FILES + L6_NATIVE_FIXTURE_FILES,
        fixture_dir=R_FIXTURES_L6,
    )

    default_metrics = _replayed_prediction_trace_metrics(svm_mode="default")
    r_parity_metrics = _replayed_prediction_trace_metrics(svm_mode="r_parity")

    if _show_replayed_prediction_mode_comparison_metrics():
        _maybe_print_metrics(
            "Replayed prediction parity mode comparison",
            [
                f"default trace kinases used: {default_metrics['trace_kinases']}",
                f"r_parity trace kinases used: {r_parity_metrics['trace_kinases']}",
                f"trace kinases skipped: {r_parity_metrics['skipped_trace_kinases'] or 'none'}",
                f"default iteration prob class-1 Pearson correlation: {default_metrics['iteration_prob_class1_corr'] * 100:.3f}%",
                f"r_parity iteration prob class-1 Pearson correlation: {r_parity_metrics['iteration_prob_class1_corr'] * 100:.3f}%",
                f"default iteration prob mean absolute difference: {default_metrics['iteration_prob_mae']:.6g}",
                f"r_parity iteration prob mean absolute difference: {r_parity_metrics['iteration_prob_mae']:.6g}",
                f"default final top-site matches: {default_metrics['final_top_site_matches']}/{default_metrics['final_top_total']}",
                f"r_parity final top-site matches: {r_parity_metrics['final_top_site_matches']}/{r_parity_metrics['final_top_total']}",
                f"default final top class-1 mean absolute difference: {default_metrics['final_top_prob_mae']:.6g}",
                f"r_parity final top class-1 mean absolute difference: {r_parity_metrics['final_top_prob_mae']:.6g}",
            ],
        )

    assert (
        default_metrics["initial_total_rows"] == r_parity_metrics["initial_total_rows"]
    )
    assert default_metrics["sample_total_rows"] == r_parity_metrics["sample_total_rows"]
