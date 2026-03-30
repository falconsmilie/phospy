from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.prediction import (
    KinasePredictor,
    PredictionSamplingTrace,
    build_candidate_substrate_list,
    prediction_debug_trace_tables,
)
from phospy.scoring import combine_profile_and_motif_scores

ROOT = Path(__file__).resolve().parents[1]
FIXTURES_DIR = ROOT / "tests" / "fixtures"


@dataclass(frozen=True, slots=True)
class ReferenceSeamCase:
    name: str
    fixture_dir: Path
    profile_scores_file: str
    motif_scores_file: str
    motif_sizes_file: str
    profile_sizes_file: str | None
    combined_scores_file: str
    combined_weights_file: str
    candidate_substrates_file: str
    candidate_top: int
    candidate_score_threshold: float
    candidate_inclusion: int
    trace_dir: Path
    trace_kinases: tuple[str, ...]
    trace_top: int
    trace_score_threshold: float
    trace_inclusion: int
    trace_ensemble_size: int
    trace_n_iterations: int
    trace_random_state: int
    trace_debug_top_n: int


L6_CASE = ReferenceSeamCase(
    name="r_reference_l6",
    fixture_dir=FIXTURES_DIR / "r_reference_l6",
    profile_scores_file="native_profile_scores.csv",
    motif_scores_file="native_motif_scores.csv",
    motif_sizes_file="native_motif_sizes.csv",
    profile_sizes_file=None,
    combined_scores_file="native_combined_scores.csv",
    combined_weights_file="native_combined_weights.csv",
    candidate_substrates_file="native_candidate_substrates.csv",
    candidate_top=30,
    candidate_score_threshold=0.6,
    candidate_inclusion=5,
    trace_dir=FIXTURES_DIR / "r_reference_l6" / "prediction_trace",
    trace_kinases=("PRKAA1", "MAPK1", "MAPK9", "IRAK1", "TBK1", "LCK"),
    trace_top=30,
    trace_score_threshold=0.6,
    trace_inclusion=5,
    trace_ensemble_size=10,
    trace_n_iterations=5,
    trace_random_state=1,
    trace_debug_top_n=10,
)

SEAM_STRESS_CASE = ReferenceSeamCase(
    name="r_reference_l6_seam_stress",
    fixture_dir=FIXTURES_DIR / "r_reference_l6_seam_stress",
    profile_scores_file="profile_scores.csv",
    motif_scores_file="motif_scores.csv",
    motif_sizes_file="motif_sizes.csv",
    profile_sizes_file="profile_sizes.csv",
    combined_scores_file="combined_scores.csv",
    combined_weights_file="combined_weights.csv",
    candidate_substrates_file="candidate_substrates.csv",
    candidate_top=50,
    candidate_score_threshold=0.8,
    candidate_inclusion=20,
    trace_dir=FIXTURES_DIR / "r_reference_l6_seam_stress" / "prediction_trace",
    trace_kinases=("MAPK1", "PRKAA1"),
    trace_top=30,
    trace_score_threshold=0.6,
    trace_inclusion=5,
    trace_ensemble_size=10,
    trace_n_iterations=5,
    trace_random_state=1,
    trace_debug_top_n=10,
)

REFERENCE_SEAM_CASES = [L6_CASE, SEAM_STRESS_CASE]


def _read_indexed_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0)


def _read_named_series(path: Path, *, index_col: str, value_col: str) -> pd.Series:
    return pd.read_csv(path).set_index(index_col).loc[:, value_col]


def _read_grouped_mapping(path: Path) -> dict[str, list[str]]:
    frame = pd.read_csv(path)
    return {
        str(kinase): group.loc[:, "site_id"].astype(str).tolist()
        for kinase, group in frame.groupby("kinase", sort=False)
    }


def _read_profile_sizes(case: ReferenceSeamCase) -> pd.Series:
    if case.profile_sizes_file is not None:
        return _read_named_series(
            case.fixture_dir / case.profile_sizes_file,
            index_col="kinase",
            value_col="substrate_count",
        )

    substrate_map = pd.read_csv(case.fixture_dir / "native_substrate_map.csv")
    return substrate_map.groupby("kinase").size().rename("substrate_count")


def _assert_frame_close(
    actual: pd.DataFrame,
    expected: pd.DataFrame,
    *,
    atol: float = 1e-8,
    rtol: float = 1e-6,
) -> None:
    pdt.assert_frame_equal(
        actual.sort_index().sort_index(axis=1),
        expected.sort_index().sort_index(axis=1),
        check_dtype=False,
        check_index_type=False,
        check_column_type=False,
        check_names=False,
        atol=atol,
        rtol=rtol,
    )


@pytest.mark.parity
@pytest.mark.parametrize(
    "case",
    REFERENCE_SEAM_CASES,
    ids=[case.name for case in REFERENCE_SEAM_CASES],
)
def test_reference_combined_score_seam_matches_expected(
    case: ReferenceSeamCase,
) -> None:
    profile_scores = _read_indexed_table(case.fixture_dir / case.profile_scores_file)
    motif_scores = _read_indexed_table(case.fixture_dir / case.motif_scores_file)
    motif_sizes = _read_named_series(
        case.fixture_dir / case.motif_sizes_file,
        index_col="kinase",
        value_col="motif_size",
    )
    profile_sizes = _read_profile_sizes(case).astype(float)

    actual_combined, actual_weights = combine_profile_and_motif_scores(
        motif_scores=motif_scores,
        profile_scores=profile_scores,
        motif_sizes=motif_sizes,
        profile_sizes=profile_sizes,
    )

    expected_combined = _read_indexed_table(
        case.fixture_dir / case.combined_scores_file
    )
    expected_weights = (
        pd.read_csv(case.fixture_dir / case.combined_weights_file)
        .set_index("kinase")
        .sort_index()
    )

    _assert_frame_close(actual_combined, expected_combined, atol=1e-8, rtol=1e-6)
    _assert_frame_close(actual_weights, expected_weights, atol=1e-8, rtol=1e-6)


@pytest.mark.parity
@pytest.mark.parametrize(
    "case",
    REFERENCE_SEAM_CASES,
    ids=[case.name for case in REFERENCE_SEAM_CASES],
)
def test_reference_candidate_selection_seam_matches_expected(
    case: ReferenceSeamCase,
) -> None:
    combined_scores = _read_indexed_table(case.fixture_dir / case.combined_scores_file)

    actual = build_candidate_substrate_list(
        combined_scores=combined_scores,
        top=case.candidate_top,
        score_threshold=case.candidate_score_threshold,
        inclusion=case.candidate_inclusion,
    )
    expected = _read_grouped_mapping(case.fixture_dir / case.candidate_substrates_file)

    assert actual == expected


@pytest.mark.parity
@pytest.mark.parametrize(
    "case",
    REFERENCE_SEAM_CASES,
    ids=[case.name for case in REFERENCE_SEAM_CASES],
)
def test_reference_replayed_prediction_trace_matches_expected_sampling_path(
    case: ReferenceSeamCase,
    tmp_path: Path,
) -> None:
    combined_scores = _read_indexed_table(case.fixture_dir / case.combined_scores_file)
    expected_candidates = pd.read_csv(case.trace_dir / "trace_candidates.csv")
    expected_initial = pd.read_csv(case.trace_dir / "trace_initial_negatives.csv")
    expected_samples = pd.read_csv(case.trace_dir / "trace_iteration_samples.csv")

    trace_selected_counts = (
        expected_candidates.loc[expected_candidates.loc[:, "selected_candidate"]]
        .groupby("kinase")
        .size()
        .to_dict()
    )

    current_candidates = build_candidate_substrate_list(
        combined_scores,
        top=case.trace_top,
        score_threshold=case.trace_score_threshold,
        inclusion=case.trace_inclusion,
    )

    eligible_kinases: list[str] = []
    for kinase in case.trace_kinases:
        trace_count = int(trace_selected_counts.get(kinase, 0))
        current_count = len(current_candidates.get(kinase, []))
        initial_counts = (
            expected_initial.loc[
                expected_initial.loc[:, "kinase"].astype(str) == kinase
            ]
            .groupby("ensemble")
            .size()
        )
        sample_counts = (
            expected_samples.loc[
                expected_samples.loc[:, "kinase"].astype(str) == kinase
            ]
            .groupby(["ensemble", "iteration", "class_label"])
            .size()
        )
        if (
            trace_count > 0
            and current_count == trace_count
            and not initial_counts.empty
            and bool((initial_counts == trace_count).all())
            and not sample_counts.empty
            and bool((sample_counts == trace_count).all())
        ):
            eligible_kinases.append(kinase)

    assert eligible_kinases, f"No replay-aligned kinases found for fixture {case.name}"

    sampling_trace = PredictionSamplingTrace.from_trace_directory(
        case.trace_dir
    ).subset_kinases(eligible_kinases)

    result = KinasePredictor(svm_mode="r_parity").predict(
        combined_scores=combined_scores,
        ensemble_size=case.trace_ensemble_size,
        top=case.trace_top,
        score_threshold=case.trace_score_threshold,
        inclusion=case.trace_inclusion,
        n_iterations=case.trace_n_iterations,
        random_state=case.trace_random_state,
        capture_debug_trace=True,
        debug_kinases=eligible_kinases,
        debug_top_n=case.trace_debug_top_n,
        sampling_trace=sampling_trace,
        trace_level="full",
        trace_sink=tmp_path / "python_trace_output",
    )
    actual_tables = prediction_debug_trace_tables(result)

    actual_initial = actual_tables["trace_initial_negatives"].copy()
    actual_samples = actual_tables["trace_iteration_samples"].copy()
    expected_initial = expected_initial.loc[
        expected_initial.loc[:, "kinase"].astype(str).isin(eligible_kinases)
    ].copy()
    expected_samples = expected_samples.loc[
        expected_samples.loc[:, "kinase"].astype(str).isin(eligible_kinases)
    ].copy()

    actual_samples = actual_samples.astype({"class_label": str})
    expected_samples = expected_samples.astype({"class_label": str})
    for table in (actual_initial, expected_initial):
        table.loc[:, "kinase"] = table.loc[:, "kinase"].astype(str)
        table.loc[:, "site"] = table.loc[:, "site"].astype(str)
    for table in (actual_samples, expected_samples):
        table.loc[:, "kinase"] = table.loc[:, "kinase"].astype(str)
        table.loc[:, "site"] = table.loc[:, "site"].astype(str)

    actual_initial = actual_initial.sort_values(
        ["kinase", "ensemble", "draw", "site"], kind="mergesort"
    ).reset_index(drop=True)
    expected_initial = expected_initial.sort_values(
        ["kinase", "ensemble", "draw", "site"], kind="mergesort"
    ).reset_index(drop=True)
    actual_samples = actual_samples.sort_values(
        ["kinase", "ensemble", "iteration", "class_label", "draw", "site"],
        kind="mergesort",
    ).reset_index(drop=True)
    expected_samples = expected_samples.sort_values(
        ["kinase", "ensemble", "iteration", "class_label", "draw", "site"],
        kind="mergesort",
    ).reset_index(drop=True)

    pdt.assert_frame_equal(actual_initial, expected_initial, check_dtype=False)
    pdt.assert_frame_equal(actual_samples, expected_samples, check_dtype=False)

    actual_probabilities = actual_tables["trace_iteration_probabilities"].copy()
    actual_decisions = actual_tables["trace_iteration_decision_values"].copy()
    actual_final_top = actual_tables["trace_final_ensemble_top"].copy()

    expected_probabilities = pd.read_csv(
        case.trace_dir / "trace_iteration_probabilities.csv"
    )
    expected_probabilities = expected_probabilities.loc[
        expected_probabilities.loc[:, "kinase"].astype(str).isin(eligible_kinases)
    ].copy()
    expected_decisions = pd.read_csv(
        case.trace_dir / "trace_iteration_decision_values.csv"
    )
    expected_decisions = expected_decisions.loc[
        expected_decisions.loc[:, "kinase"].astype(str).isin(eligible_kinases)
    ].copy()
    expected_final_top = pd.read_csv(case.trace_dir / "trace_final_ensemble_top.csv")
    expected_final_top = expected_final_top.loc[
        expected_final_top.loc[:, "kinase"].astype(str).isin(eligible_kinases)
    ].copy()

    for table in (
        actual_probabilities,
        expected_probabilities,
        actual_decisions,
        expected_decisions,
        actual_final_top,
        expected_final_top,
    ):
        table.loc[:, "kinase"] = table.loc[:, "kinase"].astype(str)
        table.loc[:, "site"] = table.loc[:, "site"].astype(str)
    actual_probabilities = actual_probabilities.astype({"label": str})
    expected_probabilities = expected_probabilities.astype({"label": str})
    actual_decisions = actual_decisions.astype({"label": str})
    expected_decisions = expected_decisions.astype({"label": str})

    merged_probabilities = actual_probabilities.merge(
        expected_probabilities,
        on=["kinase", "ensemble", "iteration", "site", "label"],
        suffixes=("_py", "_r"),
    )
    merged_decisions = actual_decisions.merge(
        expected_decisions,
        on=["kinase", "ensemble", "iteration", "site", "label"],
        suffixes=("_py", "_r"),
    )
    merged_final_top = actual_final_top.merge(
        expected_final_top,
        on=["kinase", "ensemble", "rank"],
        suffixes=("_py", "_r"),
    )

    assert not merged_probabilities.empty
    assert not merged_decisions.empty
    assert not merged_final_top.empty

    prob_class1_corr = float(
        merged_probabilities.loc[:, "prob_class_1_py"].corr(
            merged_probabilities.loc[:, "prob_class_1_r"]
        )
    )
    prob_mae = float(
        (
            merged_probabilities.loc[:, "prob_class_1_py"]
            - merged_probabilities.loc[:, "prob_class_1_r"]
        )
        .abs()
        .mean()
    )
    decision_corr = float(
        merged_decisions.loc[:, "decision_value_class_1_py"].corr(
            merged_decisions.loc[:, "decision_value_class_1_r"]
        )
    )
    decision_mae = float(
        (
            merged_decisions.loc[:, "decision_value_class_1_py"]
            - merged_decisions.loc[:, "decision_value_class_1_r"]
        )
        .abs()
        .mean()
    )
    final_top_site_matches = int(
        (merged_final_top.loc[:, "site_py"] == merged_final_top.loc[:, "site_r"]).sum()
    )
    final_top_total = int(len(merged_final_top))
    final_top_prob_mae = float(
        (
            merged_final_top.loc[:, "prob_class_1_py"]
            - merged_final_top.loc[:, "prob_class_1_r"]
        )
        .abs()
        .mean()
    )

    assert prob_class1_corr >= 0.998
    assert prob_mae <= 0.015
    assert decision_corr >= 0.999999
    assert decision_mae <= 1e-12
    assert final_top_site_matches >= int(final_top_total * 0.95)
    assert final_top_prob_mae <= 0.02
