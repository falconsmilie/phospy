#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
import tempfile
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import pandas as pd
import pandas.testing as pdt

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from phospy.prediction import (  # noqa: E402
    KinasePredictor,
    PredictionSamplingTrace,
    build_candidate_substrate_list,
    prediction_debug_trace_tables,
)
from phospy.workflow import PredMatWorkflow, SignalomeWorkflow  # noqa: E402

ModeName = Literal["default", "r_parity"]
MODES: tuple[ModeName, ModeName] = ("default", "r_parity")

EXAMPLE_DATA = REPO_ROOT / "examples" / "data"
FIXTURES_ROOT = REPO_ROOT / "tests" / "fixtures"
R_FIXTURES_L6 = FIXTURES_ROOT / "r_reference_l6"
WORKFLOW_FIXTURES = FIXTURES_ROOT / "public_workflow_reference"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "benchmarks" / "reports" / "latest"

PREDMAT_BENCHMARKS = {
    "default": "predmat_default.csv",
    "r_parity": "predmat_r_parity.csv",
}
SIGNALOME_BENCHMARKS = {
    "modules": "signalome_modules.csv",
    "map_modules": "signalome_map_modules.csv",
    "network_nodes": "signalome_network_nodes.csv",
    "network_edges": "signalome_network_edges.csv",
}


@dataclass(frozen=True)
class RuntimeStats:
    runs: list[float]
    mean_seconds: float
    median_seconds: float
    min_seconds: float
    max_seconds: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "runs_seconds": self.runs,
            "mean_seconds": self.mean_seconds,
            "median_seconds": self.median_seconds,
            "min_seconds": self.min_seconds,
            "max_seconds": self.max_seconds,
        }


@dataclass(frozen=True)
class ModeBenchmarkResult:
    mode: ModeName
    metrics: dict[str, Any]
    thresholds: dict[str, bool]
    runtime: RuntimeStats

    def to_dict(self) -> dict[str, Any]:
        return {
            "mode": self.mode,
            "metrics": self.metrics,
            "thresholds": self.thresholds,
            "runtime": self.runtime.to_dict(),
        }


@dataclass(frozen=True)
class BenchmarkCaseResult:
    benchmark: str
    fixture_family: str
    protected_seam: str
    metric_class: str
    modes: dict[ModeName, ModeBenchmarkResult]
    comparisons: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return {
            "benchmark": self.benchmark,
            "fixture_family": self.fixture_family,
            "protected_seam": self.protected_seam,
            "metric_class": self.metric_class,
            "modes": {mode: result.to_dict() for mode, result in self.modes.items()},
            "comparisons": self.comparisons,
        }


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    fixture_family: str
    protected_seam: str
    metric_class: str
    runner: Callable[[ModeName], dict[str, Any]]
    threshold_evaluator: Callable[[dict[str, Any]], dict[str, bool]]
    comparison_builder: Callable[[dict[ModeName, ModeBenchmarkResult]], dict[str, Any]]


def _require_paths(paths: list[Path]) -> None:
    missing = [str(path.relative_to(REPO_ROOT)) for path in paths if not path.exists()]
    if missing:
        joined = ", ".join(missing)
        raise FileNotFoundError(
            "Required benchmark fixtures are missing. Generate the parity fixtures "
            f"before running the benchmark. Missing: {joined}"
        )


def _read_indexed_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0)


def _read_unindexed_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def _read_ranked_sites_by_kinase(path: Path) -> dict[str, list[str]]:
    frame = pd.read_csv(path)
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
    present = [column for column in columns if column in df.columns]
    if not present:
        return df.reset_index(drop=True)
    return df.sort_values(present).reset_index(drop=True)


def _normalize_numeric_frame(df: pd.DataFrame) -> pd.DataFrame:
    normalized = df.copy()
    for column in normalized.columns:
        try:
            normalized[column] = pd.to_numeric(normalized[column])
        except (TypeError, ValueError):
            continue
    return normalized


def _max_abs_diff(actual: pd.DataFrame, expected: pd.DataFrame) -> float:
    aligned_actual = actual.sort_index().sort_index(axis=1)
    aligned_expected = expected.sort_index().sort_index(axis=1)
    diff = (aligned_actual - aligned_expected).abs().to_numpy().ravel()
    if diff.size == 0:
        return 0.0
    return float(pd.Series(diff).max())


def _mean_abs_diff(actual: pd.DataFrame, expected: pd.DataFrame) -> float:
    aligned_actual = actual.sort_index().sort_index(axis=1)
    aligned_expected = expected.sort_index().sort_index(axis=1)
    diff = (aligned_actual - aligned_expected).abs().to_numpy().ravel()
    if diff.size == 0:
        return 0.0
    return float(pd.Series(diff).mean())


def _stable_json(value: dict[str, Any]) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _format_runtime_stats(runs: list[float]) -> RuntimeStats:
    return RuntimeStats(
        runs=[round(value, 6) for value in runs],
        mean_seconds=round(statistics.fmean(runs), 6),
        median_seconds=round(statistics.median(runs), 6),
        min_seconds=round(min(runs), 6),
        max_seconds=round(max(runs), 6),
    )


def _run_with_timing(
    runner: Callable[[ModeName], dict[str, Any]],
    *,
    mode: ModeName,
    repeats: int,
) -> tuple[dict[str, Any], RuntimeStats]:
    if repeats < 1:
        raise ValueError("repeats must be at least 1")

    runtimes: list[float] = []
    baseline_signature: str | None = None
    baseline_metrics: dict[str, Any] | None = None

    for _ in range(repeats):
        start = time.perf_counter()
        metrics = runner(mode)
        duration = time.perf_counter() - start
        runtimes.append(duration)

        signature = _stable_json(metrics)
        if baseline_signature is None:
            baseline_signature = signature
            baseline_metrics = metrics
            continue
        if signature != baseline_signature:
            raise RuntimeError(
                "Benchmark metrics changed across repeated runs. "
                "The harness expects deterministic benchmark outputs for the "
                f"same mode and fixture set. Mode: {mode}"
            )

    assert baseline_metrics is not None
    return baseline_metrics, _format_runtime_stats(runtimes)


def _load_demo_inputs() -> tuple[
    pd.DataFrame,
    dict[str, list[str]],
    dict[str, str],
    dict[str, list[str]],
]:
    phospho_matrix = pd.read_csv(
        EXAMPLE_DATA / "predmat_phospho_matrix.csv",
        index_col=0,
    )
    phospho_matrix.index = phospho_matrix.index.map(str)
    substrate_map = json.loads(
        (EXAMPLE_DATA / "predmat_substrate_map.json").read_text(encoding="utf-8")
    )
    site_sequences = json.loads(
        (EXAMPLE_DATA / "predmat_site_sequences.json").read_text(encoding="utf-8")
    )
    motif_sequences = json.loads(
        (EXAMPLE_DATA / "predmat_motif_sequences.json").read_text(encoding="utf-8")
    )
    return phospho_matrix, substrate_map, site_sequences, motif_sequences


def _run_prediction_ranking_l6(mode: ModeName) -> dict[str, Any]:
    _require_paths(
        [
            R_FIXTURES_L6 / "native_combined_scores.csv",
            R_FIXTURES_L6 / "predMat.csv",
            R_FIXTURES_L6 / "native_prediction_top30.csv",
        ]
    )
    combined_scores = _read_indexed_table(R_FIXTURES_L6 / "native_combined_scores.csv")
    expected_pred = _read_indexed_table(R_FIXTURES_L6 / "predMat.csv")
    expected_top30 = _read_ranked_sites_by_kinase(
        R_FIXTURES_L6 / "native_prediction_top30.csv"
    )

    result = KinasePredictor(svm_mode=mode).predict(
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
    if not expected_kinases:
        raise RuntimeError("The L6 ranking fixture does not contain any kinases")

    common_kinases = [
        kinase for kinase in expected_kinases if kinase in expected_pred.columns
    ]
    if len(common_kinases) != len(expected_kinases):
        raise RuntimeError(
            "The expected ranking fixtures are missing kinases present in the top-30 fixture"
        )

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
            ascending=False,
            method="average",
        )
        actual_ranks = actual_pred.loc[:, kinase].rank(
            ascending=False,
            method="average",
        )
        rank_correlations.append(
            float(expected_ranks.corr(actual_ranks, method="spearman"))
        )

    return {
        "n_kinases": len(common_kinases),
        "mean_spearman": float(pd.Series(rank_correlations).mean()),
        "mean_top10_overlap": float(pd.Series(top10_overlaps).mean()),
        "mean_top20_overlap": float(pd.Series(top20_overlaps).mean()),
        "mean_top30_overlap": float(pd.Series(top30_overlaps).mean()),
        "n_good_top10": int(sum(overlap >= 0.7 for overlap in top10_overlaps)),
    }


def _read_prediction_trace_table(name: str) -> pd.DataFrame:
    return pd.read_csv(R_FIXTURES_L6 / "prediction_trace" / name)


def _run_replayed_trace_l6(mode: ModeName) -> dict[str, Any]:
    _require_paths(
        [
            R_FIXTURES_L6 / "native_combined_scores.csv",
            R_FIXTURES_L6 / "prediction_trace" / "trace_candidates.csv",
            R_FIXTURES_L6 / "prediction_trace" / "trace_initial_negatives.csv",
            R_FIXTURES_L6 / "prediction_trace" / "trace_iteration_samples.csv",
            R_FIXTURES_L6 / "prediction_trace" / "trace_iteration_probabilities.csv",
            R_FIXTURES_L6 / "prediction_trace" / "trace_iteration_decision_values.csv",
            R_FIXTURES_L6 / "prediction_trace" / "trace_final_ensemble_top.csv",
        ]
    )
    combined_scores = _read_indexed_table(R_FIXTURES_L6 / "native_combined_scores.csv")
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
        raise RuntimeError(
            "No replay-aligned trace kinases are available. Regenerate the R "
            "prediction trace fixtures so candidate counts and sampling rows "
            "match the current native candidate set."
        )

    sampling_trace = PredictionSamplingTrace.from_trace_directory(
        R_FIXTURES_L6 / "prediction_trace"
    ).subset_kinases(eligible_kinases)

    with tempfile.TemporaryDirectory(
        prefix="phospy_trace_benchmark_"
    ) as trace_output_dir:
        result = KinasePredictor(svm_mode=mode).predict(
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
        actual_tables = prediction_debug_trace_tables(result)

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
    sample_exact_matches = int(actual_samples.eq(expected_samples).all(axis=1).sum())

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

    final_top_total = int(len(merged_top))
    final_top_site_matches = int((merged_top["site_py"] == merged_top["site_r"]).sum())
    initial_total_rows = int(len(expected_initial))
    sample_total_rows = int(len(expected_samples))

    return {
        "n_kinases": len(eligible_kinases),
        "trace_kinases": ", ".join(eligible_kinases),
        "skipped_trace_kinases": ", ".join(skipped_kinases),
        "initial_exact_matches": initial_exact_matches,
        "initial_total_rows": initial_total_rows,
        "initial_exact_match_rate": (
            float(initial_exact_matches) / float(initial_total_rows)
            if initial_total_rows
            else 1.0
        ),
        "sample_exact_matches": int(
            actual_samples.eq(expected_samples).all(axis=1).sum()
        ),
        "sample_total_rows": sample_total_rows,
        "sample_exact_match_rate": (
            float(sample_exact_matches) / float(sample_total_rows)
            if sample_total_rows
            else 1.0
        ),
        "iteration_prob_class1_corr": float(
            merged_prob["prob_class_1_py"].corr(
                merged_prob["prob_class_1_r"],
                method="pearson",
            )
        ),
        "iteration_prob_class2_corr": float(
            merged_prob["prob_class_2_py"].corr(
                merged_prob["prob_class_2_r"],
                method="pearson",
            )
        ),
        "iteration_decision_class1_corr": float(
            merged_decision["decision_value_class_1_py"].corr(
                merged_decision["decision_value_class_1_r"],
                method="pearson",
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
        "final_top_match_rate": (
            float(final_top_site_matches) / float(final_top_total)
            if final_top_total
            else 1.0
        ),
        "final_top_prob_mae": float(
            (merged_top["prob_class_1_py"] - merged_top["prob_class_1_r"]).abs().mean()
        ),
    }


def _run_public_predmat_workflow(mode: ModeName) -> dict[str, Any]:
    _require_paths(
        [
            EXAMPLE_DATA / "predmat_phospho_matrix.csv",
            EXAMPLE_DATA / "predmat_substrate_map.json",
            EXAMPLE_DATA / "predmat_site_sequences.json",
            EXAMPLE_DATA / "predmat_motif_sequences.json",
            WORKFLOW_FIXTURES / PREDMAT_BENCHMARKS[mode],
        ]
    )
    phospho_matrix, substrate_map, site_sequences, motif_sequences = _load_demo_inputs()
    actual = (
        PredMatWorkflow(flank_size=2, svm_mode=mode)
        .run(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            motif_sequences=motif_sequences,
            min_substrates=2,
            min_motif_size=2,
            ensemble_size=3,
            top=4,
            score_threshold=0.75,
            inclusion=3,
            n_iterations=2,
            random_state=17,
        )
        .pred_mat_result.to_frame(copy=False)
    )
    expected = _read_indexed_table(WORKFLOW_FIXTURES / PREDMAT_BENCHMARKS[mode])

    exact_match = True
    try:
        pdt.assert_frame_equal(actual, expected)
    except AssertionError:
        exact_match = False

    dominant = actual.idxmax(axis=1).to_dict()
    expected_dominant = {
        "SITE_1": "KINASE_A",
        "SITE_2": "KINASE_A",
        "SITE_3": "KINASE_A",
        "SITE_4": "KINASE_A",
        "SITE_5": "KINASE_B",
        "SITE_6": "KINASE_B",
        "SITE_7": "KINASE_B",
        "SITE_8": "KINASE_B",
    }
    return {
        "exact_match": exact_match,
        "rows": int(actual.shape[0]),
        "kinases": int(actual.shape[1]),
        "mean_abs_diff": _mean_abs_diff(actual, expected),
        "max_abs_diff": _max_abs_diff(actual, expected),
        "dominant_assignment_matches_expected": dominant == expected_dominant,
    }


def _run_public_signalome_workflow(mode: ModeName) -> dict[str, Any]:
    _require_paths(
        [
            EXAMPLE_DATA / "predmat_phospho_matrix.csv",
            EXAMPLE_DATA / "predmat_substrate_map.json",
            EXAMPLE_DATA / "predmat_site_sequences.json",
            EXAMPLE_DATA / "predmat_motif_sequences.json",
            WORKFLOW_FIXTURES / SIGNALOME_BENCHMARKS["modules"],
            WORKFLOW_FIXTURES / SIGNALOME_BENCHMARKS["map_modules"],
            WORKFLOW_FIXTURES / SIGNALOME_BENCHMARKS["network_nodes"],
            WORKFLOW_FIXTURES / SIGNALOME_BENCHMARKS["network_edges"],
        ]
    )
    phospho_matrix, substrate_map, site_sequences, motif_sequences = _load_demo_inputs()
    pred_mat_result = PredMatWorkflow(flank_size=2, svm_mode=mode).run(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        site_sequences=site_sequences,
        motif_sequences=motif_sequences,
        min_substrates=2,
        min_motif_size=2,
        ensemble_size=3,
        top=4,
        score_threshold=0.75,
        inclusion=3,
        n_iterations=2,
        random_state=17,
    )
    signalome_result = SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A", "KINASE_B"],
        site_to_protein={
            str(site_id): str(site_id).split(";", 1)[0]
            for site_id in phospho_matrix.index
        },
        signalome_cutoff=0.5,
    )
    map_data = signalome_result.to_map_data()
    network_data = signalome_result.to_network_data()

    actual_modules = signalome_result.modules.to_frame()
    actual_map_modules = map_data.modules()
    actual_nodes = network_data.nodes()
    actual_edges = network_data.edges()

    expected_modules = _read_indexed_table(
        WORKFLOW_FIXTURES / SIGNALOME_BENCHMARKS["modules"]
    )
    expected_map_modules = _read_indexed_table(
        WORKFLOW_FIXTURES / SIGNALOME_BENCHMARKS["map_modules"]
    )
    expected_nodes = _read_indexed_table(
        WORKFLOW_FIXTURES / SIGNALOME_BENCHMARKS["network_nodes"]
    )
    expected_edges = _read_unindexed_table(
        WORKFLOW_FIXTURES / SIGNALOME_BENCHMARKS["network_edges"]
    )

    expected_modules.index.name = actual_modules.index.name
    expected_modules.columns.name = actual_modules.columns.name
    expected_map_modules.index.name = actual_map_modules.index.name
    expected_nodes.index.name = actual_nodes.index.name

    modules_exact_match = True
    map_modules_exact_match = True
    nodes_exact_match = True
    edges_exact_match = True
    try:
        pdt.assert_frame_equal(actual_modules, expected_modules)
    except AssertionError:
        modules_exact_match = False
    try:
        pdt.assert_frame_equal(actual_map_modules, expected_map_modules)
    except AssertionError:
        map_modules_exact_match = False
    try:
        pdt.assert_frame_equal(actual_nodes, expected_nodes)
    except AssertionError:
        nodes_exact_match = False
    try:
        pdt.assert_frame_equal(actual_edges, expected_edges, check_dtype=False)
    except AssertionError:
        edges_exact_match = False

    return {
        "modules_exact_match": modules_exact_match,
        "map_modules_exact_match": map_modules_exact_match,
        "network_nodes_exact_match": nodes_exact_match,
        "network_edges_exact_match": edges_exact_match,
        "overall_exact_match": (
            modules_exact_match
            and map_modules_exact_match
            and nodes_exact_match
            and edges_exact_match
        ),
        "n_modules": int(actual_modules.shape[0]),
        "n_network_nodes": int(actual_nodes.shape[0]),
        "n_network_edges": int(actual_edges.shape[0]),
    }


def _evaluate_prediction_ranking_thresholds(metrics: dict[str, Any]) -> dict[str, bool]:
    return {
        "mean_spearman_ge_0_96": bool(metrics["mean_spearman"] >= 0.96),
        "mean_top20_overlap_ge_0_85": bool(metrics["mean_top20_overlap"] >= 0.85),
        "mean_top30_overlap_ge_0_88": bool(metrics["mean_top30_overlap"] >= 0.88),
        "n_good_top10_ge_20": bool(metrics["n_good_top10"] >= 20),
    }


def _evaluate_trace_replay_thresholds(metrics: dict[str, Any]) -> dict[str, bool]:
    return {
        "initial_exact_match_rate_eq_1": bool(
            math.isclose(metrics["initial_exact_match_rate"], 1.0)
        ),
        "sample_exact_match_rate_eq_1": bool(
            math.isclose(metrics["sample_exact_match_rate"], 1.0)
        ),
        "iteration_decision_class1_corr_ge_0_999999": bool(
            metrics["iteration_decision_class1_corr"] >= 0.999999
        ),
        "iteration_decision_mae_le_1e_minus_12": bool(
            metrics["iteration_decision_mae"] <= 1e-12
        ),
        "iteration_prob_class1_corr_ge_0_998": bool(
            metrics["iteration_prob_class1_corr"] >= 0.998
        ),
        "iteration_prob_mae_le_0_015": bool(metrics["iteration_prob_mae"] <= 0.015),
        "final_top_match_rate_eq_1": bool(
            math.isclose(metrics["final_top_match_rate"], 1.0)
        ),
    }


def _evaluate_public_predmat_thresholds(metrics: dict[str, Any]) -> dict[str, bool]:
    return {
        "exact_match": bool(metrics["exact_match"]),
        "dominant_assignment_matches_expected": bool(
            metrics["dominant_assignment_matches_expected"]
        ),
    }


def _evaluate_public_signalome_thresholds(metrics: dict[str, Any]) -> dict[str, bool]:
    return {
        "overall_exact_match": bool(metrics["overall_exact_match"]),
        "modules_exact_match": bool(metrics["modules_exact_match"]),
        "map_modules_exact_match": bool(metrics["map_modules_exact_match"]),
        "network_nodes_exact_match": bool(metrics["network_nodes_exact_match"]),
        "network_edges_exact_match": bool(metrics["network_edges_exact_match"]),
    }


def _numeric_delta(default_value: Any, r_parity_value: Any) -> float | None:
    if isinstance(default_value, bool) or isinstance(r_parity_value, bool):
        return None
    if isinstance(default_value, (int, float)) and isinstance(
        r_parity_value, (int, float)
    ):
        return float(r_parity_value) - float(default_value)
    return None


def _build_standard_comparison(
    results: dict[ModeName, ModeBenchmarkResult],
) -> dict[str, Any]:
    default_result = results["default"]
    parity_result = results["r_parity"]

    deltas: dict[str, float] = {}
    for key, default_value in default_result.metrics.items():
        if key not in parity_result.metrics:
            continue
        delta = _numeric_delta(default_value, parity_result.metrics[key])
        if delta is not None:
            deltas[key] = delta

    runtime_ratio = (
        parity_result.runtime.median_seconds / default_result.runtime.median_seconds
        if default_result.runtime.median_seconds
        else float("nan")
    )

    return {
        "metric_deltas_r_parity_minus_default": deltas,
        "median_runtime_seconds": {
            "default": default_result.runtime.median_seconds,
            "r_parity": parity_result.runtime.median_seconds,
        },
        "r_parity_over_default_runtime_ratio": round(runtime_ratio, 6),
    }


def _build_prediction_comparison(
    results: dict[ModeName, ModeBenchmarkResult],
) -> dict[str, Any]:
    comparison = _build_standard_comparison(results)
    default_metrics = results["default"].metrics
    r_parity_metrics = results["r_parity"].metrics
    comparison["protected_mode_checks"] = {
        "same_n_kinases": default_metrics["n_kinases"] == r_parity_metrics["n_kinases"],
        "r_parity_mean_spearman_ge_default": (
            r_parity_metrics["mean_spearman"] >= default_metrics["mean_spearman"]
        ),
        "r_parity_mean_top10_overlap_ge_default": (
            r_parity_metrics["mean_top10_overlap"]
            >= default_metrics["mean_top10_overlap"]
        ),
        "r_parity_mean_top10_overlap_ge_0_82": (
            r_parity_metrics["mean_top10_overlap"] >= 0.82
        ),
        "r_parity_mean_top20_overlap_ge_default": (
            r_parity_metrics["mean_top20_overlap"]
            >= default_metrics["mean_top20_overlap"]
        ),
        "r_parity_mean_top30_overlap_ge_default": (
            r_parity_metrics["mean_top30_overlap"]
            >= default_metrics["mean_top30_overlap"]
        ),
    }
    return comparison


def _build_trace_comparison(
    results: dict[ModeName, ModeBenchmarkResult],
) -> dict[str, Any]:
    comparison = _build_standard_comparison(results)
    default_metrics = results["default"].metrics
    r_parity_metrics = results["r_parity"].metrics
    comparison["protected_mode_checks"] = {
        "same_initial_total_rows": (
            default_metrics["initial_total_rows"]
            == r_parity_metrics["initial_total_rows"]
        ),
        "same_sample_total_rows": (
            default_metrics["sample_total_rows"]
            == r_parity_metrics["sample_total_rows"]
        ),
        "r_parity_iteration_decision_corr_ge_default": (
            r_parity_metrics["iteration_decision_class1_corr"]
            >= default_metrics["iteration_decision_class1_corr"]
        ),
        "r_parity_final_top_matches_ge_default": (
            r_parity_metrics["final_top_site_matches"]
            >= default_metrics["final_top_site_matches"]
        ),
    }
    return comparison


def _build_workflow_comparison(
    results: dict[ModeName, ModeBenchmarkResult],
) -> dict[str, Any]:
    comparison = _build_standard_comparison(results)
    comparison["protected_mode_checks"] = {
        "default_thresholds_all_pass": all(results["default"].thresholds.values()),
        "r_parity_thresholds_all_pass": all(results["r_parity"].thresholds.values()),
    }
    return comparison


def _benchmark_cases() -> list[BenchmarkCase]:
    return [
        BenchmarkCase(
            name="prediction_ranking_l6",
            fixture_family="tests/fixtures/r_reference_l6",
            protected_seam=(
                "Prediction ranking agreement against the main L6 parity dataset"
            ),
            metric_class=(
                "Mean Spearman rank agreement, top-N overlap, and top-10 support counts"
            ),
            runner=_run_prediction_ranking_l6,
            threshold_evaluator=_evaluate_prediction_ranking_thresholds,
            comparison_builder=_build_prediction_comparison,
        ),
        BenchmarkCase(
            name="replayed_prediction_trace_l6",
            fixture_family="tests/fixtures/r_reference_l6",
            protected_seam=(
                "Replay fidelity against the committed L6 R sampling trace"
            ),
            metric_class=(
                "Exact trace-row agreement, probability correlation, mean absolute error, and final top-site agreement"
            ),
            runner=_run_replayed_trace_l6,
            threshold_evaluator=_evaluate_trace_replay_thresholds,
            comparison_builder=_build_trace_comparison,
        ),
        BenchmarkCase(
            name="public_predmat_workflow",
            fixture_family="tests/fixtures/public_workflow_reference",
            protected_seam="Public PredMatWorkflow demo benchmark outputs",
            metric_class="Exact benchmark equality plus dominant-kinase assignment checks",
            runner=_run_public_predmat_workflow,
            threshold_evaluator=_evaluate_public_predmat_thresholds,
            comparison_builder=_build_workflow_comparison,
        ),
        BenchmarkCase(
            name="public_signalome_workflow",
            fixture_family="tests/fixtures/public_workflow_reference",
            protected_seam="Public SignalomeWorkflow demo benchmark outputs",
            metric_class=(
                "Exact benchmark equality across modules, map data, and network outputs"
            ),
            runner=_run_public_signalome_workflow,
            threshold_evaluator=_evaluate_public_signalome_thresholds,
            comparison_builder=_build_workflow_comparison,
        ),
    ]


def _run_case(case: BenchmarkCase, *, repeats: int) -> BenchmarkCaseResult:
    mode_results: dict[ModeName, ModeBenchmarkResult] = {}
    print(f"Running benchmark: {case.name}", file=sys.stderr, flush=True)
    for mode in MODES:
        print(
            f"  - mode={mode}, repeats={repeats}",
            file=sys.stderr,
            flush=True,
        )
        metrics, runtime = _run_with_timing(case.runner, mode=mode, repeats=repeats)
        thresholds = case.threshold_evaluator(metrics)
        mode_results[mode] = ModeBenchmarkResult(
            mode=mode,
            metrics=metrics,
            thresholds=thresholds,
            runtime=runtime,
        )
    comparisons = case.comparison_builder(mode_results)
    return BenchmarkCaseResult(
        benchmark=case.name,
        fixture_family=case.fixture_family,
        protected_seam=case.protected_seam,
        metric_class=case.metric_class,
        modes=mode_results,
        comparisons=comparisons,
    )


def _build_report(
    case_results: list[BenchmarkCaseResult], *, repeats: int
) -> dict[str, Any]:
    return {
        "report": "prediction_mode_comparison",
        "repeats_per_mode": repeats,
        "benchmarks": [case_result.to_dict() for case_result in case_results],
        "summary": {
            "benchmarks": [case_result.benchmark for case_result in case_results],
            "all_default_thresholds_pass": all(
                all(case_result.modes["default"].thresholds.values())
                for case_result in case_results
            ),
            "all_r_parity_thresholds_pass": all(
                all(case_result.modes["r_parity"].thresholds.values())
                for case_result in case_results
            ),
        },
    }


def _format_value(value: Any) -> str:
    if isinstance(value, bool):
        return "pass" if value else "fail"
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        return f"{value:.6f}"
    return str(value)


def _case_verdict(case_result: BenchmarkCaseResult) -> str:
    default_pass = all(case_result.modes["default"].thresholds.values())
    r_parity_pass = all(case_result.modes["r_parity"].thresholds.values())
    runtime_ratio = case_result.comparisons["r_parity_over_default_runtime_ratio"]
    return (
        f"default={'pass' if default_pass else 'fail'}, "
        f"r_parity={'pass' if r_parity_pass else 'fail'}, "
        f"runtime_ratio={runtime_ratio:.3f}"
    )


def _render_summary_table(case_results: list[BenchmarkCaseResult]) -> list[str]:
    lines = [
        "| Benchmark | Fixture family | Default median (s) | r_parity median (s) | Runtime ratio | Verdict |",
        "| --- | --- | ---: | ---: | ---: | --- |",
    ]
    for case_result in case_results:
        default_runtime = case_result.modes["default"].runtime.median_seconds
        r_parity_runtime = case_result.modes["r_parity"].runtime.median_seconds
        runtime_ratio = case_result.comparisons["r_parity_over_default_runtime_ratio"]
        lines.append(
            "| "
            f"`{case_result.benchmark}` | `{case_result.fixture_family}` | "
            f"{default_runtime:.6f} | {r_parity_runtime:.6f} | {runtime_ratio:.6f} | "
            f"{_case_verdict(case_result)} |"
        )
    return lines


def _render_metrics_table(case_result: BenchmarkCaseResult) -> list[str]:
    default_metrics = case_result.modes["default"].metrics
    r_parity_metrics = case_result.modes["r_parity"].metrics
    lines = [
        "| Metric | default | r_parity | Delta (r_parity - default) |",
        "| --- | ---: | ---: | ---: |",
    ]
    all_keys = list(default_metrics.keys())
    for key in all_keys:
        default_value = default_metrics[key]
        r_parity_value = r_parity_metrics.get(key)
        delta = _numeric_delta(default_value, r_parity_value)
        delta_text = _format_value(delta) if delta is not None else "n/a"
        lines.append(
            f"| `{key}` | {_format_value(default_value)} | {_format_value(r_parity_value)} | {delta_text} |"
        )
    return lines


def _render_threshold_table(case_result: BenchmarkCaseResult) -> list[str]:
    lines = [
        "| Threshold | default | r_parity |",
        "| --- | --- | --- |",
    ]
    keys = list(case_result.modes["default"].thresholds.keys())
    for key in keys:
        lines.append(
            f"| `{key}` | {_format_value(case_result.modes['default'].thresholds[key])} | "
            f"{_format_value(case_result.modes['r_parity'].thresholds[key])} |"
        )
    return lines


def _render_comparison_list(case_result: BenchmarkCaseResult) -> list[str]:
    lines = ["Protected comparison checks:"]
    checks = case_result.comparisons.get("protected_mode_checks", {})
    for key, value in checks.items():
        lines.append(f"- `{key}`: {_format_value(value)}")
    lines.append(
        "- `default_median_runtime_seconds`: "
        f"{case_result.modes['default'].runtime.median_seconds:.6f}"
    )
    lines.append(
        "- `r_parity_median_runtime_seconds`: "
        f"{case_result.modes['r_parity'].runtime.median_seconds:.6f}"
    )
    lines.append(
        "- `r_parity_over_default_runtime_ratio`: "
        f"{case_result.comparisons['r_parity_over_default_runtime_ratio']:.6f}"
    )
    return lines


def render_markdown(report: dict[str, Any]) -> str:
    case_results = [
        BenchmarkCaseResult(
            benchmark=case["benchmark"],
            fixture_family=case["fixture_family"],
            protected_seam=case["protected_seam"],
            metric_class=case["metric_class"],
            modes={
                mode: ModeBenchmarkResult(
                    mode=mode,  # type: ignore[arg-type]
                    metrics=payload["metrics"],
                    thresholds=payload["thresholds"],
                    runtime=RuntimeStats(
                        runs=payload["runtime"]["runs_seconds"],
                        mean_seconds=payload["runtime"]["mean_seconds"],
                        median_seconds=payload["runtime"]["median_seconds"],
                        min_seconds=payload["runtime"]["min_seconds"],
                        max_seconds=payload["runtime"]["max_seconds"],
                    ),
                )
                for mode, payload in case["modes"].items()
            },
            comparisons=case["comparisons"],
        )
        for case in report["benchmarks"]
    ]

    lines = [
        "# Prediction Mode Benchmark Report",
        "",
        'This report compares `svm_mode="default"` and `svm_mode="r_parity"` '
        "on the selected parity fixture families.",
        "",
        f"- repeats per mode: {report['repeats_per_mode']}",
        f"- all default thresholds pass: {_format_value(report['summary']['all_default_thresholds_pass'])}",
        f"- all r_parity thresholds pass: {_format_value(report['summary']['all_r_parity_thresholds_pass'])}",
        "",
        "## Benchmark Summary",
        "",
        *(_render_summary_table(case_results)),
    ]

    for case_result in case_results:
        lines.extend(
            [
                "",
                f"## `{case_result.benchmark}`",
                "",
                f"- fixture family: `{case_result.fixture_family}`",
                f"- protected seam: {case_result.protected_seam}",
                f"- metric class: {case_result.metric_class}",
                "",
                "### Metrics",
                "",
                *(_render_metrics_table(case_result)),
                "",
                "### Threshold checks",
                "",
                *(_render_threshold_table(case_result)),
                "",
                "### Mode comparison",
                "",
                *(_render_comparison_list(case_result)),
            ]
        )

    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark svm_mode='default' versus svm_mode='r_parity' on the "
            "selected parity fixture families."
        )
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="Number of repeated timed runs per mode and benchmark case (default: 3)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=(
            "Directory for the generated JSON and Markdown reports "
            f"(default: {DEFAULT_OUTPUT_DIR.relative_to(REPO_ROOT)})"
        ),
    )
    parser.add_argument(
        "--stdout-only",
        action="store_true",
        help="Print the Markdown report to stdout without writing files.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    case_results = [
        _run_case(case, repeats=args.repeats) for case in _benchmark_cases()
    ]
    report = _build_report(case_results, repeats=args.repeats)
    markdown = render_markdown(report)

    if args.stdout_only:
        print(markdown)
        return 0

    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "compare_prediction_modes.json"
    markdown_path = output_dir / "compare_prediction_modes.md"

    json_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    markdown_path.write_text(markdown, encoding="utf-8")

    print(f"Wrote JSON report: {json_path.relative_to(REPO_ROOT)}")
    print(f"Wrote Markdown report: {markdown_path.relative_to(REPO_ROOT)}")
    print()
    print(markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
