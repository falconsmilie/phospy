#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd


def _load_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _time_call(repeats: int, func, *args, **kwargs) -> tuple[Any, float]:
    durations: list[float] = []
    result: Any = None
    for _ in range(repeats):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        durations.append(time.perf_counter() - start)
    assert result is not None
    return result, sum(durations) / len(durations)


def _frame_matches(actual: pd.DataFrame, expected: pd.DataFrame) -> bool:
    try:
        pd.testing.assert_frame_equal(actual, expected, check_dtype=False)
    except AssertionError:
        return False
    return True


def _build_markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Prediction mode comparison",
        "",
        "This report compares the supported public prediction presets on the protected parity evidence.",
        "",
        "## Threshold summary",
        "",
        f"- all default thresholds pass: `{report['thresholds']['all_default_thresholds_pass']}`",
        f"- all r_parity thresholds pass: `{report['thresholds']['all_r_parity_thresholds_pass']}`",
        "",
        "## Ranking parity",
        "",
        f"- default mean Spearman: `{report['ranking']['default']['metrics']['mean_spearman']:.6f}`",
        f"- r_parity mean Spearman: `{report['ranking']['r_parity']['metrics']['mean_spearman']:.6f}`",
        f"- default mean top-10 overlap: `{report['ranking']['default']['metrics']['mean_top10_overlap']:.6f}`",
        f"- r_parity mean top-10 overlap: `{report['ranking']['r_parity']['metrics']['mean_top10_overlap']:.6f}`",
        f"- default mean runtime (s): `{report['ranking']['default']['mean_runtime_seconds']:.6f}`",
        f"- r_parity mean runtime (s): `{report['ranking']['r_parity']['mean_runtime_seconds']:.6f}`",
        "",
        "## Replay-trace parity",
        "",
        f"- default iteration prob corr: `{report['replay']['default']['metrics']['iteration_prob_class1_corr']:.6f}`",
        f"- r_parity iteration prob corr: `{report['replay']['r_parity']['metrics']['iteration_prob_class1_corr']:.6f}`",
        f"- default final top-site matches: `{report['replay']['default']['metrics']['final_top_site_matches']}/{report['replay']['default']['metrics']['final_top_total']}`",
        f"- r_parity final top-site matches: `{report['replay']['r_parity']['metrics']['final_top_site_matches']}/{report['replay']['r_parity']['metrics']['final_top_total']}`",
        f"- default mean runtime (s): `{report['replay']['default']['mean_runtime_seconds']:.6f}`",
        f"- r_parity mean runtime (s): `{report['replay']['r_parity']['mean_runtime_seconds']:.6f}`",
        "",
        "## Public workflow fixtures",
        "",
        f"- predMat default benchmark match: `{report['public_workflows']['predmat_default_matches_fixture']}`",
        f"- predMat r_parity benchmark match: `{report['public_workflows']['predmat_r_parity_matches_fixture']}`",
        f"- signalome default benchmark match: `{report['public_workflows']['signalome_default_matches_fixture']}`",
        f"- signalome r_parity benchmark match: `{report['public_workflows']['signalome_r_parity_matches_fixture']}`",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Compare default and r_parity on the protected parity datasets and write a reviewable report."
        )
    )
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("benchmarks/reports/latest"),
        help="Directory for compare_prediction_modes.json and compare_prediction_modes.md",
    )
    parser.add_argument(
        "--stdout-only",
        action="store_true",
        help="Print the JSON report to stdout instead of writing files.",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root / "src"))

    parity_module = _load_module(
        repo_root / "tests" / "test_parity-with_metrics.py",
        "phospy_test_parity_with_metrics",
    )
    e2e_module = _load_module(
        repo_root / "tests" / "test_end_to_end_parity.py",
        "phospy_test_end_to_end_parity",
    )

    ranking_default_metrics, ranking_default_runtime = _time_call(
        args.repeats,
        parity_module._prediction_parity_metrics,
        svm_mode="default",
    )
    ranking_r_parity_metrics, ranking_r_parity_runtime = _time_call(
        args.repeats,
        parity_module._prediction_parity_metrics,
        svm_mode="r_parity",
    )
    replay_default_metrics, replay_default_runtime = _time_call(
        args.repeats,
        parity_module._replayed_prediction_trace_metrics,
        svm_mode="default",
    )
    replay_r_parity_metrics, replay_r_parity_runtime = _time_call(
        args.repeats,
        parity_module._replayed_prediction_trace_metrics,
        svm_mode="r_parity",
    )

    predmat_default, predmat_default_runtime = _time_call(
        args.repeats,
        e2e_module._run_public_predmat_workflow,
        svm_mode="default",
    )
    predmat_r_parity, predmat_r_parity_runtime = _time_call(
        args.repeats,
        e2e_module._run_public_predmat_workflow,
        svm_mode="r_parity",
    )
    signalome_default, signalome_default_runtime = _time_call(
        args.repeats,
        e2e_module._run_public_signalome_workflow,
        svm_mode="default",
    )
    signalome_r_parity, signalome_r_parity_runtime = _time_call(
        args.repeats,
        e2e_module._run_public_signalome_workflow,
        svm_mode="r_parity",
    )

    expected_predmat_default = e2e_module._read_indexed_fixture(
        e2e_module.PREDMAT_BENCHMARKS["default"]
    )
    expected_predmat_r_parity = e2e_module._read_indexed_fixture(
        e2e_module.PREDMAT_BENCHMARKS["r_parity"]
    )
    expected_modules = e2e_module._read_indexed_fixture(
        e2e_module.SIGNALOME_BENCHMARKS["modules"]
    )
    expected_map_modules = e2e_module._read_indexed_fixture(
        e2e_module.SIGNALOME_BENCHMARKS["map_modules"]
    )
    expected_nodes = e2e_module._read_indexed_fixture(
        e2e_module.SIGNALOME_BENCHMARKS["network_nodes"]
    )
    expected_edges = e2e_module._read_unindexed_fixture(
        e2e_module.SIGNALOME_BENCHMARKS["network_edges"]
    )

    (
        signalome_default_modules,
        signalome_default_map,
        signalome_default_nodes,
        signalome_default_edges,
    ) = signalome_default
    (
        signalome_r_parity_modules,
        signalome_r_parity_map,
        signalome_r_parity_nodes,
        signalome_r_parity_edges,
    ) = signalome_r_parity

    default_thresholds = parity_module.DEFAULT_RANKING_RELEASE_THRESHOLDS
    r_parity_ranking_thresholds = parity_module.R_PARITY_RANKING_RELEASE_THRESHOLDS
    r_parity_replay_thresholds = parity_module.R_PARITY_REPLAY_RELEASE_THRESHOLDS

    all_default_thresholds_pass = (
        ranking_default_metrics["mean_spearman"] >= default_thresholds["mean_spearman"]
        and ranking_default_metrics["mean_top20_overlap"]
        >= default_thresholds["mean_top20_overlap"]
        and ranking_default_metrics["mean_top30_overlap"]
        >= default_thresholds["mean_top30_overlap"]
        and ranking_default_metrics["n_good_top10"]
        >= default_thresholds["n_good_top10"]
    )
    all_r_parity_thresholds_pass = (
        all(
            ranking_r_parity_metrics[metric_name]
            >= ranking_default_metrics[metric_name]
            for metric_name in r_parity_ranking_thresholds[
                "must_match_or_exceed_default"
            ]
        )
        and ranking_r_parity_metrics["mean_top10_overlap"]
        >= r_parity_ranking_thresholds["mean_top10_overlap_floor"]
        and replay_r_parity_metrics["initial_exact_matches"]
        == replay_r_parity_metrics["initial_total_rows"]
        and replay_r_parity_metrics["sample_exact_matches"]
        == replay_r_parity_metrics["sample_total_rows"]
        and replay_r_parity_metrics["iteration_decision_class1_corr"]
        >= r_parity_replay_thresholds["iteration_decision_class1_corr"]
        and replay_r_parity_metrics["iteration_decision_mae"]
        <= r_parity_replay_thresholds["iteration_decision_mae"]
        and replay_r_parity_metrics["iteration_prob_class1_corr"]
        >= r_parity_replay_thresholds["iteration_prob_class1_corr"]
        and replay_r_parity_metrics["iteration_prob_mae"]
        <= r_parity_replay_thresholds["iteration_prob_mae"]
        and replay_r_parity_metrics["final_top_site_matches"]
        == replay_r_parity_metrics["final_top_total"]
    )

    report = {
        "repeats": args.repeats,
        "thresholds": {
            "all_default_thresholds_pass": all_default_thresholds_pass,
            "all_r_parity_thresholds_pass": all_r_parity_thresholds_pass,
            "default_ranking": default_thresholds,
            "r_parity_ranking": {
                "mean_top10_overlap_floor": r_parity_ranking_thresholds[
                    "mean_top10_overlap_floor"
                ],
                "must_match_or_exceed_default": list(
                    r_parity_ranking_thresholds["must_match_or_exceed_default"]
                ),
            },
            "r_parity_replay": r_parity_replay_thresholds,
        },
        "ranking": {
            "default": {
                "metrics": ranking_default_metrics,
                "mean_runtime_seconds": ranking_default_runtime,
            },
            "r_parity": {
                "metrics": ranking_r_parity_metrics,
                "mean_runtime_seconds": ranking_r_parity_runtime,
            },
        },
        "replay": {
            "default": {
                "metrics": replay_default_metrics,
                "mean_runtime_seconds": replay_default_runtime,
            },
            "r_parity": {
                "metrics": replay_r_parity_metrics,
                "mean_runtime_seconds": replay_r_parity_runtime,
            },
        },
        "public_workflows": {
            "predmat_default_matches_fixture": _frame_matches(
                predmat_default, expected_predmat_default
            ),
            "predmat_r_parity_matches_fixture": _frame_matches(
                predmat_r_parity, expected_predmat_r_parity
            ),
            "signalome_default_matches_fixture": all(
                (
                    _frame_matches(signalome_default_modules, expected_modules),
                    _frame_matches(signalome_default_map, expected_map_modules),
                    _frame_matches(signalome_default_nodes, expected_nodes),
                    _frame_matches(signalome_default_edges, expected_edges),
                )
            ),
            "signalome_r_parity_matches_fixture": all(
                (
                    _frame_matches(signalome_r_parity_modules, expected_modules),
                    _frame_matches(signalome_r_parity_map, expected_map_modules),
                    _frame_matches(signalome_r_parity_nodes, expected_nodes),
                    _frame_matches(signalome_r_parity_edges, expected_edges),
                )
            ),
            "predmat_default_mean_runtime_seconds": predmat_default_runtime,
            "predmat_r_parity_mean_runtime_seconds": predmat_r_parity_runtime,
            "signalome_default_mean_runtime_seconds": signalome_default_runtime,
            "signalome_r_parity_mean_runtime_seconds": signalome_r_parity_runtime,
        },
    }

    if args.stdout_only:
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "compare_prediction_modes.json"
    markdown_path = output_dir / "compare_prediction_modes.md"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    markdown_path.write_text(_build_markdown(report), encoding="utf-8")
    print(f"Wrote {json_path}")
    print(f"Wrote {markdown_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
