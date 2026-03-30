#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from phospy.prediction import PredictionSamplingTrace


def ensure_src_on_path() -> None:
    root = Path(__file__).resolve().parents[1]
    src = root / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export Python prediction debug traces in a CSV layout comparable to the R trace fixtures."
    )
    parser.add_argument(
        "--combined-scores",
        default="tests/fixtures/r_reference_l6/native_combined_scores.csv",
        help="Path to the combined score matrix CSV.",
    )
    parser.add_argument(
        "--outdir",
        default="tests/fixtures/python_reference_l6/prediction_trace",
        help="Directory where Python trace CSVs will be written.",
    )
    parser.add_argument(
        "--trace-kinases",
        default="PRKAA1,MAPK1",
        help="Comma-separated list of kinases to trace.",
    )
    parser.add_argument("--top", type=int, default=30)
    parser.add_argument("--score-threshold", type=float, default=0.6)
    parser.add_argument("--inclusion", type=int, default=5)
    parser.add_argument("--ensemble-size", type=int, default=10)
    parser.add_argument("--n-iterations", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=1)
    parser.add_argument("--debug-top-n", type=int, default=10)
    parser.add_argument(
        "--svm-mode",
        choices=["default", "r_parity"],
        default="r_parity",
        help="SVM configuration mode to use for Python trace export.",
    )
    parser.add_argument(
        "--sampling-trace-dir",
        default=None,
        help=(
            "Optional directory containing R-exported trace_initial_negatives.csv "
            "and trace_iteration_samples.csv. When provided, Python replays "
            "those sampled rows so the remaining delta is model-only."
        ),
    )
    return parser.parse_args()


def parse_csv_values(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def build_trace_candidates(
    combined_scores: pd.DataFrame,
    trace_kinases: list[str],
    top: int,
    score_threshold: float,
    inclusion: int,
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    rows: list[dict[str, object]] = []
    substrate_list: dict[str, list[str]] = {}
    for kinase in trace_kinases:
        if kinase not in combined_scores.columns:
            continue
        ordered = combined_scores.loc[:, kinase].sort_values(
            ascending=False, kind="mergesort"
        )
        top_sites = ordered.head(top)
        selected_sites = top_sites.loc[top_sites > score_threshold].index.tolist()
        if len(selected_sites) >= inclusion:
            substrate_list[kinase] = selected_sites
        for rank, (site, score) in enumerate(ordered.items(), start=1):
            rows.append(
                {
                    "kinase": kinase,
                    "rank": rank,
                    "site": site,
                    "combined_score": float(score),
                    "within_top": rank <= top,
                    "above_threshold": float(score) > score_threshold,
                    "selected_candidate": site in selected_sites,
                }
            )
    return pd.DataFrame(rows), substrate_list


def resolve_replay_aligned_trace_kinases(
    *,
    combined_scores: pd.DataFrame,
    requested_trace_kinases: list[str],
    sampling_trace_dir: str | Path,
    top: int,
    score_threshold: float,
    inclusion: int,
) -> tuple[list[str], list[str], PredictionSamplingTrace]:
    from phospy.prediction import PredictionSamplingTrace

    trace_dir = Path(sampling_trace_dir)
    trace_candidates_path = trace_dir / "trace_candidates.csv"
    trace_initial_path = trace_dir / "trace_initial_negatives.csv"
    trace_samples_path = trace_dir / "trace_iteration_samples.csv"

    if not trace_candidates_path.exists():
        msg = (
            "sampling trace replay requires trace_candidates.csv so the exporter "
            "can verify candidate-set alignment before applying the override"
        )
        raise ValueError(msg)

    _, substrate_list = build_trace_candidates(
        combined_scores,
        requested_trace_kinases,
        top=top,
        score_threshold=score_threshold,
        inclusion=inclusion,
    )
    trace_candidates = pd.read_csv(trace_candidates_path)
    trace_initial = (
        pd.read_csv(trace_initial_path)
        if trace_initial_path.exists()
        else pd.DataFrame(columns=["kinase", "ensemble", "draw", "site"])
    )
    trace_samples = (
        pd.read_csv(trace_samples_path)
        if trace_samples_path.exists()
        else pd.DataFrame(
            columns=["kinase", "ensemble", "iteration", "class_label", "draw", "site"]
        )
    )

    eligible_kinases: list[str] = []
    skipped_kinases: list[str] = []
    for kinase in requested_trace_kinases:
        current_count = len(substrate_list.get(kinase, []))
        trace_count = int(
            trace_candidates.loc[
                (trace_candidates["kinase"].astype(str) == kinase)
                & (trace_candidates["selected_candidate"]),
                "site",
            ].shape[0]
        )
        initial_counts = (
            trace_initial.loc[trace_initial["kinase"].astype(str) == kinase]
            .groupby("ensemble")
            .size()
        )
        sample_counts = (
            trace_samples.loc[trace_samples["kinase"].astype(str) == kinase]
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
            eligible_kinases.append(kinase)
        else:
            skipped_kinases.append(kinase)

    if not eligible_kinases:
        msg = (
            "No replay-aligned trace kinases are available in the requested set. "
            "Regenerate the R prediction trace fixtures so candidate counts and "
            "sampling rows match the current native candidate set."
        )
        raise ValueError(msg)

    sampling_trace = PredictionSamplingTrace.from_trace_directory(
        trace_dir
    ).subset_kinases(eligible_kinases)
    return eligible_kinases, skipped_kinases, sampling_trace


def export_traces(
    result,
    combined_scores: pd.DataFrame,
    outdir: Path,
    *,
    top: int,
    score_threshold: float,
    inclusion: int,
    trace_kinases: list[str],
    svm_mode: str,
    sampling_trace_dir: str | None,
    skipped_trace_kinases: list[str] | None = None,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    trace_candidates, _ = build_trace_candidates(
        combined_scores,
        trace_kinases,
        top=top,
        score_threshold=score_threshold,
        inclusion=inclusion,
    )
    trace_candidates.to_csv(outdir / "trace_candidates.csv", index=False)

    from phospy.prediction import prediction_debug_trace_tables

    for name, table in prediction_debug_trace_tables(result).items():
        table.to_csv(outdir / f"{name}.csv", index=False)

    readme = [
        "# Python prediction trace fixtures",
        "",
        "These files are generated from the native Python predictor with debug tracing enabled.",
        "They are intended for direct comparison with the R trace fixtures generated from the PhosR L6 example path.",
        "",
        f"Trace kinases: {', '.join(trace_kinases)}",
        f"SVM mode: {svm_mode}",
        f"Sampling trace override: {sampling_trace_dir or 'none'}",
        f"Skipped trace kinases: {', '.join(skipped_trace_kinases or []) or 'none'}",
        "",
        "Files:",
        "- trace_candidates.csv: ranked combined-score candidates for the traced kinases",
        "- trace_initial_negatives.csv: initial negative draw for each ensemble member",
        "- trace_iteration_probabilities.csv: per-iteration class probabilities on the base train set",
        "- trace_iteration_probability_parameters.csv: per-iteration libsvm probability-calibration parameters",
        "- trace_iteration_decision_values.csv: per-iteration binary decision values aligned to class 1",
        "- trace_iteration_resampling_weights.csv: per-iteration class-specific resampling weights",
        "- trace_iteration_samples.csv: resampled site identities for each iteration and class",
        "- trace_final_ensemble_predictions.csv: final per-ensemble prediction probabilities for all sites",
        "- trace_final_ensemble_decision_values.csv: final per-ensemble binary decision values aligned to class 1",
        "- trace_final_ensemble_top.csv: final per-ensemble top-ranked sites",
    ]
    (outdir / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")


def main() -> None:

    from phospy.prediction import KinasePredictor

    args = parse_args()
    combined_scores = pd.read_csv(args.combined_scores, index_col=0)
    trace_kinases = parse_csv_values(args.trace_kinases)

    sampling_trace = None
    skipped_trace_kinases: list[str] = []
    if args.sampling_trace_dir:
        trace_kinases, skipped_trace_kinases, sampling_trace = (
            resolve_replay_aligned_trace_kinases(
                combined_scores=combined_scores,
                requested_trace_kinases=trace_kinases,
                sampling_trace_dir=args.sampling_trace_dir,
                top=args.top,
                score_threshold=args.score_threshold,
                inclusion=args.inclusion,
            )
        )

    predictor = KinasePredictor()
    result = predictor.predict(
        combined_scores=combined_scores,
        ensemble_size=args.ensemble_size,
        top=args.top,
        score_threshold=args.score_threshold,
        inclusion=args.inclusion,
        n_iterations=args.n_iterations,
        random_state=args.random_state,
        capture_debug_trace=True,
        debug_kinases=trace_kinases,
        debug_top_n=args.debug_top_n,
        svm_mode=args.svm_mode,
        sampling_trace=sampling_trace,
        trace_level="full",
        trace_sink=Path(args.outdir),
    )

    export_traces(
        result,
        combined_scores,
        Path(args.outdir),
        top=args.top,
        score_threshold=args.score_threshold,
        inclusion=args.inclusion,
        trace_kinases=trace_kinases,
        svm_mode=args.svm_mode,
        sampling_trace_dir=args.sampling_trace_dir,
        skipped_trace_kinases=skipped_trace_kinases,
    )

    if skipped_trace_kinases:
        print(
            "Replay trace skipped kinases with stale or incompatible fixture rows: "
            + ", ".join(skipped_trace_kinases)
        )
    print(f"Done. Python prediction traces written to: {args.outdir}")


if __name__ == "__main__":
    main()
