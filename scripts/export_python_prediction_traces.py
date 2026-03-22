#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from phospy import KinasePredictor, PredictionSamplingTrace

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


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


def export_traces(
    result,
    combined_scores: pd.DataFrame,
    outdir: Path,
    *,
    top: int,
    score_threshold: float,
    inclusion: int,
    trace_kinases: list[str],
    sampling_trace_dir: str | None,
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

    initial_rows: list[dict[str, object]] = []
    probability_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []
    final_prediction_rows: list[dict[str, object]] = []
    final_top_rows: list[dict[str, object]] = []

    debug_traces = result.debug_traces or {}
    for kinase, trace in debug_traces.items():
        for ensemble_trace in trace.ensemble_traces:
            for draw, site in enumerate(ensemble_trace.initial_negative_sites, start=1):
                initial_rows.append(
                    {
                        "kinase": kinase,
                        "ensemble": ensemble_trace.ensemble_index,
                        "draw": draw,
                        "site": site,
                    }
                )

            for iteration_trace in ensemble_trace.iterations:
                labels = iteration_trace.labels
                probs = iteration_trace.probabilities
                for site in probs.index:
                    probability_rows.append(
                        {
                            "kinase": kinase,
                            "ensemble": ensemble_trace.ensemble_index,
                            "iteration": iteration_trace.iteration_index,
                            "site": site,
                            "label": str(int(labels.loc[site])),
                            "prob_class_1": float(probs.loc[site, "1"])
                            if "1" in probs.columns
                            else float("nan"),
                            "prob_class_2": float(probs.loc[site, "2"])
                            if "2" in probs.columns
                            else float("nan"),
                        }
                    )
                for draw, site in enumerate(
                    iteration_trace.sampled_positive_sites, start=1
                ):
                    sample_rows.append(
                        {
                            "kinase": kinase,
                            "ensemble": ensemble_trace.ensemble_index,
                            "iteration": iteration_trace.iteration_index,
                            "class_label": "1",
                            "draw": draw,
                            "site": site,
                        }
                    )
                for draw, site in enumerate(
                    iteration_trace.sampled_negative_sites, start=1
                ):
                    sample_rows.append(
                        {
                            "kinase": kinase,
                            "ensemble": ensemble_trace.ensemble_index,
                            "iteration": iteration_trace.iteration_index,
                            "class_label": "2",
                            "draw": draw,
                            "site": site,
                        }
                    )

            final_probs = ensemble_trace.final_prediction_probabilities
            for site in final_probs.index:
                final_prediction_rows.append(
                    {
                        "kinase": kinase,
                        "ensemble": ensemble_trace.ensemble_index,
                        "site": site,
                        "prob_class_1": float(final_probs.loc[site, "1"])
                        if "1" in final_probs.columns
                        else float("nan"),
                        "prob_class_2": float(final_probs.loc[site, "2"])
                        if "2" in final_probs.columns
                        else float("nan"),
                    }
                )
            for rank, site in enumerate(ensemble_trace.final_top_sites, start=1):
                final_top_rows.append(
                    {
                        "kinase": kinase,
                        "ensemble": ensemble_trace.ensemble_index,
                        "rank": rank,
                        "site": site,
                        "prob_class_1": float(final_probs.loc[site, "1"])
                        if "1" in final_probs.columns
                        else float("nan"),
                    }
                )

    pd.DataFrame(initial_rows).to_csv(
        outdir / "trace_initial_negatives.csv", index=False
    )
    pd.DataFrame(probability_rows).to_csv(
        outdir / "trace_iteration_probabilities.csv", index=False
    )
    pd.DataFrame(sample_rows).to_csv(
        outdir / "trace_iteration_samples.csv", index=False
    )
    pd.DataFrame(final_prediction_rows).to_csv(
        outdir / "trace_final_ensemble_predictions.csv", index=False
    )
    pd.DataFrame(final_top_rows).to_csv(
        outdir / "trace_final_ensemble_top.csv", index=False
    )

    readme = [
        "# Python prediction trace fixtures",
        "",
        "These files are generated from the native Python predictor with debug tracing enabled.",
        "They are intended for direct comparison with the R trace fixtures generated from the PhosR L6 example path.",
        "",
        f"Trace kinases: {', '.join(trace_kinases)}",
        "",
        f"Sampling trace override: {sampling_trace_dir or 'none'}",
        "",
        "Files:",
        "- trace_candidates.csv: ranked combined-score candidates for the traced kinases",
        "- trace_initial_negatives.csv: initial negative draw for each ensemble member",
        "- trace_iteration_probabilities.csv: per-iteration class probabilities on the base train set",
        "- trace_iteration_samples.csv: resampled site identities for each iteration and class",
        "- trace_final_ensemble_predictions.csv: final per-ensemble prediction probabilities for all sites",
        "- trace_final_ensemble_top.csv: final per-ensemble top-ranked sites",
    ]
    (outdir / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    combined_scores = pd.read_csv(args.combined_scores, index_col=0)
    trace_kinases = parse_csv_values(args.trace_kinases)

    predictor = KinasePredictor()
    sampling_trace = None
    if args.sampling_trace_dir:
        sampling_trace = PredictionSamplingTrace.from_trace_directory(
            args.sampling_trace_dir
        ).subset_kinases(trace_kinases)
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
        sampling_trace=sampling_trace,
    )

    export_traces(
        result,
        combined_scores,
        Path(args.outdir),
        top=args.top,
        score_threshold=args.score_threshold,
        inclusion=args.inclusion,
        trace_kinases=trace_kinases,
        sampling_trace_dir=args.sampling_trace_dir,
    )

    print(f"Done. Python prediction traces written to: {args.outdir}")


if __name__ == "__main__":
    main()
