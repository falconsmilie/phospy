#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

FIXTURE_DIR = Path("tests/fixtures/synthetic_adaptive_sampling_edge")
TOP = 4
SCORE_THRESHOLD = 0.85
INCLUSION = 2
ENSEMBLE_SIZE = 1
N_ITERATIONS = 2
RANDOM_STATE = 1
DEBUG_TOP_N = 4
TRACE_KINASES = ["KINASE_A", "KINASE_B"]


def build_combined_scores() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "KINASE_A": [0.95, 0.95, 0.88, 0.20, 0.20, 0.10],
            "KINASE_B": [0.10, 0.20, 0.20, 0.88, 0.95, 0.95],
        },
        index=["SITE_A", "SITE_B", "SITE_C", "SITE_D", "SITE_E", "SITE_F"],
    )


def build_trace_candidates(combined_scores: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for kinase in TRACE_KINASES:
        ordered = combined_scores.loc[:, kinase].sort_values(
            ascending=False,
            kind="mergesort",
        )
        top_sites = ordered.head(TOP)
        selected_sites = top_sites.loc[top_sites > SCORE_THRESHOLD].index.tolist()
        for rank, (site, score) in enumerate(ordered.items(), start=1):
            rows.append(
                {
                    "kinase": kinase,
                    "rank": rank,
                    "site": site,
                    "combined_score": float(score),
                    "within_top": rank <= TOP,
                    "above_threshold": float(score) > SCORE_THRESHOLD,
                    "selected_candidate": site in selected_sites,
                }
            )
    return pd.DataFrame(rows)


def write_sampling_override(trace_dir: Path) -> None:
    trace_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 1, "site": "SITE_D"},
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 2, "site": "SITE_E"},
            {"kinase": "KINASE_A", "ensemble": 1, "draw": 3, "site": "SITE_F"},
            {"kinase": "KINASE_B", "ensemble": 1, "draw": 1, "site": "SITE_C"},
            {"kinase": "KINASE_B", "ensemble": 1, "draw": 2, "site": "SITE_B"},
            {"kinase": "KINASE_B", "ensemble": 1, "draw": 3, "site": "SITE_A"},
        ]
    ).to_csv(trace_dir / "trace_initial_negatives.csv", index=False)
    pd.DataFrame(
        [
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 1,
                "draw": 1,
                "site": "SITE_A",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 1,
                "draw": 2,
                "site": "SITE_B",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 1,
                "draw": 3,
                "site": "SITE_C",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 2,
                "draw": 1,
                "site": "SITE_D",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 2,
                "draw": 2,
                "site": "SITE_E",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 2,
                "draw": 3,
                "site": "SITE_F",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 1,
                "draw": 1,
                "site": "SITE_B",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 1,
                "draw": 2,
                "site": "SITE_B",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 1,
                "draw": 3,
                "site": "SITE_A",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 2,
                "draw": 1,
                "site": "SITE_E",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 2,
                "draw": 2,
                "site": "SITE_D",
            },
            {
                "kinase": "KINASE_A",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 2,
                "draw": 3,
                "site": "SITE_F",
            },
            {
                "kinase": "KINASE_B",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 1,
                "draw": 1,
                "site": "SITE_E",
            },
            {
                "kinase": "KINASE_B",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 1,
                "draw": 2,
                "site": "SITE_F",
            },
            {
                "kinase": "KINASE_B",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 1,
                "draw": 3,
                "site": "SITE_D",
            },
            {
                "kinase": "KINASE_B",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 2,
                "draw": 1,
                "site": "SITE_C",
            },
            {
                "kinase": "KINASE_B",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 2,
                "draw": 2,
                "site": "SITE_B",
            },
            {
                "kinase": "KINASE_B",
                "ensemble": 1,
                "iteration": 1,
                "class_label": 2,
                "draw": 3,
                "site": "SITE_A",
            },
            {
                "kinase": "KINASE_B",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 1,
                "draw": 1,
                "site": "SITE_F",
            },
            {
                "kinase": "KINASE_B",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 1,
                "draw": 2,
                "site": "SITE_E",
            },
            {
                "kinase": "KINASE_B",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 1,
                "draw": 3,
                "site": "SITE_E",
            },
            {
                "kinase": "KINASE_B",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 2,
                "draw": 1,
                "site": "SITE_B",
            },
            {
                "kinase": "KINASE_B",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 2,
                "draw": 2,
                "site": "SITE_C",
            },
            {
                "kinase": "KINASE_B",
                "ensemble": 1,
                "iteration": 2,
                "class_label": 2,
                "draw": 3,
                "site": "SITE_A",
            },
        ]
    ).to_csv(trace_dir / "trace_iteration_samples.csv", index=False)


def write_readme(outdir: Path) -> None:
    readme = [
        "# Synthetic adaptive-sampling decision seam fixtures",
        "",
        "These fixtures are intentionally small and deterministic.",
        "They complement the R-backed L6 parity traces by pinning edge-case replay behaviour:",
        "",
        "- tied candidate scores with stable mergesort ordering",
        "- tiny candidate and negative-pool sizes",
        "- explicit per-iteration sampling overrides",
        "- exact final top-site decisions on a deterministic replay path",
        "",
        "They are not a standalone claim of PhosR parity.",
    ]
    (outdir / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")


def main() -> None:
    from phospy.prediction import (
        KinasePredictor,
        PredictionSamplingTrace,
        prediction_debug_trace_tables,
    )

    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    combined_scores = build_combined_scores()
    combined_scores.to_csv(FIXTURE_DIR / "combined_scores.csv")

    trace_dir = FIXTURE_DIR
    write_sampling_override(trace_dir)
    build_trace_candidates(combined_scores).to_csv(
        FIXTURE_DIR / "trace_candidates.csv",
        index=False,
    )

    sampling_trace = PredictionSamplingTrace.from_trace_directory(trace_dir)
    result = KinasePredictor(svm_mode="r_parity").predict(
        combined_scores=combined_scores,
        ensemble_size=ENSEMBLE_SIZE,
        top=TOP,
        score_threshold=SCORE_THRESHOLD,
        inclusion=INCLUSION,
        n_iterations=N_ITERATIONS,
        random_state=RANDOM_STATE,
        capture_debug_trace=True,
        debug_kinases=TRACE_KINASES,
        debug_top_n=DEBUG_TOP_N,
        sampling_trace=sampling_trace,
    )
    for name, table in prediction_debug_trace_tables(result).items():
        table.to_csv(FIXTURE_DIR / f"{name}.csv", index=False)

    write_readme(FIXTURE_DIR)


if __name__ == "__main__":
    main()
