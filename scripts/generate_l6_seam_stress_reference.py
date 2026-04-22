#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
SOURCE_DIR = ROOT / "tests" / "fixtures" / "rewrite_parity" / "r_reference_l6"
DEFAULT_OUTPUT_DIR = (
    ROOT / "tests" / "fixtures" / "rewrite_parity" / "r_reference_l6_seam_stress"
)
TRACE_SOURCE_DIR = SOURCE_DIR / "prediction_trace"

CANDIDATE_TOP = 50
CANDIDATE_SCORE_THRESHOLD = 0.8
CANDIDATE_INCLUSION = 20
TRACE_TOP = 30
TRACE_SCORE_THRESHOLD = 0.6
TRACE_INCLUSION = 5
TRACE_ENSEMBLE_SIZE = 10
TRACE_N_ITERATIONS = 5
TRACE_RANDOM_STATE = 1
TRACE_DEBUG_TOP_N = 10
PREFERRED_TRACE_KINASES: tuple[str, ...] = ("MAPK1", "IRAK1")

TRACE_TABLES = (
    "trace_candidates.csv",
    "trace_selected_candidates.csv",
    "trace_negative_pool.csv",
    "trace_initial_negatives.csv",
    "trace_iteration_labels.csv",
    "trace_iteration_probabilities.csv",
    "trace_iteration_probability_parameters.csv",
    "trace_iteration_decision_values.csv",
    "trace_iteration_resampling_weights.csv",
    "trace_iteration_samples.csv",
    "trace_final_ensemble_predictions.csv",
    "trace_final_ensemble_decision_values.csv",
    "trace_final_ensemble_top.csv",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a smaller R-backed seam-stress fixture by slicing the committed "
            "L6 R reference outputs and trace rows."
        )
    )
    parser.add_argument(
        "--source-dir",
        default=str(SOURCE_DIR),
        help="Directory containing the committed L6 R reference outputs.",
    )
    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory where the seam-stress fixture will be written.",
    )
    return parser.parse_args()


def read_indexed_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0)


def read_grouped_mapping(path: Path) -> dict[str, list[str]]:
    frame = pd.read_csv(path)
    return {
        str(kinase): group.loc[:, "site_id"].astype(str).tolist()
        for kinase, group in frame.groupby("kinase", sort=False)
    }


def flatten_mapping(
    grouped_map: dict[str, list[str]], *, key_col: str, value_col: str
) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for key, values in grouped_map.items():
        for value in values:
            rows.append({key_col: str(key), value_col: str(value)})
    return pd.DataFrame(rows, columns=[key_col, value_col])


def build_prediction_top_table(pred_mat: pd.DataFrame, top_n: int = 30) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for kinase in pred_mat.columns:
        ordered = pred_mat.loc[:, kinase].sort_values(ascending=False, kind="mergesort")
        for rank, (site, score) in enumerate(ordered.head(top_n).items(), start=1):
            rows.append(
                {
                    "kinase": kinase,
                    "site_id": site,
                    "pred_score": float(score),
                    "rank": rank,
                }
            )
    return pd.DataFrame(rows)


def candidate_status(candidate_count: int) -> str:
    if candidate_count >= 30:
        return "robust"
    if candidate_count >= CANDIDATE_INCLUSION:
        return "just_above_inclusion"
    if candidate_count >= 3:
        return "just_below_inclusion"
    return "dropped"


def select_row_index(trace_dir: Path, trace_kinases: tuple[str, ...]) -> list[str]:
    trace_candidates = pd.read_csv(trace_dir / "trace_candidates.csv")
    trace_initial = pd.read_csv(trace_dir / "trace_initial_negatives.csv")
    trace_samples = pd.read_csv(trace_dir / "trace_iteration_samples.csv")
    trace_final_top = pd.read_csv(trace_dir / "trace_final_ensemble_top.csv")

    selected_rows: set[str] = set()
    for kinase in trace_kinases:
        selected_rows.update(
            trace_candidates.loc[
                (trace_candidates.loc[:, "kinase"].astype(str) == kinase)
                & (trace_candidates.loc[:, "selected_candidate"]),
                "site",
            ].astype(str)
        )
        for table in (trace_initial, trace_samples, trace_final_top):
            selected_rows.update(
                table.loc[
                    table.loc[:, "kinase"].astype(str) == kinase,
                    "site",
                ].astype(str)
            )
    return sorted(selected_rows)


def filter_trace_table(table: pd.DataFrame, row_index: set[str]) -> pd.DataFrame:
    filtered = table.copy()
    if "site" in filtered.columns:
        filtered = filtered.loc[filtered.loc[:, "site"].astype(str).isin(row_index)]
    return filtered.reset_index(drop=True)


def resolve_trace_kinases(trace_dir: Path) -> tuple[str, ...]:
    trace_candidates = pd.read_csv(trace_dir / "trace_candidates.csv")
    available_kinases = (
        trace_candidates.loc[:, "kinase"].astype(str).drop_duplicates().tolist()
    )
    if not available_kinases:
        msg = "No traced kinases were found in the source trace_candidates.csv table."
        raise ValueError(msg)

    preferred = [
        kinase for kinase in PREFERRED_TRACE_KINASES if kinase in available_kinases
    ]
    ordered_kinases = list(dict.fromkeys(preferred + available_kinases))
    return tuple(ordered_kinases)


def available_trace_tables(trace_dir: Path) -> tuple[str, ...]:
    return tuple(name for name in TRACE_TABLES if (trace_dir / name).exists())


def main() -> None:
    args = parse_args()
    source_dir = Path(args.source_dir)
    outdir = Path(args.outdir)
    trace_source_dir = source_dir / "prediction_trace"

    if str(SRC) not in sys.path:
        sys.path.insert(0, str(SRC))

    from phospy.prediction import (
        build_candidate_substrate_list,
        combine_profile_and_motif_scores,
    )

    outdir.mkdir(parents=True, exist_ok=True)
    trace_outdir = outdir / "prediction_trace"
    trace_outdir.mkdir(parents=True, exist_ok=True)

    trace_kinases = resolve_trace_kinases(trace_source_dir)
    trace_tables = available_trace_tables(trace_source_dir)
    if not trace_tables:
        msg = "No source prediction trace tables were found for seam-stress generation."
        raise ValueError(msg)

    for existing_trace_file in trace_outdir.glob("trace_*.csv"):
        existing_trace_file.unlink()

    row_index = select_row_index(trace_source_dir, trace_kinases)
    row_index_set = set(row_index)

    profile_scores_full = read_indexed_csv(source_dir / "native_profile_scores.csv")
    motif_scores_full = read_indexed_csv(source_dir / "native_motif_scores.csv")
    combined_scores_full = read_indexed_csv(source_dir / "native_combined_scores.csv")
    pred_mat_full = read_indexed_csv(source_dir / "predMat.csv")
    combined_weights = (
        pd.read_csv(source_dir / "native_combined_weights.csv")
        .set_index("kinase")
        .sort_index()
    )
    motif_sizes = (
        pd.read_csv(source_dir / "native_motif_sizes.csv")
        .set_index("kinase")
        .loc[:, "motif_size"]
        .sort_index()
    )
    substrate_map = read_grouped_mapping(source_dir / "native_substrate_map.csv")
    profile_sizes = pd.Series(
        {kinase: len(sites) for kinase, sites in substrate_map.items()},
        name="substrate_count",
    ).sort_index()

    profile_scores = profile_scores_full.loc[row_index].sort_index().sort_index(axis=1)
    motif_scores = motif_scores_full.loc[row_index].sort_index().sort_index(axis=1)
    combined_scores = (
        combined_scores_full.loc[row_index].sort_index().sort_index(axis=1)
    )
    pred_mat = pred_mat_full.loc[row_index].sort_index().sort_index(axis=1)

    actual_combined, actual_weights = combine_profile_and_motif_scores(
        motif_scores=motif_scores,
        profile_scores=profile_scores,
        motif_sizes=motif_sizes,
        profile_sizes=profile_sizes.astype(float),
    )
    actual_combined = actual_combined.sort_index().sort_index(axis=1)
    actual_weights = actual_weights.sort_index()
    combined_max_abs = float((actual_combined - combined_scores).abs().max().max())
    if combined_max_abs > 1e-12:
        raise ValueError(
            "Selected seam-stress slice no longer reproduces the committed L6 combined scores. "
            f"Max abs diff: {combined_max_abs}"
        )
    weight_max_abs = float((actual_weights - combined_weights).abs().max().max())
    if weight_max_abs > 1e-12:
        raise ValueError(
            "Selected seam-stress slice no longer reproduces the committed L6 combined weights. "
            f"Max abs diff: {weight_max_abs}"
        )

    candidate_substrates = build_candidate_substrate_list(
        combined_scores=combined_scores,
        top=CANDIDATE_TOP,
        score_threshold=CANDIDATE_SCORE_THRESHOLD,
        inclusion=CANDIDATE_INCLUSION,
    )
    prediction_top30 = build_prediction_top_table(pred_mat, top_n=30)

    screening_summary = pd.DataFrame(
        {
            "kinase": combined_scores.columns,
            "motif_size": motif_sizes.reindex(combined_scores.columns)
            .astype(int)
            .to_list(),
            "substrate_count": profile_sizes.reindex(combined_scores.columns)
            .astype(int)
            .to_list(),
            "candidate_count": [
                int(
                    (
                        combined_scores.loc[:, kinase]
                        .sort_values(ascending=False, kind="mergesort")
                        .head(CANDIDATE_TOP)
                        > CANDIDATE_SCORE_THRESHOLD
                    ).sum()
                )
                for kinase in combined_scores.columns
            ],
        }
    )
    screening_summary["included_in_candidate_list"] = screening_summary.loc[
        :, "kinase"
    ].isin(candidate_substrates)
    screening_summary["candidate_status"] = screening_summary.loc[
        :, "candidate_count"
    ].map(candidate_status)
    screening_summary = screening_summary.sort_values(
        ["candidate_count", "kinase"], kind="mergesort"
    )

    candidate_counts_for_trace: dict[str, int] = {}
    trace_candidates_full = pd.read_csv(trace_source_dir / "trace_candidates.csv")
    for kinase in trace_kinases:
        trace_count = int(
            trace_candidates_full.loc[
                (trace_candidates_full.loc[:, "kinase"].astype(str) == kinase)
                & (trace_candidates_full.loc[:, "selected_candidate"]),
                :,
            ].shape[0]
        )
        current_count = len(
            build_candidate_substrate_list(
                combined_scores.loc[:, [kinase]],
                top=TRACE_TOP,
                score_threshold=TRACE_SCORE_THRESHOLD,
                inclusion=TRACE_INCLUSION,
            ).get(kinase, [])
        )
        candidate_counts_for_trace[kinase] = current_count
        if current_count != trace_count:
            raise ValueError(
                f"Trace candidate alignment failed for {kinase}: current_count={current_count}, trace_count={trace_count}"
            )

    readme_lines = [
        "# L6 seam-stress reference dataset",
        "",
        "This directory is a smaller R-backed seam fixture derived by slicing the committed L6 reference outputs.",
        "",
        "It is intentionally not a second independent biological source dataset. Instead, it keeps the reference provenance R-backed while narrowing the row set to stress different native seam decisions:",
        "- thinner candidate pools for selected kinases under stricter candidate-selection settings",
        "- retained exact R sampling replay for a smaller traced kinase subset",
        "- full overlap-kinase weight behaviour preserved by keeping the committed L6 score columns intact",
        "",
        f"Row count: {len(row_index)}",
        f"Kinase count: {combined_scores.shape[1]}",
        f"Trace kinases: {', '.join(trace_kinases)}",
        f"Candidate-selection settings: top={CANDIDATE_TOP}, score_threshold={CANDIDATE_SCORE_THRESHOLD}, inclusion={CANDIDATE_INCLUSION}",
        f"Trace replay settings: top={TRACE_TOP}, score_threshold={TRACE_SCORE_THRESHOLD}, inclusion={TRACE_INCLUSION}, ensemble_size={TRACE_ENSEMBLE_SIZE}, n_iterations={TRACE_N_ITERATIONS}, random_state={TRACE_RANDOM_STATE}, debug_top_n={TRACE_DEBUG_TOP_N}",
        "",
        "Provenance:",
        "- profile_scores.csv / motif_scores.csv / combined_scores.csv are direct row slices of the committed L6 R reference outputs",
        "- motif_sizes.csv / profile_sizes.csv / combined_weights.csv are direct L6 R-backed seam metadata tables",
        "- predMat.csv and prediction_top30.csv are R-backed prediction outputs sliced to the seam-stress row set",
        "- prediction_trace/* is filtered from the committed L6 R trace for the traced kinases only",
        "",
        "Files:",
        "- profile_scores.csv / motif_scores.csv / combined_scores.csv: seam-score reference tables",
        "- motif_sizes.csv / profile_sizes.csv / combined_weights.csv: score-combination weight inputs and outputs",
        "- candidate_substrates.csv / screening_summary.csv: stricter candidate-selection seam references",
        "- predMat.csv / prediction_top30.csv: R-backed prediction ranking references on the seam-stress row set",
        "- prediction_trace/*: filtered R sampling and debug-trace reference tables for replay checks",
    ]
    (outdir / "README.md").write_text("\n".join(readme_lines) + "\n", encoding="utf-8")

    profile_scores.to_csv(outdir / "profile_scores.csv")
    motif_scores.to_csv(outdir / "motif_scores.csv")
    combined_scores.to_csv(outdir / "combined_scores.csv")
    pred_mat.to_csv(outdir / "predMat.csv")
    prediction_top30.to_csv(outdir / "prediction_top30.csv", index=False)
    screening_summary.to_csv(outdir / "screening_summary.csv", index=False)
    (
        motif_sizes.rename("motif_size")
        .rename_axis("kinase")
        .reset_index()
        .to_csv(outdir / "motif_sizes.csv", index=False)
    )
    (
        profile_sizes.rename("substrate_count")
        .rename_axis("kinase")
        .reset_index()
        .to_csv(outdir / "profile_sizes.csv", index=False)
    )
    combined_weights.reset_index().to_csv(outdir / "combined_weights.csv", index=False)
    flatten_mapping(
        candidate_substrates,
        key_col="kinase",
        value_col="site_id",
    ).to_csv(outdir / "candidate_substrates.csv", index=False)

    for table_name in trace_tables:
        filtered = filter_trace_table(
            pd.read_csv(trace_source_dir / table_name).loc[
                lambda frame: frame.loc[:, "kinase"].astype(str).isin(trace_kinases)
            ],
            row_index_set,
        )
        filtered.to_csv(trace_outdir / table_name, index=False)

    print(f"Done. L6 seam-stress reference fixtures written to: {outdir}")


if __name__ == "__main__":
    main()
