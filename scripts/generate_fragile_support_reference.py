#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"

DEFAULT_SOURCE_DIR = ROOT / "tests" / "fixtures" / "r_reference_l6"
DEFAULT_OUTPUT_DIR = ROOT / "tests" / "fixtures" / "fragile_support_reference"
SELECTED_KINASES: tuple[str, ...] = (
    "MAPK1",
    "AKT1",
    "IRAK1",
    "PRKAA1",
    "PRKAA2",
    "LCK",
)
TARGET_TOP_THRESHOLD_ROWS: dict[str, int] = {
    "MAPK1": 24,
    "AKT1": 25,
    "IRAK1": 20,
    "PRKAA1": 20,
    "PRKAA2": 20,
    "LCK": 2,
}
SCORE_THRESHOLD = 0.8
INCLUSION = 20
TOP = 50


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate the curated fragile-support reference dataset derived from "
            "the committed L6 R reference outputs."
        )
    )
    parser.add_argument(
        "--source-dir",
        default=str(DEFAULT_SOURCE_DIR),
        help="Directory containing the source L6 reference outputs.",
    )
    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory where the fragile-support fixture family will be written.",
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
    grouped_map: Mapping[str, Sequence[str]],
    *,
    key_col: str,
    value_col: str,
) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for key, values in grouped_map.items():
        for value in values:
            rows.append({key_col: str(key), value_col: str(value)})
    return pd.DataFrame(rows, columns=[key_col, value_col])


def candidate_status(candidate_count: int) -> str:
    if candidate_count >= 30:
        return "robust"
    if candidate_count >= INCLUSION:
        return "just_above_inclusion"
    if candidate_count >= 3:
        return "just_below_inclusion"
    return "dropped"


def main() -> None:
    args = parse_args()
    source_dir = Path(args.source_dir)
    output_dir = Path(args.outdir)

    if str(SRC) not in sys.path:
        sys.path.insert(0, str(SRC))

    from phospy.profiles import build_kinase_substrate_profiles

    from phospy.prediction import (
        KinaseScorer,
        build_candidate_substrate_list,
        combine_profile_and_motif_scores,
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    phospho_matrix = read_indexed_csv(source_dir / "l6_phospho_matrix.csv")
    site_sequence_frame = pd.read_csv(source_dir / "l6_site_sequences.csv")
    site_sequence_frame["site_id"] = site_sequence_frame["site_id"].astype(str)
    site_sequence_frame["short_site_id"] = (
        site_sequence_frame["site_id"]
        .str.split(";")
        .map(lambda parts: f"{parts[1]};{parts[2]};" if len(parts) >= 3 else "")
    )
    site_sequences = site_sequence_frame.set_index("short_site_id")[
        "centralized_sequence"
    ]
    substrate_map_full = read_grouped_mapping(source_dir / "native_substrate_map.csv")
    motif_scores_full = read_indexed_csv(source_dir / "native_motif_scores.csv")
    motif_sizes_full = pd.read_csv(source_dir / "native_motif_sizes.csv").set_index(
        "kinase"
    )["motif_size"]
    combined_scores_full = read_indexed_csv(source_dir / "native_combined_scores.csv")

    selected_rows: set[str] = set()
    for kinase in SELECTED_KINASES:
        ranked = combined_scores_full.loc[:, kinase].sort_values(
            ascending=False,
            kind="mergesort",
        )
        selected_rows.update(
            ranked.loc[ranked > SCORE_THRESHOLD]
            .head(TARGET_TOP_THRESHOLD_ROWS[kinase])
            .index.tolist()
        )
        selected_rows.update(substrate_map_full[kinase])

    row_index = sorted(selected_rows)
    phospho_subset = phospho_matrix.loc[row_index].sort_index().sort_index(axis=1)
    sequence_subset = site_sequences.loc[row_index].sort_index()
    substrate_map_subset = {
        kinase: [site for site in substrate_map_full[kinase] if site in selected_rows]
        for kinase in SELECTED_KINASES
    }

    profile_result = build_kinase_substrate_profiles(
        substrate_map=substrate_map_subset,
        phospho_matrix=phospho_subset,
        min_substrates=1,
    )
    profile_scores = (
        KinaseScorer(profile_result.profile_matrix)
        .score_phosphosite_profiles(phospho_subset)
        .sort_index()
        .sort_index(axis=1)
    )
    motif_scores = (
        motif_scores_full.loc[row_index, list(SELECTED_KINASES)]
        .sort_index()
        .sort_index(axis=1)
    )
    motif_sizes = motif_sizes_full.loc[list(SELECTED_KINASES)].sort_index()
    combined_scores, combined_weights = combine_profile_and_motif_scores(
        motif_scores=motif_scores,
        profile_scores=profile_scores,
        motif_sizes=motif_sizes,
        profile_sizes=profile_result.substrate_counts.astype(float),
    )
    combined_scores = combined_scores.sort_index().sort_index(axis=1)
    combined_weights = combined_weights.sort_index()

    candidate_substrates = build_candidate_substrate_list(
        combined_scores=combined_scores,
        top=TOP,
        score_threshold=SCORE_THRESHOLD,
        inclusion=INCLUSION,
    )

    motif_sequence_rows: list[dict[str, str]] = []
    for kinase, sites in substrate_map_subset.items():
        for site in sites:
            motif_sequence_rows.append(
                {
                    "kinase": kinase,
                    "site_id": site,
                    "sequence": str(sequence_subset.loc[site]),
                }
            )
    motif_sequence_frame = pd.DataFrame(
        motif_sequence_rows,
        columns=["kinase", "site_id", "sequence"],
    ).sort_values(["kinase", "site_id"], kind="mergesort")

    summary_rows: list[dict[str, object]] = []
    for kinase in SELECTED_KINASES:
        candidate_count = int((combined_scores.loc[:, kinase] > SCORE_THRESHOLD).sum())
        summary_rows.append(
            {
                "kinase": kinase,
                "substrate_count": int(profile_result.substrate_counts.loc[kinase]),
                "motif_size": int(motif_sizes.loc[kinase]),
                "candidate_count": candidate_count,
                "included_in_candidate_list": kinase in candidate_substrates,
                "candidate_status": candidate_status(candidate_count),
                "motif_profile_spearman": float(
                    profile_scores.loc[:, kinase].corr(
                        motif_scores.loc[:, kinase],
                        method="spearman",
                    )
                ),
                "motif_weight": float(combined_weights.loc[kinase, "motif_weight"]),
                "profile_weight": float(combined_weights.loc[kinase, "profile_weight"]),
            }
        )
    summary = pd.DataFrame(summary_rows).sort_values("candidate_count")

    readme_lines = [
        "# Fragile-support reference dataset",
        "",
        "This dataset is a curated L6-derived reference family designed to stress decision fragility rather than broad coverage.",
        "",
        "It is intentionally smaller and more uneven than the full L6 reference family:",
        "- mixed kinase support counts",
        "- uneven motif/profile agreement",
        "- smaller and more fragile candidate pools",
        "- at least one dropped kinase, one kinase below inclusion, one just above inclusion, and multiple robust kinases",
        "",
        f"Selected kinases: {', '.join(SELECTED_KINASES)}",
        f"Row count: {len(row_index)}",
        f"Candidate selection settings: top={TOP}, score_threshold={SCORE_THRESHOLD}, inclusion={INCLUSION}",
        "",
        "This directory is not a blanket parity claim. It is a committed curated dataset for widening evidence beyond the main L6 path and for future seam expansion.",
        "",
        "Files:",
        "- phospho_matrix.csv: curated phosphosite matrix",
        "- site_sequences.csv: centralized sequences for the curated site index",
        "- substrate_map.csv: selected kinase-to-site mapping",
        "- motif_sequences.csv: flattened kinase motif sequences for the selected substrate sites",
        "- profile_matrix.csv / profile_sizes.csv: Python-built substrate profile reference outputs",
        "- profile_scores.csv / motif_scores.csv: deterministic scoring seam outputs",
        "- motif_sizes.csv / combined_scores.csv / combined_weights.csv: score-combination seam outputs",
        "- candidate_substrates.csv: candidate selection output under the configured threshold settings",
        "- screening_summary.csv: selection-summary table for the curated dataset",
    ]
    (output_dir / "README.md").write_text(
        "\n".join(readme_lines) + "\n", encoding="utf-8"
    )

    phospho_subset.to_csv(output_dir / "phospho_matrix.csv")
    (
        sequence_subset.rename("centralized_sequence")
        .rename_axis("site_id")
        .reset_index()
        .to_csv(output_dir / "site_sequences.csv", index=False)
    )
    flatten_mapping(
        substrate_map_subset,
        key_col="kinase",
        value_col="site_id",
    ).to_csv(output_dir / "substrate_map.csv", index=False)
    motif_sequence_frame.to_csv(output_dir / "motif_sequences.csv", index=False)
    profile_result.profile_matrix.sort_index().sort_index(axis=1).to_csv(
        output_dir / "profile_matrix.csv"
    )
    (
        profile_result.substrate_counts.rename("substrate_count")
        .rename_axis("kinase")
        .reset_index()
        .to_csv(output_dir / "profile_sizes.csv", index=False)
    )
    profile_scores.to_csv(output_dir / "profile_scores.csv")
    motif_scores.to_csv(output_dir / "motif_scores.csv")
    (
        motif_sizes.rename("motif_size")
        .rename_axis("kinase")
        .reset_index()
        .to_csv(output_dir / "motif_sizes.csv", index=False)
    )
    combined_scores.to_csv(output_dir / "combined_scores.csv")
    combined_weights.reset_index().to_csv(
        output_dir / "combined_weights.csv", index=False
    )
    flatten_mapping(
        candidate_substrates,
        key_col="kinase",
        value_col="site_id",
    ).to_csv(output_dir / "candidate_substrates.csv", index=False)
    summary.to_csv(output_dir / "screening_summary.csv", index=False)

    print(f"Wrote fragile-support reference dataset to: {output_dir}")


if __name__ == "__main__":
    main()
