from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from legacy_archive.phospy_legacy.activities.scoring import (
    build_kinase_target_table as legacy_build_kinase_target_table,
)
from legacy_archive.phospy_legacy.activities.scoring import (
    compute_activity_from_inputs as legacy_compute_activity_from_inputs,
)
from legacy_archive.phospy_legacy.activities.scoring import (
    count_predicted_targets as legacy_count_predicted_targets,
)
from legacy_archive.phospy_legacy.validation.requests.analysis import (
    validate_analysis_request as legacy_validate_analysis_request,
)

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FIXTURE_DIR = ROOT / "tests" / "fixtures" / "rewrite_parity" / "r_reference_l6"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Archival utility: regenerate activity donor snapshot tables for "
            "rewrite parity fixtures."
        )
    )
    parser.add_argument(
        "--fixture-dir",
        type=Path,
        default=DEFAULT_FIXTURE_DIR,
        help="Directory containing predMat.csv and l6_phospho_matrix.csv.",
    )
    parser.add_argument("--threshold", type=float, default=0.6)
    parser.add_argument("--min-substrates", type=int, default=3)
    parser.add_argument("--top-n-substrates", type=int, default=20)
    return parser.parse_args()


def _materialize_activity_outputs(
    *,
    fixture_dir: Path,
    threshold: float,
    min_substrates: int,
    top_n_substrates: int,
) -> dict[str, pd.DataFrame | pd.Series]:
    pred_mat = pd.read_csv(fixture_dir / "predMat.csv", index_col=0)
    phospho_matrix = pd.read_csv(fixture_dir / "l6_phospho_matrix.csv", index_col=0)

    validated = legacy_validate_analysis_request(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        threshold=threshold,
        min_substrates=min_substrates,
        top_n_substrates=top_n_substrates,
    )
    weighted_activity, ksea_scores, ksea_counts = legacy_compute_activity_from_inputs(
        validated
    )
    weighted_activity.index.name = "kinase"
    ksea_scores.index.name = "kinase"

    target_counts = legacy_count_predicted_targets(
        pred_mat=validated.pred_mat,
        threshold=validated.threshold,
    )
    target_counts.index.name = "kinase"
    target_counts.name = "n_targets"
    target_table = legacy_build_kinase_target_table(
        pred_mat=validated.pred_mat,
        threshold=validated.threshold,
    )

    return {
        "kinase_activity_matrix": weighted_activity,
        "ksea_scores": ksea_scores,
        "ksea_counts": ksea_counts.rename("n_substrates"),
        "kinase_target_counts": target_counts,
        "kinase_target_table": target_table,
    }


def _write_outputs(
    outputs: dict[str, pd.DataFrame | pd.Series],
    *,
    fixture_dir: Path,
) -> None:
    weighted_activity = outputs["kinase_activity_matrix"]
    assert isinstance(weighted_activity, pd.DataFrame)
    weighted_activity.to_csv(fixture_dir / "kinase_activity_matrix.csv")

    ksea_scores = outputs["ksea_scores"]
    assert isinstance(ksea_scores, pd.DataFrame)
    ksea_scores.to_csv(fixture_dir / "ksea_scores.csv")

    ksea_counts = outputs["ksea_counts"]
    assert isinstance(ksea_counts, pd.Series)
    ksea_counts.to_frame(name="n_substrates").to_csv(fixture_dir / "ksea_counts.csv")

    target_counts = outputs["kinase_target_counts"]
    assert isinstance(target_counts, pd.Series)
    target_counts.to_frame(name="n_targets").to_csv(
        fixture_dir / "kinase_target_counts.csv"
    )

    target_table = outputs["kinase_target_table"]
    assert isinstance(target_table, pd.DataFrame)
    target_table.to_csv(fixture_dir / "kinase_target_table.csv", index=False)


def main() -> int:
    args = _parse_args()
    fixture_dir = args.fixture_dir.resolve()
    outputs = _materialize_activity_outputs(
        fixture_dir=fixture_dir,
        threshold=args.threshold,
        min_substrates=args.min_substrates,
        top_n_substrates=args.top_n_substrates,
    )
    _write_outputs(outputs, fixture_dir=fixture_dir)
    print(f"wrote activity donor snapshot outputs under: {fixture_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
