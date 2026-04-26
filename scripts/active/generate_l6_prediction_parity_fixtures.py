#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
DEFAULT_OUTPUT_DIR = (
    ROOT / "tests" / "fixtures" / "rewrite_parity" / "r_reference_l6_prediction"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate rewrite-owned L6 scoring/prediction parity fixtures for "
            "the supported strict motif-sequence-validation lane."
        )
    )
    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory where native_* and predMat.csv fixtures are written.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if str(SRC) not in sys.path:
        sys.path.insert(0, str(SRC))
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))

    from phospy.api import (
        KinasePredictionConfig,
        KinaseScoringConfig,
        KinaseWorkflow,
        KinaseWorkflowRequest,
        ReferencePreset,
    )
    from tests.support.rewrite_fixture_data import build_rat_l6_dataset

    dataset = build_rat_l6_dataset(n_sites=None)
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                include_diagnostic_scoring_tables=True,
            ),
            prediction_config=KinasePredictionConfig(
                top_k=30,
                deterministic_max_selected_kinases=10,
                adaptive_ensemble_runs=10,
                mode="adaptive_ensemble",
                adaptive_policy="stable",
                n_iterations=5,
                random_state=1,
            ),
            activity_config=None,
        )
    )

    profile_scores = result.scoring_result.profile_scores.copy(deep=True)
    rank_weighted_fusion_scores = result.scoring_result.rank_weighted_fusion_scores
    if rank_weighted_fusion_scores is None:
        raise RuntimeError(
            "rank_weighted_fusion_scores missing in L6 fixture generation"
        )
    score_fusion_weights = result.scoring_result.score_fusion_weights
    if score_fusion_weights is None:
        raise RuntimeError("score_fusion_weights missing in L6 fixture generation")
    pred_mat = result.prediction_result.pred_mat.copy(deep=True)
    substrate_list = result.prediction_result.substrate_list
    if substrate_list is None:
        raise RuntimeError("substrate_list missing in L6 fixture generation")

    profile_scores.to_csv(output_dir / "native_profile_scores.csv")
    rank_weighted_fusion_scores.to_csv(
        output_dir / "native_rank_weighted_fusion_scores.csv"
    )
    score_fusion_weights.to_csv(output_dir / "native_score_fusion_weights.csv")
    pred_mat.to_csv(output_dir / "predMat.csv")

    top30 = (
        substrate_list.rename(
            columns={
                "substrate_site": "site_id",
                "score": "pred_score",
            }
        )
        .loc[:, ["kinase", "site_id", "pred_score", "rank"]]
        .sort_values(["kinase", "rank"], kind="mergesort")
        .reset_index(drop=True)
    )
    top30.to_csv(output_dir / "native_prediction_top30.csv", index=False)

    candidates = top30.loc[:, ["kinase", "site_id"]].copy(deep=True)
    candidates.to_csv(output_dir / "native_candidate_substrates.csv", index=False)


if __name__ == "__main__":
    main()
