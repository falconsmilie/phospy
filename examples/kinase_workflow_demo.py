#!/usr/bin/env python3
"""Run the supported first-run dataset-builder + kinase workflow lane."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    KinaseWorkflow,
)
from phospy.api import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    KinaseWorkflowResult,
    Organism,
    ReferencePreset,
)


def build_demo_dataset() -> AnalysisReadyPhosphoDataset:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.00, 0.70],
            "sample_b": [1.10, 0.80],
            "sample_c": [0.95, 0.75],
        },
        index=["TSC2;S939;", "GSK3B;S9;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["TSC2", "GSK3B"],
            "site": ["S939", "S9"],
            "site_sequence": [
                "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
                "_______MSGRPRTTSFAESCKPVQQPSAFG",
            ],
            "protein_id": ["TSC2", "GSK3B"],
        },
        index=phospho.index.copy(),
    )
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )


def run_demo(outdir: Path) -> tuple[KinaseWorkflowResult, dict[str, Path]]:
    outdir.mkdir(parents=True, exist_ok=True)
    dataset = build_demo_dataset()
    result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            activity_config=None,
        )
    )

    pred_mat_path = outdir / "pred_mat.csv"
    written: dict[str, Path] = {"pred_mat": pred_mat_path}
    result.prediction_result.pred_mat.to_csv(pred_mat_path)
    if result.activity_result is not None:
        weighted_activity_path = outdir / "weighted_activity.csv"
        ksea_scores_path = outdir / "ksea_scores.csv"
        result.activity_result.weighted_activity.to_csv(weighted_activity_path)
        result.activity_result.ksea_scores.to_csv(ksea_scores_path)
        written["weighted_activity"] = weighted_activity_path
        written["ksea_scores"] = ksea_scores_path

    return result, written


def main() -> None:
    with TemporaryDirectory(prefix="phospy-kinase-workflow-") as tmp_dir:
        result, written = run_demo(Path(tmp_dir))
        print("Kinase workflow demo")
        print("Resolved reference organism:", result.references.organism.value)
        print()
        print("Profile score shape:", result.scoring_result.profile_scores.shape)
        if result.scoring_result.motif_scores is not None:
            print("Motif score shape:", result.scoring_result.motif_scores.shape)
        else:
            print("Motif score table: unavailable")
        if result.scoring_result.combined_scores is not None:
            print("Combined score shape:", result.scoring_result.combined_scores.shape)
        else:
            print("Combined score table: unavailable")
        if result.scoring_result.weights is not None:
            print("Weight table shape:", result.scoring_result.weights.shape)
        else:
            print("Weight table: unavailable")
        print("Prediction matrix")
        print(result.prediction_result.pred_mat.round(4))
        if result.prediction_result.substrate_list is not None:
            print()
            print("Prediction substrate list")
            print(result.prediction_result.substrate_list.round(4))
        if result.activity_result is not None:
            print()
            print("Weighted activity")
            print(result.activity_result.weighted_activity.round(4))
            print()
            print("KSEA scores")
            print(result.activity_result.ksea_scores.round(4))
        print()
        print("Written files")
        print("\n".join(str(path) for path in written.values()))


if __name__ == "__main__":
    main()
