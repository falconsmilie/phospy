#!/usr/bin/env python3
"""Run the supported builder + simple kinase workflow route."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DatasetBuildRequest,
    KinaseWorkflow,
    Organism,
    ReferencePreset,
    SimpleKinaseWorkflowRequest,
    SimpleKinaseWorkflowResult,
)


def build_demo_dataset() -> AnalysisReadyPhosphoDataset:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.00, 0.70],
            "sample_b": [1.15, 0.80],
            "sample_c": [0.95, 0.75],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B"],
            "site": ["Y182", "S9"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
            ],
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


def run_demo(outdir: Path) -> tuple[SimpleKinaseWorkflowResult, dict[str, Path]]:
    outdir.mkdir(parents=True, exist_ok=True)
    dataset = build_demo_dataset()
    result = KinaseWorkflow().run(
        SimpleKinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
        )
    )

    if result.activity_result is None:
        raise RuntimeError("activity_result is expected in this demo")

    pred_mat_path = outdir / "pred_mat.csv"
    activity_path = outdir / "activity_scores.csv"
    result.prediction_result.pred_mat.to_csv(pred_mat_path)
    result.activity_result.activity_scores.to_csv(activity_path)

    return result, {"pred_mat": pred_mat_path, "activity_scores": activity_path}


def main() -> None:
    with TemporaryDirectory(prefix="phospy-simple-workflow-") as tmp_dir:
        result, written = run_demo(Path(tmp_dir))
        print("Simple kinase workflow demo")
        print("Resolved reference organism:", result.references.organism.value)
        print()
        print("Profile score shape:", result.scoring_result.profile_scores.shape)
        print("Prediction matrix")
        print(result.prediction_result.pred_mat.round(4))
        print()
        print("Activity scores")
        print(result.activity_result.activity_scores.round(4))
        print()
        print("Written files")
        print("\n".join(str(path) for path in written.values()))


if __name__ == "__main__":
    main()
