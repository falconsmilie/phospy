#!/usr/bin/env python3
"""Run the preferred 1.5.0 dataset-builder to kinase workflow lane."""

from __future__ import annotations

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
            "localisation_confidence": [0.95] * phospho.shape[0],
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


def run_demo() -> KinaseWorkflowResult:
    dataset = build_demo_dataset()
    return KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            activity_config=None,
            site_sequence_conflict_policy="prefer_reference",
        )
    )


def main() -> None:
    result = run_demo()
    print("Preferred 1.5.0 kinase workflow lane")
    print("Dataset organism:", result.dataset.organism.value)
    print("Reference input: ReferencePreset.AUTO")
    print("Resolved reference organism:", result.references.organism.value)
    print("Profile score shape:", result.scoring_result.profile_scores.shape)
    print("Prediction matrix")
    print(result.prediction_result.pred_mat.round(4))


if __name__ == "__main__":
    main()
