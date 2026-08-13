#!/usr/bin/env python3
"""Run the supported kinase workflow from a dataset-builder input."""

from __future__ import annotations

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    KinaseWorkflow,
)
from phospy.advanced import (
    KinaseReliabilityProfile,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api import (
    DatasetBuildRequest,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
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
        index=["TSC2;S939;", "GSK3A;S21;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["TSC2", "GSK3A"],
            "site": ["S939", "S21"],
            "site_sequence": [
                "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
                "PSGGGPGGSGRARTSSFAEPGGGGGGGGGGP",
            ],
            "display_id": ["TSC2;S939;", "GSK3A;S21;"],
            "organism": ["rat", "rat"],
            "protein_namespace": ["protein_id", "protein_id"],
            "protein_identifier": ["TSC2", "GSK3A"],
            "localisation_confidence": [0.95] * phospho.shape[0],
            # Required later by signalome; kinase identity itself uses site_key.
            "protein_group_id": ["TSC2", "GSK3A"],
        },
        index=phospho.index.copy(),
    )
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                # Fail fast on missing/below-threshold localisation confidence
                # to avoid ambiguous site-level kinase interpretation.
                localisation=DatasetLocalisationConfig(
                    mode="require_threshold",
                    confidence_column="localisation_confidence",
                    min_confidence=0.75,
                )
            ),
        )
    )


def run_demo() -> KinaseWorkflowResult:
    dataset = build_demo_dataset()
    return KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(
                reliability_profile=KinaseReliabilityProfile.CUSTOM,
                reference_context_compatibility_policy=(
                    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                ),
            ),
            activity_config=None,
            site_sequence_conflict_policy="prefer_reference",
        )
    )


def main() -> None:
    result = run_demo()
    print("Supported kinase workflow")
    print("Dataset organism:", result.dataset.organism.value)
    print(
        result.dataset.site_metadata.loc[
            :,
            [
                "site_key",
                "display_id",
                "gene_symbol",
                "site",
                "protein_namespace",
                "protein_identifier",
                "protein_group_id",
            ],
        ]
    )
    print("Reference input: ReferencePreset.AUTO")
    print("Resolved reference organism:", result.references.organism.value)
    print("Profile score shape:", result.scoring_result.profile_scores.shape)
    print("Prediction matrix")
    print(result.prediction_result.pred_mat.round(4))


if __name__ == "__main__":
    main()
