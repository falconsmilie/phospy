#!/usr/bin/env python3
"""Build an analysis-ready dataset in the preferred 1.5.0 first-run lane."""

from __future__ import annotations

import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset
from phospy.api import (
    DatasetBuildRequest,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    Organism,
)


def _example_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 0.7],
            "sample_b": [1.2, 0.8],
            "sample_c": [0.9, 0.75],
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
            # Signalome has a separate explicit protein_id requirement.
            "protein_id": ["TSC2", "GSK3A"],
        },
        index=phospho.index.copy(),
    )
    return phospho, site_metadata


def build_demo_dataset() -> AnalysisReadyPhosphoDataset:
    phospho, site_metadata = _example_tables()
    request = DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        input_intensity_scale="linear",
        preprocessing_config=DatasetPreprocessingConfig(
            # Fail fast on missing/low localisation confidence so downstream
            # site-level interpretation does not rely on ambiguous site mapping.
            localisation=DatasetLocalisationConfig(
                mode="require_threshold",
                confidence_column="localisation_confidence",
                min_confidence=0.75,
            )
        ),
    )
    return AnalysisReadyDatasetBuilder().run(request)


def main() -> None:
    dataset = build_demo_dataset()
    print("Preferred 1.5.0 dataset builder lane")
    print("Input route: pandas DataFrame values")
    print("Bundled reference-compatible organism:", dataset.organism.value)
    print("Phospho shape:", dataset.phospho.shape)
    print("Analysis-ready row index:", dataset.phospho.index.name)
    print(
        dataset.site_metadata.loc[
            :,
            [
                "site_key",
                "display_id",
                "gene_symbol",
                "site",
                "organism",
                "protein_namespace",
                "protein_identifier",
                "protein_id",
                "site_sequence",
            ],
        ]
    )
    print(
        "protein_id present for all sites:",
        bool(
            dataset.site_metadata["protein_id"]
            .astype("string")
            .str.strip()
            .ne("")
            .all()
        ),
    )
    print(
        "Intensity state:",
        (
            f"scale={dataset.intensity_scale_state.label}, "
            f"quantitative_meaning={dataset.intensity_scale_state.quantity.value} "
            "(builder-established pass-through lane)"
        ),
    )
    print(
        "Processing state:",
        dataset.processing_state,
    )


if __name__ == "__main__":
    main()
