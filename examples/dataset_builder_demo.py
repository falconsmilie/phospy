#!/usr/bin/env python3
"""Build an analysis-ready dataset in the preferred 1.5.0 first-run lane."""

from __future__ import annotations

import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset
from phospy.api import DatasetBuildRequest, Organism


def _example_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 0.7],
            "sample_b": [1.2, 0.8],
            "sample_c": [0.9, 0.75],
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
    return phospho, site_metadata


def build_demo_dataset() -> AnalysisReadyPhosphoDataset:
    phospho, site_metadata = _example_tables()
    request = DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
    )
    return AnalysisReadyDatasetBuilder().run(request)


def main() -> None:
    dataset = build_demo_dataset()
    print("Preferred 1.5.0 dataset builder lane")
    print("Input route: pandas DataFrame values")
    print("Bundled reference-compatible organism:", dataset.organism.value)
    print("Phospho shape:", dataset.phospho.shape)
    print("Site metadata columns:", list(dataset.site_metadata.columns))
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
        "Transformation state:",
        f"{dataset.transformation_state.label} (builder-established pass-through lane)",
    )


if __name__ == "__main__":
    main()
