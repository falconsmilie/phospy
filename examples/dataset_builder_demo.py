#!/usr/bin/env python3
"""Build an analysis-ready dataset with the supported rewrite builder API."""

from __future__ import annotations

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DatasetBuildRequest,
    Organism,
)


def run_demo() -> AnalysisReadyPhosphoDataset:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 0.7],
            "sample_b": [1.2, 0.8],
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
    request = DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
    )
    return AnalysisReadyDatasetBuilder().run(request)


def main() -> None:
    dataset = run_demo()
    print("Dataset builder demo")
    print("Phospho shape:", dataset.phospho.shape)
    print("Site metadata columns:", list(dataset.site_metadata.columns))
    print("Organism:", None if dataset.organism is None else dataset.organism.value)
    print("Transformation state:", dataset.transformation_state.label)


if __name__ == "__main__":
    main()
