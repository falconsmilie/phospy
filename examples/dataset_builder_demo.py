#!/usr/bin/env python3
"""Build analysis-ready datasets through both supported builder input routes."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DatasetBuildRequest,
    Organism,
)


def _example_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
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
    return phospho, site_metadata


def build_from_dataframes() -> AnalysisReadyPhosphoDataset:
    phospho, site_metadata = _example_tables()
    request = DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
    )
    return AnalysisReadyDatasetBuilder().run(request)


def build_from_file_paths() -> AnalysisReadyPhosphoDataset:
    phospho, site_metadata = _example_tables()
    with TemporaryDirectory(prefix="phospy-builder-demo-") as tmp_dir:
        root = Path(tmp_dir)
        phospho_path = root / "phospho.csv"
        site_metadata_path = root / "site_metadata.csv"
        phospho.to_csv(phospho_path)
        site_metadata.to_csv(site_metadata_path)
        request = DatasetBuildRequest(
            phospho=phospho_path,
            site_metadata=str(site_metadata_path),
            organism=Organism.RAT,
        )
        return AnalysisReadyDatasetBuilder().run(request)


def main() -> None:
    df_dataset = build_from_dataframes()
    path_dataset = build_from_file_paths()
    print("Dataset builder demo")
    print("DataFrame route phospho shape:", df_dataset.phospho.shape)
    print("File-path route phospho shape:", path_dataset.phospho.shape)
    print("Site metadata columns:", list(path_dataset.site_metadata.columns))
    print(
        "Organism:",
        None if path_dataset.organism is None else path_dataset.organism.value,
    )
    print(
        "Transformation state:",
        f"{path_dataset.transformation_state.label} (builder pass-through lane)",
    )


if __name__ == "__main__":
    main()
