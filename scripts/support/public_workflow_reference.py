"""Shared dataset builders for active public-workflow reference generators."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from phospy.api import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DatasetBuildRequest,
    Organism,
)

ROOT = Path(__file__).resolve().parents[2]
RAT_L6_PHOSPHO = (
    ROOT
    / "tests"
    / "fixtures"
    / "rewrite_parity"
    / "r_reference_l6"
    / "l6_phospho_matrix.csv"
)
RAT_L6_SITE_SEQUENCES = (
    ROOT
    / "src"
    / "phospy"
    / "data"
    / "reference_bundles"
    / "rat"
    / "l6_native"
    / "site_sequences.csv"
)


def _site_metadata_for(phospho: pd.DataFrame) -> pd.DataFrame:
    split = phospho.index.to_series().astype(str).str.split(";", expand=True)
    sequence_frame = pd.read_csv(RAT_L6_SITE_SEQUENCES)
    site_sequences = sequence_frame.set_index("site_id").loc[:, "centralized_sequence"]
    site_sequences = site_sequences.reindex(phospho.index)

    if site_sequences.isna().any():
        missing = int(site_sequences.isna().sum())
        raise ValueError(
            f"fixture missing site sequences for {missing} phosphosites in {RAT_L6_SITE_SEQUENCES}"
        )

    return pd.DataFrame(
        {
            "gene_symbol": split.loc[:, 0].values,
            "site": split.loc[:, 1].values,
            "site_sequence": site_sequences.values,
            "protein_id": split.loc[:, 0].values,
        },
        index=phospho.index.copy(),
    )


def build_rat_l6_dataset(
    *,
    n_sites: int | None = 260,
) -> AnalysisReadyPhosphoDataset:
    phospho = pd.read_csv(RAT_L6_PHOSPHO, index_col=0)
    if n_sites is not None:
        phospho = phospho.head(n_sites)

    request = DatasetBuildRequest(
        phospho=phospho,
        site_metadata=_site_metadata_for(phospho),
        organism=Organism.RAT,
    )
    return AnalysisReadyDatasetBuilder().run(request)
