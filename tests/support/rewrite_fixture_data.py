from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DatasetBuildRequest,
    Organism,
)

ROOT = Path(__file__).resolve().parents[2]
RAT_L6_PHOSPHO = (
    ROOT / "tests_legacy" / "fixtures" / "r_reference_l6" / "l6_phospho_matrix.csv"
)
RAT_L6_EXPECTED_PROFILE = (
    ROOT / "tests_legacy" / "fixtures" / "r_reference_l6" / "native_profile_scores.csv"
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


@lru_cache(maxsize=1)
def load_rat_l6_phospho() -> pd.DataFrame:
    return pd.read_csv(RAT_L6_PHOSPHO, index_col=0)


@lru_cache(maxsize=1)
def load_rat_l6_sequence_table() -> pd.Series:
    sequence_frame = pd.read_csv(RAT_L6_SITE_SEQUENCES)
    return sequence_frame.set_index("site_id").loc[:, "centralized_sequence"]


@lru_cache(maxsize=1)
def load_expected_profile_scores() -> pd.DataFrame:
    return pd.read_csv(RAT_L6_EXPECTED_PROFILE, index_col=0)


def site_metadata_for(phospho: pd.DataFrame) -> pd.DataFrame:
    split = phospho.index.to_series().astype(str).str.split(";", expand=True)
    site_sequences = load_rat_l6_sequence_table().reindex(phospho.index)
    if site_sequences.isna().any():
        missing = int(site_sequences.isna().sum())
        raise AssertionError(
            f"fixture missing site sequences for {missing} phosphosites"
        )
    return pd.DataFrame(
        {
            "gene_symbol": split.loc[:, 0].values,
            "site": split.loc[:, 1].values,
            "site_sequence": site_sequences.values,
        },
        index=phospho.index.copy(),
    )


def build_rat_l6_dataset(*, n_sites: int | None = 220) -> AnalysisReadyPhosphoDataset:
    phospho = load_rat_l6_phospho().copy(deep=True)
    if n_sites is not None:
        phospho = phospho.head(n_sites)
    request = DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata_for(phospho),
        organism=Organism.RAT,
    )
    return AnalysisReadyDatasetBuilder().run(request)
