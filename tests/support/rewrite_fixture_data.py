from __future__ import annotations

import json
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
PUBLIC_WORKFLOW_REFERENCE = (
    ROOT / "tests_legacy" / "fixtures" / "public_workflow_reference"
)
REWRITE_PUBLIC_WORKFLOW_REFERENCE = (
    ROOT / "tests" / "fixtures" / "public_workflow_reference"
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
SIGNALOME_REWRITE_L6_ASSIGNMENTS_SELECTED = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE
    / "signalome_rewrite_l6_module_assignments_selected.csv"
)
SIGNALOME_REWRITE_L6_MODULES = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "signalome_rewrite_l6_modules.csv"
)
SIGNALOME_REWRITE_L6_NETWORK_NODES = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "signalome_rewrite_l6_network_nodes.csv"
)
SIGNALOME_REWRITE_L6_NETWORK_EDGES_SELECTED = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE
    / "signalome_rewrite_l6_network_edges_selected.csv"
)
SIGNALOME_REWRITE_L6_CONTRACT = (
    REWRITE_PUBLIC_WORKFLOW_REFERENCE / "signalome_rewrite_l6_contract.json"
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


@lru_cache(maxsize=1)
def load_signalome_rewrite_l6_contract() -> dict[str, object]:
    return json.loads(SIGNALOME_REWRITE_L6_CONTRACT.read_text(encoding="utf-8"))


@lru_cache(maxsize=1)
def load_signalome_rewrite_l6_module_assignments_selected() -> pd.DataFrame:
    frame = pd.read_csv(SIGNALOME_REWRITE_L6_ASSIGNMENTS_SELECTED, index_col=0)
    frame.index = pd.Index(frame.index.astype(str), name="site_id")
    return frame.astype(
        {
            "protein_id": str,
            "module_id": "int64",
            "top_kinase": str,
            "top_score": float,
            "top_kinase_tie_count": "int64",
            "top_kinase_is_ambiguous": bool,
        }
    )


@lru_cache(maxsize=1)
def load_signalome_rewrite_l6_modules() -> pd.DataFrame:
    frame = pd.read_csv(SIGNALOME_REWRITE_L6_MODULES, index_col=0)
    frame.index = pd.Index(frame.index.astype("int64"), name="module_id")
    frame.columns = pd.Index(frame.columns.astype(str), name="kinase")
    return frame.astype(float)


@lru_cache(maxsize=1)
def load_signalome_rewrite_l6_network_nodes() -> pd.DataFrame:
    frame = pd.read_csv(SIGNALOME_REWRITE_L6_NETWORK_NODES, index_col=0)
    frame.index = pd.Index(frame.index.astype(str), name="kinase")
    return frame.astype({"degree": "int64", "n_substrates": "int64"})


@lru_cache(maxsize=1)
def load_signalome_rewrite_l6_network_edges_selected() -> pd.DataFrame:
    return pd.read_csv(SIGNALOME_REWRITE_L6_NETWORK_EDGES_SELECTED).astype(
        {
            "source_kinase": str,
            "target_kinase": str,
            "correlation": float,
        }
    )


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
