from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd
from pandas.api.types import (
    is_bool_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_object_dtype,
)

from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    KinasePredictionConfig,
    KinaseScoringConfig,
    Organism,
    ReferencePreset,
    SignalomeConfig,
    SignalomeWorkflow,
    SignalomeWorkflowRequest,
    SimpleKinaseWorkflow,
    SimpleKinaseWorkflowRequest,
)

ROOT = Path(__file__).resolve().parents[2]
RAT_L6_PHOSPHO = (
    ROOT / "tests_legacy" / "fixtures" / "r_reference_l6" / "l6_phospho_matrix.csv"
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
def _load_rat_l6_phospho() -> pd.DataFrame:
    return pd.read_csv(RAT_L6_PHOSPHO, index_col=0)


@lru_cache(maxsize=1)
def _load_rat_l6_sequence_table() -> pd.Series:
    sequence_frame = pd.read_csv(RAT_L6_SITE_SEQUENCES)
    return sequence_frame.set_index("site_id").loc[:, "centralized_sequence"]


def _site_metadata_for(phospho: pd.DataFrame) -> pd.DataFrame:
    split = phospho.index.to_series().astype(str).str.split(";", expand=True)
    site_sequences = _load_rat_l6_sequence_table().reindex(phospho.index)
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


def _build_dataset(*, n_sites: int = 260) -> object:
    phospho = _load_rat_l6_phospho().head(n_sites).copy(deep=True)
    request = DatasetBuildRequest(
        phospho=phospho,
        site_metadata=_site_metadata_for(phospho),
        organism=Organism.RAT,
    )
    return AnalysisReadyDatasetBuilder().run(request)


def test_signalome_workflow_runs_first_real_vertical_slice() -> None:
    dataset = _build_dataset(n_sites=260)
    kinase_result = SimpleKinaseWorkflow().run(
        SimpleKinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=1),
            prediction_config=KinasePredictionConfig(top_k=6, ensemble_size=12),
            activity_config=None,
        )
    )
    result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(signalome_cutoff=0.5),
        )
    )

    assignments = result.module_assignments.table
    assert not assignments.empty
    assert assignments.index.name == "site_id"
    assert {
        "protein_id",
        "module_id",
        "top_kinase",
        "top_score",
        "top_kinase_candidates",
        "top_kinase_tie_count",
        "top_kinase_is_ambiguous",
    }.issubset(set(assignments.columns))
    assert is_object_dtype(assignments.loc[:, "protein_id"])
    assert is_integer_dtype(assignments.loc[:, "module_id"])
    assert is_object_dtype(assignments.loc[:, "top_kinase"])
    assert is_float_dtype(assignments.loc[:, "top_score"])
    assert is_integer_dtype(assignments.loc[:, "top_kinase_tie_count"])
    assert is_bool_dtype(assignments.loc[:, "top_kinase_is_ambiguous"])

    modules = result.signalome_modules.table
    assert not modules.empty
    assert modules.index.name == "module_id"
    assert modules.columns.name == "kinase"
    assert is_float_dtype(modules.to_numpy(dtype=float))

    network_nodes = result.kinase_network.nodes
    assert network_nodes is not None
    assert not network_nodes.empty
    assert network_nodes.index.name == "kinase"
    assert {"degree", "n_substrates"} == set(network_nodes.columns)
    assert is_integer_dtype(network_nodes.loc[:, "degree"])
    assert is_integer_dtype(network_nodes.loc[:, "n_substrates"])

    network_edges = result.kinase_network.edges
    assert not network_edges.empty
    assert {"source_kinase", "target_kinase", "correlation"} == set(
        network_edges.columns
    )
    assert is_object_dtype(network_edges.loc[:, "source_kinase"])
    assert is_object_dtype(network_edges.loc[:, "target_kinase"])
    assert is_float_dtype(network_edges.loc[:, "correlation"])

    assert result.expanded_signalome is None
