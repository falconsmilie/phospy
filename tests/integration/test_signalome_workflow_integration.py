from __future__ import annotations

import pytest
from pandas.api.types import (
    is_bool_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_object_dtype,
    is_string_dtype,
)

from phospy import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    ReferencePreset,
    SignalomeConfig,
    SignalomeWorkflow,
    SignalomeWorkflowRequest,
    SimpleKinaseWorkflowRequest,
)
from tests.support.rewrite_fixture_data import build_rat_l6_dataset

pytestmark = pytest.mark.integration


def _is_text_dtype(values: object) -> bool:
    return is_object_dtype(values) or is_string_dtype(values)


def test_signalome_workflow_runs_dataset_to_kinase_to_signalome_path() -> None:
    dataset = build_rat_l6_dataset(n_sites=260)
    kinase_result = KinaseWorkflow().run(
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
    assert _is_text_dtype(assignments.loc[:, "protein_id"])
    assert is_integer_dtype(assignments.loc[:, "module_id"])
    assert _is_text_dtype(assignments.loc[:, "top_kinase"])
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
    assert _is_text_dtype(network_edges.loc[:, "source_kinase"])
    assert _is_text_dtype(network_edges.loc[:, "target_kinase"])
    assert is_float_dtype(network_edges.loc[:, "correlation"])

    assert result.expanded_signalome is None
