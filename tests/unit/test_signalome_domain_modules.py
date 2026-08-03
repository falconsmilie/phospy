from __future__ import annotations

import pandas as pd
import pytest

from phospy.api.configs import (
    SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
    SIGNALOME_KINASE_NETWORK_POLICY_SIGNED,
)
from phospy.errors import WorkflowStageError
from phospy.science.signalomes.expanded import build_expanded_signalome_table
from phospy.science.signalomes.modules import build_signalome_module_table
from phospy.science.signalomes.network import (
    build_kinase_network_with_diagnostics,
)


def test_domain_modules_weighted_top_requires_weight_column() -> None:
    module_assignments = pd.DataFrame(
        {
            "protein_group_id": ["P1", "P2"],
            "module_id": [1, 2],
        },
        index=pd.Index(["S1", "S2"], name="site_id"),
    )

    with pytest.raises(WorkflowStageError, match="top_kinase_weights"):
        build_signalome_module_table(
            module_assignments=module_assignments,
            kinase_substrates={"K1": (), "K2": ()},
            kinase_order=["K1", "K2"],
            assignment_policy=SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
        )


def test_domain_network_reports_non_finite_scores_without_edge_creation() -> None:
    downstream_scores = pd.DataFrame(
        {
            "K1": [1.0, float("inf")],
            "K2": [2.0, 3.0],
        },
        index=pd.Index(["S1", "S2"], name="site_id"),
    )
    edges, _, candidates, diagnostics = build_kinase_network_with_diagnostics(
        downstream_score_matrix=downstream_scores,
        kinase_order=["K1", "K2"],
        kinase_substrates={"K1": (), "K2": ()},
        threshold=0.0,
        network_policy=SIGNALOME_KINASE_NETWORK_POLICY_SIGNED,
    )

    assert edges.empty
    assert candidates.shape[0] == 1
    assert candidates.at[0, "correlation_status"] == "non_finite_values"
    assert pd.isna(candidates.at[0, "correlation"])
    assert diagnostics.non_finite_value_correlations == 1
    assert diagnostics.edges_skipped_non_finite_score == 1


def test_domain_expanded_requires_network_edge_columns() -> None:
    site_key = "S1"
    module_assignments = pd.DataFrame(
        {
            "site_key": [site_key],
            "display_id": ["S1"],
            "gene_symbol": ["P1"],
            "site": ["S1"],
            "protein_group_id": ["P1"],
            "protein_accession": [""],
            "isoform_id": [""],
            "module_id": [1],
            "top_kinase": ["K1"],
            "top_score": [0.9],
            "top_kinase_weights": [(("K1", 1.0),)],
        },
        index=pd.Index([site_key], name="site_key"),
    )
    signalome_modules = pd.DataFrame(
        {"K1": [100.0]},
        index=pd.Index([1], name="module_id"),
    )
    invalid_edges = pd.DataFrame({"source": ["K1"], "target": ["K2"]})

    with pytest.raises(
        WorkflowStageError,
        match="missing required columns for expanded signalome",
    ):
        build_expanded_signalome_table(
            module_assignments=module_assignments,
            signalome_modules=signalome_modules,
            kinase_network_edges=invalid_edges,
            kinase_substrates={"K1": ("S1",)},
            assignment_policy=SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
        )
