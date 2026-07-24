from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.errors import PhosPyValidationError
from phospy.science.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
)
from phospy.science.tables.signalome import SignalomeAssignmentsTable
from tests.support.site_keys import (
    protein_site_key_index,
    site_key_index_from_display_ids,
)


def _valid_assignments() -> pd.DataFrame:
    display_ids = ["P1;S1;", "P2;S2;"]
    site_index = site_key_index_from_display_ids(display_ids)
    return pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": display_ids,
            "gene_symbol": ["P1", "P2"],
            "site": ["S1", "S2"],
            "protein_id": ["P1", "P2"],
            "protein_accession": ["", ""],
            "isoform_id": ["", ""],
            "module_id": [1, 0],
            "top_kinase": ["K1", "__UNSUPPORTED__"],
            "top_score": [0.9, np.nan],
            "top_kinase_candidates": [("K1",), ()],
            "top_kinase_weights": [(("K1", 1.0),), ()],
            "top_kinase_tie_count": [1, 0],
            "top_kinase_is_ambiguous": [False, False],
            "top_kinase_selection_policy": [
                "max_score_then_lexicographic_tiebreak",
                "no_support",
            ],
            "module_top_kinase": ["K1", "__UNSUPPORTED__"],
            "module_top_kinase_candidates": [("K1",), ()],
            "module_top_kinase_tie_count": [1, 0],
            "module_top_kinase_is_ambiguous": [False, False],
            "module_top_kinase_selection_policy": [
                "max_score_then_lexicographic_tiebreak",
                "no_support",
            ],
        },
        index=site_index,
    )


def _valid_modules() -> pd.DataFrame:
    return pd.DataFrame(
        {"K1": [60.0, 0.0], "K2": [40.0, 0.0]},
        index=pd.Index([1, 2], name="module_id"),
    )


def _valid_edges() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source_kinase": ["K1"],
            "target_kinase": ["K2"],
            "correlation": [0.8],
            "valid_observations": [4],
        }
    )


def _valid_nodes() -> pd.DataFrame:
    return pd.DataFrame(
        {"degree": [1, 1], "n_substrates": [3, 4]},
        index=pd.Index(["K1", "K2"], name="kinase"),
    )


def _valid_candidate_correlations() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source_kinase": ["K1", "K2"],
            "target_kinase": ["K2", "K3"],
            "correlation": [0.8, np.nan],
            "correlation_status": ["finite", "missing_values"],
            "valid_observations": [4, 1],
            "correlation_reason": [
                np.nan,
                "missing values reduced paired observations below minimum",
            ],
        }
    )


def test_signalome_assignments_table_accepts_valid_encoded_site_key_identity() -> None:
    table = _valid_assignments()

    wrapper = SignalomeAssignmentsTable(frame=table)

    assert wrapper.frame.index.name == "site_key"
    assert wrapper.frame.index.astype(str).tolist() == (
        wrapper.frame.loc[:, "site_key"].astype(str).tolist()
    )
    assert "display_id" in wrapper.frame.columns


def test_signalome_assignments_table_rejects_display_style_site_key_index() -> None:
    table = _valid_assignments()
    display_ids = ["AKT1;T308;", "GSK3B;S9;"]
    table.index = pd.Index(display_ids, name="site_key")
    table.loc[:, "site_key"] = display_ids
    table.loc[:, "display_id"] = display_ids

    with pytest.raises(PhosPyValidationError, match="display-indexed"):
        SignalomeAssignmentsTable(frame=table)


def test_signalome_assignments_table_rejects_arbitrary_site_key_strings() -> None:
    table = _valid_assignments()
    site_keys = ["not-a-site-key", "also-not-a-site-key"]
    table.index = pd.Index(site_keys, name="site_key")
    table.loc[:, "site_key"] = site_keys

    with pytest.raises(PhosPyValidationError, match="valid PhosPy site_key"):
        SignalomeAssignmentsTable(frame=table)


def test_signalome_assignments_table_rejects_valid_site_key_column_mismatch() -> None:
    table = _valid_assignments()
    table.loc[:, "site_key"] = list(reversed(table.index.astype(str).tolist()))

    with pytest.raises(PhosPyValidationError, match="site_key must exactly match"):
        SignalomeAssignmentsTable(frame=table)


def test_signalome_assignments_table_allows_duplicate_display_ids() -> None:
    table = _valid_assignments()
    site_index = protein_site_key_index(
        protein_identifiers=["P31749", "Q9Y243"],
        sites=["T308", "T308"],
    )
    table.index = site_index
    table.loc[:, "site_key"] = site_index.astype(str).tolist()
    table.loc[:, "display_id"] = ["AKT1;T308;", "AKT1;T308;"]
    table.loc[:, "gene_symbol"] = ["AKT1", "AKT1"]
    table.loc[:, "site"] = ["T308", "T308"]
    table.loc[:, "protein_id"] = ["P31749", "Q9Y243"]

    wrapper = SignalomeAssignmentsTable(frame=table)

    assert wrapper.frame.loc[:, "display_id"].duplicated().any()
    assert wrapper.frame.loc[:, "site_key"].is_unique


def test_signalome_assignments_reject_missing_required_columns() -> None:
    table = _valid_assignments().drop(columns=["top_kinase"])
    with pytest.raises(PhosPyValidationError, match="missing required columns"):
        SignalomeAssignments(table=table)


@pytest.mark.parametrize("column_name", ["site_key", "display_id"])
def test_signalome_assignments_reject_missing_identity_columns(
    column_name: str,
) -> None:
    table = _valid_assignments().drop(columns=[column_name])
    with pytest.raises(PhosPyValidationError, match=column_name):
        SignalomeAssignments(table=table)


def test_signalome_assignments_reject_duplicate_site_index() -> None:
    table = _valid_assignments()
    table.index = pd.Index([table.index[0], table.index[0]], name="site_key")
    with pytest.raises(PhosPyValidationError, match="index must be unique"):
        SignalomeAssignments(table=table)


def test_signalome_assignments_reject_negative_module_id() -> None:
    table = _valid_assignments()
    table.loc[table.index[0], "module_id"] = -1
    with pytest.raises(
        PhosPyValidationError, match="module_id must contain non-negative"
    ):
        SignalomeAssignments(table=table)


def test_signalome_assignments_reject_invalid_boolean_dtype() -> None:
    table = _valid_assignments().astype({"top_kinase_is_ambiguous": object})
    table.loc[table.index[0], "top_kinase_is_ambiguous"] = "yes"
    with pytest.raises(PhosPyValidationError, match="must contain boolean values"):
        SignalomeAssignments(table=table)


def test_signalome_assignments_reject_missing_supported_top_score() -> None:
    table = _valid_assignments()
    table.loc[table.index[0], "top_score"] = np.nan
    with pytest.raises(PhosPyValidationError, match="top_score must be finite"):
        SignalomeAssignments(table=table)


def test_signalome_assignments_reject_malformed_weight_shape() -> None:
    table = _valid_assignments()
    table.at[table.index[0], "top_kinase_weights"] = ("K1", 1.0)
    with pytest.raises(
        PhosPyValidationError, match="entries must be \\(kinase, weight\\) pairs"
    ):
        SignalomeAssignments(table=table)


def test_signalome_modules_reject_invalid_percentages() -> None:
    table = _valid_modules()
    table.loc[1, "K1"] = 150.0
    with pytest.raises(
        PhosPyValidationError, match="values must be between 0.0 and 100.0"
    ):
        SignalomeModules(table=table)


def test_signalome_modules_reject_invalid_row_totals() -> None:
    table = _valid_modules()
    table.loc[1, "K1"] = 30.0
    with pytest.raises(
        PhosPyValidationError, match="row totals must be approximately 0.0 or 100.0"
    ):
        SignalomeModules(table=table)


def test_kinase_network_rejects_malformed_edges() -> None:
    edges = _valid_edges().rename(columns={"source_kinase": "source"})
    with pytest.raises(PhosPyValidationError, match="missing required columns"):
        KinaseNetwork(edges=edges)


def test_kinase_network_edges_require_valid_observation_counts() -> None:
    edges = _valid_edges().drop(columns=["valid_observations"])
    with pytest.raises(PhosPyValidationError, match="missing required columns"):
        KinaseNetwork(edges=edges)


def test_kinase_network_edges_reject_negative_valid_observation_counts() -> None:
    edges = _valid_edges()
    edges.loc[0, "valid_observations"] = -1
    with pytest.raises(PhosPyValidationError, match="non-negative integer"):
        KinaseNetwork(edges=edges)


def test_kinase_network_rejects_out_of_bounds_edge_correlation() -> None:
    edges = _valid_edges()
    edges.loc[0, "correlation"] = 1.5
    with pytest.raises(PhosPyValidationError, match="must be between -1.0 and 1.0"):
        KinaseNetwork(edges=edges)


def test_kinase_network_rejects_malformed_candidate_correlation_columns() -> None:
    candidates = _valid_candidate_correlations().drop(columns=["valid_observations"])
    with pytest.raises(PhosPyValidationError, match="missing required columns"):
        KinaseNetwork(
            edges=_valid_edges(),
            nodes=_valid_nodes(),
            candidate_correlations=candidates,
        )


def test_kinase_network_rejects_invalid_candidate_correlation_semantics() -> None:
    candidates = _valid_candidate_correlations()
    candidates.loc[0, "correlation_status"] = "finite"
    candidates.loc[0, "correlation"] = np.nan
    with pytest.raises(
        PhosPyValidationError, match="must be present when correlation_status='finite'"
    ):
        KinaseNetwork(
            edges=_valid_edges(),
            nodes=_valid_nodes(),
            candidate_correlations=candidates,
        )
