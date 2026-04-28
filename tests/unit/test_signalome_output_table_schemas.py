from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.errors import PhosPyValidationError
from phospy.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
)


def _valid_assignments() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "protein_id": ["P1", "P2"],
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
        index=pd.Index(["P1;S1;", "P2;S2;"], name="site_id"),
    )


def _valid_modules() -> pd.DataFrame:
    return pd.DataFrame(
        {"K1": [60.0, 0.0], "K2": [40.0, 0.0]},
        index=pd.Index([1, 2], name="module_id"),
    )


def _valid_edges() -> pd.DataFrame:
    return pd.DataFrame(
        {"source_kinase": ["K1"], "target_kinase": ["K2"], "correlation": [0.8]}
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


def test_signalome_assignments_reject_missing_required_columns() -> None:
    table = _valid_assignments().drop(columns=["top_kinase"])
    with pytest.raises(PhosPyValidationError, match="missing required columns"):
        SignalomeAssignments(table=table)


def test_signalome_assignments_reject_duplicate_site_index() -> None:
    table = _valid_assignments()
    table.index = pd.Index(["P1;S1;", "P1;S1;"], name="site_id")
    with pytest.raises(PhosPyValidationError, match="index must be unique"):
        SignalomeAssignments(table=table)


def test_signalome_assignments_reject_negative_module_id() -> None:
    table = _valid_assignments()
    table.loc["P1;S1;", "module_id"] = -1
    with pytest.raises(
        PhosPyValidationError, match="module_id must contain non-negative"
    ):
        SignalomeAssignments(table=table)


def test_signalome_assignments_reject_invalid_boolean_dtype() -> None:
    table = _valid_assignments().astype({"top_kinase_is_ambiguous": object})
    table.loc["P1;S1;", "top_kinase_is_ambiguous"] = "yes"
    with pytest.raises(PhosPyValidationError, match="must contain boolean values"):
        SignalomeAssignments(table=table)


def test_signalome_assignments_reject_missing_supported_top_score() -> None:
    table = _valid_assignments()
    table.loc["P1;S1;", "top_score"] = np.nan
    with pytest.raises(PhosPyValidationError, match="top_score must be finite"):
        SignalomeAssignments(table=table)


def test_signalome_assignments_reject_malformed_weight_shape() -> None:
    table = _valid_assignments()
    table.at["P1;S1;", "top_kinase_weights"] = ("K1", 1.0)
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
