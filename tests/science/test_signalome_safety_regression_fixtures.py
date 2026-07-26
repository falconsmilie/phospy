from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.errors.validation import ContractValidationError
from phospy.io.bundles.signalome import SignalomeWorkflowConfigSnapshot
from phospy.science.signalomes.clustering import cluster_sites_with_diagnostics
from phospy.science.signalomes.clustering.tree_building import (
    prepare_signalome_clustering_matrix,
)
from phospy.science.signalomes.science import build_kinase_network_with_diagnostics
from tests.support.signalome_config import build_signalome_config

pytestmark = pytest.mark.release_gate

FIXTURE_DIR = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "release_validation_regression"
    / "signalome_safety"
)


def _manifest() -> dict[str, Any]:
    return json.loads((FIXTURE_DIR / "MANIFEST.json").read_text(encoding="utf-8"))


def _read_matrix(name: str) -> pd.DataFrame:
    return pd.read_csv(FIXTURE_DIR / name).set_index("site_id")


def test_signalome_safety_fixture_manifest_hashes_match_files() -> None:
    manifest = _manifest()

    assert manifest["classification"] == "regression"
    assert manifest["fixture_family"] == "signalome_safety"
    assert "not external parity" in manifest["source_policy"]
    assert manifest["seed"] == 20260724

    for file_entry in manifest["files"]:
        path = FIXTURE_DIR / str(file_entry["relative_path"])
        assert path.is_file()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == file_entry["sha256"]


def test_signalome_two_paired_observations_do_not_create_edge() -> None:
    edges, _, candidates, diagnostics = build_kinase_network_with_diagnostics(
        downstream_score_matrix=_read_matrix("network_two_observations.csv"),
        kinase_order=["K1", "K2"],
        kinase_substrates={"K1": (), "K2": ()},
        threshold=0.0,
        network_policy="signed",
        min_paired_observations=3,
    )

    assert edges.empty
    assert candidates.at[0, "valid_observations"] == 2
    assert candidates.at[0, "correlation_status"] == "insufficient_observations"
    assert diagnostics.edges_skipped_insufficient_paired_observations == 1


def test_signalome_threshold_boundary_three_and_default_five_observation_cases() -> (
    None
):
    three_edges, _, three_candidates, three_diagnostics = (
        build_kinase_network_with_diagnostics(
            downstream_score_matrix=_read_matrix("network_three_observations.csv"),
            kinase_order=["K1", "K2"],
            kinase_substrates={"K1": (), "K2": ()},
            threshold=0.9,
            network_policy="signed",
            min_paired_observations=3,
        )
    )
    five_edges, _, five_candidates, five_diagnostics = (
        build_kinase_network_with_diagnostics(
            downstream_score_matrix=_read_matrix("network_five_observations.csv"),
            kinase_order=["K1", "K2"],
            kinase_substrates={"K1": (), "K2": ()},
            threshold=0.9,
            network_policy="signed",
        )
    )

    assert three_candidates.at[0, "valid_observations"] == 3
    assert three_candidates.at[0, "correlation_status"] == "finite"
    assert three_edges.to_dict("records") == [
        {
            "source_kinase": "K1",
            "target_kinase": "K2",
            "correlation": pytest.approx(1.0),
            "valid_observations": 3,
        }
    ]
    assert three_diagnostics.edges_created == 1
    assert five_candidates.at[0, "valid_observations"] == 5
    assert five_candidates.at[0, "correlation_status"] == "finite"
    assert int(five_edges.shape[0]) == 1
    assert five_diagnostics.edges_created == 1


def test_signalome_missing_dimensions_are_dropped_or_median_imputed() -> None:
    prepared = prepare_signalome_clustering_matrix(
        _read_matrix("clustering_missing_dimensions.csv")
    )

    assert prepared.dropped_fully_missing_column_labels == ("K3_all_missing",)
    assert prepared.imputed_value_counts_by_column == {
        "K1": 0,
        "K2": 0,
        "K4_partial_missing": 2,
    }
    assert prepared.imputed_value_count == 2
    assert np.isfinite(prepared.values).all()
    pdt.assert_series_equal(
        prepared.prepared_matrix.loc[:, "K4_partial_missing"],
        pd.Series(
            [10.0, 0.0, 12.0, -10.0, 0.0, -12.0],
            index=prepared.prepared_matrix.index,
            name="K4_partial_missing",
        ),
    )


def test_signalome_clustering_is_invariant_to_added_all_missing_dimension() -> None:
    with_all_missing = _read_matrix("clustering_missing_dimensions.csv")
    baseline = with_all_missing.drop(columns=["K3_all_missing"])

    baseline_result = cluster_sites_with_diagnostics(
        scoring_matrix=baseline,
        requested_module_count=2,
    )
    observed = cluster_sites_with_diagnostics(
        scoring_matrix=with_all_missing,
        requested_module_count=2,
    )

    pdt.assert_series_equal(observed.site_clusters, baseline_result.site_clusters)
    assert (
        observed.clustering_preparation_diagnostics.dropped_fully_missing_dimension_labels
        == ("K3_all_missing",)
    )
    assert (
        baseline_result.clustering_preparation_diagnostics.dropped_fully_missing_dimension_labels
        == ()
    )
    assert (
        observed.module_selection_diagnostics.selected_module_count
        == baseline_result.module_selection_diagnostics.selected_module_count
    )


def test_signalome_historical_threshold2_reconstructs_but_new_execution_rejects() -> (
    None
):
    snapshot_payload = json.loads(
        (FIXTURE_DIR / "historical_threshold2_config.json").read_text(encoding="utf-8")
    )
    snapshot = SignalomeWorkflowConfigSnapshot.from_payload(snapshot_payload)

    assert snapshot.signalome_config.output.network_min_paired_finite_observations == 2
    with pytest.raises(
        ContractValidationError,
        match="network_min_paired_finite_observations",
    ):
        build_signalome_config(
            network_policy="signed",
            network_min_paired_finite_observations=2,
        )
