from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.science.signalomes.clustering import (
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    cluster_sites_with_diagnostics,
    derive_protein_modules,
    run_signalome_clustering_engine,
)
from phospy.science.signalomes.clustering import exact_python as legacy_exact

pytestmark = pytest.mark.parity


def _fixture_scoring_matrix() -> pd.DataFrame:
    values = np.asarray(
        [
            [1.0, 0.2, 0.4],
            [0.9, 0.1, 0.5],
            [-1.0, -0.2, -0.4],
            [-0.8, -0.1, -0.5],
            [0.2, 0.8, 0.3],
        ],
        dtype=float,
    )
    return pd.DataFrame(
        values,
        index=[f"P{idx};S{idx};" for idx in range(1, values.shape[0] + 1)],
        columns=["K1", "K2", "K3"],
    )


def _fixture_site_to_protein(site_index: pd.Index) -> pd.Series:
    proteins = ["P1", "P1", "P2", "P2", "P3"]
    return pd.Series(
        proteins,
        index=pd.Index(site_index.astype(str), name="site_id"),
        name="protein_id",
        dtype=str,
    )


def test_backend_facade_matches_legacy_exact_python_outputs() -> None:
    scoring_matrix = _fixture_scoring_matrix()
    legacy = legacy_exact.cluster_sites_with_diagnostics(
        scoring_matrix=scoring_matrix,
        requested_module_count=None,
        max_clusters=4,
    )
    facade = cluster_sites_with_diagnostics(
        scoring_matrix=scoring_matrix,
        requested_module_count=None,
        max_clusters=4,
    )

    pd.testing.assert_series_equal(facade.site_clusters, legacy.site_clusters)
    assert facade.module_selection_diagnostics == legacy.module_selection_diagnostics
    assert facade.candidate_scoring_mode == legacy.candidate_scoring_mode
    assert facade.candidate_scoring_evaluated == legacy.candidate_scoring_evaluated
    assert facade.candidate_scoring_skip_reason == legacy.candidate_scoring_skip_reason
    assert facade.exact_cluster_tree_built == legacy.exact_cluster_tree_built


def test_backend_protocol_matches_legacy_clustering_and_module_derivation() -> None:
    scoring_matrix = _fixture_scoring_matrix()
    site_to_protein = _fixture_site_to_protein(scoring_matrix.index)

    legacy = legacy_exact.cluster_sites_with_diagnostics(
        scoring_matrix=scoring_matrix,
        requested_module_count=None,
        max_clusters=4,
    )
    legacy_modules = derive_protein_modules(
        site_clusters=legacy.site_clusters,
        site_to_protein=site_to_protein,
    )

    backend = run_signalome_clustering_engine(
        scoring_matrix=scoring_matrix,
        site_to_protein=site_to_protein,
        requested_module_count=None,
        max_clusters=4,
        clustering_engine=SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    )

    pd.testing.assert_series_equal(backend.site_clusters, legacy.site_clusters)
    pd.testing.assert_series_equal(backend.protein_modules, legacy_modules)
    assert backend.module_selection_diagnostics == legacy.module_selection_diagnostics
    assert backend.selected_module_count == int(
        legacy.module_selection_diagnostics.selected_module_count
    )
