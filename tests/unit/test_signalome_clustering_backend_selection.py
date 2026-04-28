from __future__ import annotations

import pandas as pd
import pytest

from phospy.signalomes.clustering import (
    SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON_VERSION,
    available_clustering_backends,
    resolve_clustering_backend,
    run_signalome_clustering_backend,
)


def _small_scoring_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        [[1.0, 0.1], [0.1, 1.0], [0.7, 0.6]],
        index=["P1;S1;", "P2;S2;", "P3;S3;"],
        columns=["K1", "K2"],
        dtype=float,
    )


def _small_site_to_protein() -> pd.Series:
    return pd.Series(
        ["P1", "P2", "P3"],
        index=pd.Index(["P1;S1;", "P2;S2;", "P3;S3;"], name="site_id"),
        name="protein_id",
        dtype=str,
    )


def test_backend_registry_exposes_exact_python_backend() -> None:
    assert SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON in available_clustering_backends()
    backend = resolve_clustering_backend(SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON)
    assert backend.name == SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON
    assert backend.version == SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON_VERSION


def test_backend_selection_rejects_unsupported_backend_name() -> None:
    with pytest.raises(ValueError, match="unsupported signalome clustering backend"):
        resolve_clustering_backend("not_a_backend")


def test_backend_result_surfaces_backend_and_limit_metadata() -> None:
    result = run_signalome_clustering_backend(
        scoring_matrix=_small_scoring_matrix(),
        site_to_protein=_small_site_to_protein(),
        requested_module_count=None,
        primary_threshold=0.5,
        fallback_threshold=0.1,
        max_clusters=5,
        max_exact_cluster_tree_sites=100,
        max_full_correlation_sites=50,
        backend_name=SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON,
    )

    assert result.backend_name == SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON
    assert result.backend_version == SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON_VERSION
    assert result.threshold_metadata == {
        "primary_threshold": 0.5,
        "fallback_threshold": 0.1,
    }
    assert result.limit_metadata == {
        "max_exact_cluster_tree_sites": 100,
        "max_full_correlation_sites": 50,
        "max_clusters": 5,
    }
