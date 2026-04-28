from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.signalomes.clustering import (
    SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON_VERSION,
    SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL,
    SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL_VERSION,
    available_clustering_backends,
    resolve_clustering_backend,
    run_signalome_clustering_backend,
)


def _small_scoring_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        [
            [0.95, 0.1, 0.2],
            [0.1, 0.95, 0.3],
            [0.85, 0.7, 0.25],
            [0.2, 0.1, 0.9],
        ],
        index=["P1;S1;", "P2;S2;", "P3;S3;", "P4;S4;"],
        columns=["K1", "K2", "K3"],
        dtype=float,
    )


def _small_site_to_protein() -> pd.Series:
    return pd.Series(
        ["P1", "P2", "P3", "P4"],
        index=pd.Index(["P1;S1;", "P2;S2;", "P3;S3;", "P4;S4;"], name="site_id"),
        name="protein_id",
        dtype=str,
    )


def _required_backend_diagnostic_keys() -> set[str]:
    return {
        "backend_name",
        "uses_scipy",
        "linkage_method",
        "distance_metric",
        "selected_module_count",
        "input_site_count",
        "exact_tree_path_used",
    }


def test_backend_registry_exposes_exact_and_scipy_backends() -> None:
    names = set(available_clustering_backends())
    assert SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON in names
    assert SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL in names

    exact = resolve_clustering_backend(SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON)
    assert exact.name == SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON
    assert exact.version == SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON_VERSION

    scipy_backend = resolve_clustering_backend(
        SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL
    )
    assert scipy_backend.name == SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL
    assert (
        scipy_backend.version == SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL_VERSION
    )


def test_backend_selection_rejects_unsupported_backend_name() -> None:
    with pytest.raises(ValueError, match="unsupported signalome clustering backend"):
        resolve_clustering_backend("not_a_backend")


def test_exact_backend_result_surfaces_limit_threshold_and_backend_diagnostics() -> (
    None
):
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
    assert result.backend_diagnostics is not None
    assert _required_backend_diagnostic_keys() <= set(result.backend_diagnostics)
    assert result.backend_diagnostics["backend_name"] == result.backend_name
    assert result.backend_diagnostics["uses_scipy"] is False


def test_scipy_backend_matches_exact_backend_for_small_deterministic_fixture() -> None:
    scoring_matrix = _small_scoring_matrix()
    site_to_protein = _small_site_to_protein()
    exact = run_signalome_clustering_backend(
        scoring_matrix=scoring_matrix,
        site_to_protein=site_to_protein,
        requested_module_count=None,
        primary_threshold=0.5,
        fallback_threshold=0.1,
        max_clusters=5,
        backend_name=SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON,
    )
    scipy_backend = run_signalome_clustering_backend(
        scoring_matrix=scoring_matrix,
        site_to_protein=site_to_protein,
        requested_module_count=None,
        primary_threshold=0.5,
        fallback_threshold=0.1,
        max_clusters=5,
        backend_name=SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL,
    )

    pdt.assert_series_equal(scipy_backend.site_clusters, exact.site_clusters)
    pdt.assert_series_equal(scipy_backend.protein_modules, exact.protein_modules)
    assert (
        scipy_backend.module_selection_diagnostics == exact.module_selection_diagnostics
    )
    assert scipy_backend.selected_module_count == exact.selected_module_count
    assert scipy_backend.backend_diagnostics is not None
    assert scipy_backend.backend_diagnostics["uses_scipy"] is True
    assert scipy_backend.backend_diagnostics["backend_name"] == (
        SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL
    )


def test_scipy_backend_is_stable_across_repeated_runs() -> None:
    kwargs = {
        "scoring_matrix": _small_scoring_matrix(),
        "site_to_protein": _small_site_to_protein(),
        "requested_module_count": None,
        "primary_threshold": 0.5,
        "fallback_threshold": 0.1,
        "max_clusters": 5,
        "backend_name": SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL,
    }
    run_one = run_signalome_clustering_backend(**kwargs)
    run_two = run_signalome_clustering_backend(**kwargs)

    pdt.assert_series_equal(run_one.site_clusters, run_two.site_clusters)
    pdt.assert_series_equal(run_one.protein_modules, run_two.protein_modules)
    assert run_one.module_selection_diagnostics == run_two.module_selection_diagnostics
    assert run_one.selected_module_count == run_two.selected_module_count
    assert run_one.backend_diagnostics == run_two.backend_diagnostics
