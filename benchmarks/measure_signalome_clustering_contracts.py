#!/usr/bin/env python3
"""Benchmark signalome clustering backend contracts on stable fixtures.

Targets:
- backend parity visibility for `exact_python` and `scipy_hierarchical`
- exact-tree guard and full-correlation guard behavior
- sampled candidate-scoring activation behavior
"""

from __future__ import annotations

import json
import sys
import time
import tracemalloc
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


@dataclass(frozen=True, slots=True)
class _BackendFixture:
    name: str
    scoring_matrix: pd.DataFrame
    site_to_protein: pd.Series
    requested_module_count: int | None = None
    max_clusters: int = 8
    candidate_scoring_backend: str | None = None
    max_exact_cluster_tree_sites: int | None = None
    max_full_correlation_sites: int = 2000
    expect_error_substring: str | None = None


def _build_deterministic_scoring_matrix(
    *,
    n_sites: int,
    n_kinases: int,
    seed: int,
) -> pd.DataFrame:
    from tests.support.performance_contracts import deterministic_matrix

    matrix = deterministic_matrix(n_sites=n_sites, n_samples=n_kinases, seed=seed)
    matrix.columns = pd.Index(
        [f"KINASE_{index + 1:03d}" for index in range(int(matrix.shape[1]))],
        name="kinase",
    )
    return matrix


def _build_realistic_scoring_matrix(
    *,
    n_sites: int,
    n_kinases: int,
    seed: int,
) -> pd.DataFrame:
    from tests.support.performance_contracts import deterministic_site_ids

    rng = np.random.default_rng(seed)
    n_modules = max(4, min(20, int(n_sites) // 12))
    module_profiles = rng.normal(
        loc=10.0,
        scale=1.7,
        size=(n_modules, int(n_kinases)),
    )
    module_assignments = (np.arange(int(n_sites), dtype=int) * 5) % int(n_modules)
    noise = rng.normal(loc=0.0, scale=0.28, size=(int(n_sites), int(n_kinases)))
    values = np.round(module_profiles[module_assignments] + noise, decimals=6)
    site_ids = deterministic_site_ids(int(n_sites), start=40_000, gene_prefix="SIGSITE")
    return pd.DataFrame(
        values,
        index=site_ids,
        columns=pd.Index(
            [f"KINASE_{index + 1:03d}" for index in range(int(n_kinases))],
            name="kinase",
        ),
        dtype=float,
    )


def _build_site_to_protein(
    site_ids: pd.Index,
    *,
    sites_per_protein: int,
) -> pd.Series:
    if sites_per_protein < 1:
        raise ValueError("sites_per_protein must be >= 1")
    proteins = [
        f"PROT_{(index // int(sites_per_protein)) + 1:05d}"
        for index in range(int(site_ids.size))
    ]
    return pd.Series(
        proteins,
        index=pd.Index(site_ids.astype(str), name="site_id"),
        name="protein_id",
        dtype=str,
    )


def _measure_backend_call(
    func, *, warmup: bool = True
) -> tuple[object | None, Exception | None, float, float]:
    if warmup:
        try:
            func()
        except Exception:
            pass
    result: object | None = None
    error: Exception | None = None
    tracemalloc.start()
    started = time.perf_counter()
    try:
        result = func()
    except Exception as exc:  # pragma: no cover - exercised by manual benchmark runs
        error = exc
    runtime_seconds = time.perf_counter() - started
    _current_bytes, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    peak_mib = float(peak_bytes) / (1024.0 * 1024.0)
    return result, error, float(runtime_seconds), peak_mib


def _run_fixture_for_backend(
    *,
    fixture: _BackendFixture,
    backend_name: str,
) -> dict[str, object]:
    from phospy.errors.workflows import SignalomeScaleError
    from phospy.signalomes.clustering import (
        SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED,
        run_signalome_clustering_backend,
    )

    result, error, runtime_seconds, peak_mib = _measure_backend_call(
        lambda: run_signalome_clustering_backend(
            scoring_matrix=fixture.scoring_matrix,
            site_to_protein=fixture.site_to_protein,
            requested_module_count=fixture.requested_module_count,
            primary_threshold=0.45,
            fallback_threshold=0.15,
            max_clusters=fixture.max_clusters,
            candidate_scoring_backend=fixture.candidate_scoring_backend,
            max_exact_cluster_tree_sites=fixture.max_exact_cluster_tree_sites,
            max_full_correlation_sites=fixture.max_full_correlation_sites,
            backend_name=backend_name,
        ),
        warmup=True,
    )

    base_record: dict[str, object] = {
        "suite": "signalome_clustering_contracts_v2",
        "fixture_name": fixture.name,
        "backend_name": backend_name,
        "site_count": int(fixture.scoring_matrix.shape[0]),
        "kinase_count": int(fixture.scoring_matrix.shape[1]),
        "requested_candidate_scoring_backend": (
            "auto"
            if fixture.candidate_scoring_backend is None
            else str(fixture.candidate_scoring_backend)
        ),
        "requested_module_count": (
            None
            if fixture.requested_module_count is None
            else int(fixture.requested_module_count)
        ),
        "max_exact_cluster_tree_sites": (
            None
            if fixture.max_exact_cluster_tree_sites is None
            else int(fixture.max_exact_cluster_tree_sites)
        ),
        "max_full_correlation_sites": int(fixture.max_full_correlation_sites),
        "runtime_seconds": float(runtime_seconds),
        "peak_mib": float(peak_mib),
        "selected_module_count": None,
        "cluster_tree_backend_mode": None,
        "candidate_scoring_mode": None,
        "sampled_candidate_scoring_activated": None,
        "candidate_scoring_skipped": None,
        "exact_tree_construction_occurred": False,
        "status": "ok",
        "error_type": None,
        "error_contains_expected_token": None,
    }

    if error is not None:
        if not isinstance(error, SignalomeScaleError):
            raise error
        if fixture.expect_error_substring is None:
            raise RuntimeError(
                f"unexpected guard failure for fixture={fixture.name} backend={backend_name}: {error}"
            ) from error
        normalized_message = str(error).lower()
        base_record["status"] = "guard_error"
        base_record["error_type"] = type(error).__name__
        base_record["error_contains_expected_token"] = (
            fixture.expect_error_substring.lower() in normalized_message
        )
        return base_record

    if fixture.expect_error_substring is not None:
        raise RuntimeError(
            f"expected guard error for fixture={fixture.name} backend={backend_name}"
        )
    assert result is not None
    diagnostics = result.module_selection_diagnostics
    candidate_mode = str(result.candidate_scoring_mode)
    sampled_activated = bool(
        result.candidate_scoring_evaluated
        and candidate_mode == SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED
    )

    base_record.update(
        {
            "selected_module_count": int(diagnostics.selected_module_count),
            "cluster_tree_backend_mode": str(result.cluster_tree_backend),
            "candidate_scoring_mode": candidate_mode,
            "sampled_candidate_scoring_activated": sampled_activated,
            "candidate_scoring_skipped": bool(not result.candidate_scoring_evaluated),
            "exact_tree_construction_occurred": bool(result.exact_cluster_tree_built),
        }
    )
    return base_record


def main() -> None:
    from phospy.signalomes.clustering import (
        SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
        SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED,
        SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON,
        SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL,
    )

    backends = (
        SIGNALOME_CLUSTERING_BACKEND_EXACT_PYTHON,
        SIGNALOME_CLUSTERING_BACKEND_SCIPY_HIERARCHICAL,
    )

    small_deterministic = _build_deterministic_scoring_matrix(
        n_sites=48,
        n_kinases=14,
        seed=5101,
    )
    medium_realistic = _build_realistic_scoring_matrix(
        n_sites=240,
        n_kinases=24,
        seed=5102,
    )
    near_threshold = _build_realistic_scoring_matrix(
        n_sites=72,
        n_kinases=12,
        seed=5103,
    )
    full_guard = _build_realistic_scoring_matrix(
        n_sites=90,
        n_kinases=10,
        seed=5104,
    )
    sampled_mode = _build_realistic_scoring_matrix(
        n_sites=420,
        n_kinases=24,
        seed=5105,
    )

    fixtures = (
        _BackendFixture(
            name="signalome_small_deterministic_v1",
            scoring_matrix=small_deterministic,
            site_to_protein=_build_site_to_protein(
                small_deterministic.index,
                sites_per_protein=2,
            ),
            max_clusters=8,
            candidate_scoring_backend=SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
            max_exact_cluster_tree_sites=96,
            max_full_correlation_sites=96,
        ),
        _BackendFixture(
            name="signalome_medium_realistic_v1",
            scoring_matrix=medium_realistic,
            site_to_protein=_build_site_to_protein(
                medium_realistic.index,
                sites_per_protein=3,
            ),
            max_clusters=10,
            candidate_scoring_backend=SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
            max_exact_cluster_tree_sites=320,
            max_full_correlation_sites=320,
        ),
        _BackendFixture(
            name="signalome_near_exact_tree_limit_v1",
            scoring_matrix=near_threshold,
            site_to_protein=_build_site_to_protein(
                near_threshold.index,
                sites_per_protein=2,
            ),
            max_clusters=7,
            candidate_scoring_backend=SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
            max_exact_cluster_tree_sites=72,
            max_full_correlation_sites=72,
        ),
        _BackendFixture(
            name="signalome_full_correlation_guard_v1",
            scoring_matrix=full_guard,
            site_to_protein=_build_site_to_protein(
                full_guard.index,
                sites_per_protein=2,
            ),
            max_clusters=6,
            candidate_scoring_backend=SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
            max_exact_cluster_tree_sites=120,
            max_full_correlation_sites=80,
            expect_error_substring="full candidate-correlation scoring would evaluate",
        ),
        _BackendFixture(
            name="signalome_sampled_candidate_scoring_v1",
            scoring_matrix=sampled_mode,
            site_to_protein=_build_site_to_protein(
                sampled_mode.index,
                sites_per_protein=3,
            ),
            max_clusters=10,
            candidate_scoring_backend=SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED,
            max_exact_cluster_tree_sites=500,
            max_full_correlation_sites=180,
        ),
    )

    records: list[dict[str, object]] = []
    for fixture in fixtures:
        for backend_name in backends:
            records.append(
                _run_fixture_for_backend(fixture=fixture, backend_name=backend_name)
            )

    print("report_format=jsonl")
    print("benchmark_suite=signalome_clustering_contracts_v2")
    for record in records:
        print(
            json.dumps(
                record,
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=True,
            )
        )


if __name__ == "__main__":
    main()
