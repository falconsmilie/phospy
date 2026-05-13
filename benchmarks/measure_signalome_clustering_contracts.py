#!/usr/bin/env python3
"""Benchmark signalome clustering contracts with explicit guardrail scenarios.

Targets:
- exact clustering below `max_exact_tree_sites`
- exact-tree guard failure above `max_exact_tree_sites`
- candidate scoring behavior under `full` and `sampled` policies
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
class _Fixture:
    name: str
    scoring_matrix: pd.DataFrame
    site_to_protein: pd.Series
    candidate_scoring_policy: str
    max_exact_tree_sites: int
    max_full_candidate_scoring_sites: int
    max_clusters: int = 10


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
    return pd.DataFrame(
        values,
        index=deterministic_site_ids(int(n_sites), start=40_000, gene_prefix="SIGSITE"),
        columns=pd.Index(
            [f"KINASE_{index + 1:03d}" for index in range(int(n_kinases))],
            name="kinase",
        ),
        dtype=float,
    )


def _build_site_to_protein(site_ids: pd.Index, *, sites_per_protein: int) -> pd.Series:
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


def _measure(
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
    except Exception as exc:  # pragma: no cover - exercised by benchmark runs
        error = exc
    elapsed = time.perf_counter() - started
    _current_bytes, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result, error, float(elapsed), float(peak_bytes) / (1024.0 * 1024.0)


def _run_fixture(*, fixture: _Fixture) -> dict[str, object]:
    from phospy.errors.workflows import SignalomeScaleError
    from phospy.science.signalomes.clustering import run_signalome_clustering_engine

    result, error, runtime_seconds, peak_mib = _measure(
        lambda: run_signalome_clustering_engine(
            scoring_matrix=fixture.scoring_matrix,
            site_to_protein=fixture.site_to_protein,
            requested_module_count=None,
            primary_threshold=0.45,
            fallback_threshold=0.15,
            max_clusters=fixture.max_clusters,
            candidate_scoring_policy=fixture.candidate_scoring_policy,
            max_exact_tree_sites=fixture.max_exact_tree_sites,
            max_full_candidate_scoring_sites=fixture.max_full_candidate_scoring_sites,
        ),
        warmup=True,
    )

    record: dict[str, object] = {
        "suite": "signalome_clustering_contracts_v3",
        "fixture_name": fixture.name,
        "site_count": int(fixture.scoring_matrix.shape[0]),
        "kinase_count": int(fixture.scoring_matrix.shape[1]),
        "candidate_scoring_policy": fixture.candidate_scoring_policy,
        "max_exact_tree_sites": int(fixture.max_exact_tree_sites),
        "max_full_candidate_scoring_sites": int(
            fixture.max_full_candidate_scoring_sites
        ),
        "runtime_seconds": float(runtime_seconds),
        "peak_mib": float(peak_mib),
        "status": "ok",
        "candidate_scoring_mode": None,
        "selected_module_count": None,
        "exact_tree_construction_occurred": False,
        "error_type": None,
    }

    if error is not None:
        if not isinstance(error, SignalomeScaleError):
            raise error
        record["status"] = "guard_error"
        record["error_type"] = type(error).__name__
        record["error_message"] = str(error)
        return record

    assert result is not None
    diagnostics = result.module_selection_diagnostics
    record.update(
        {
            "candidate_scoring_mode": str(result.candidate_scoring_mode),
            "selected_module_count": int(diagnostics.selected_module_count),
            "exact_tree_construction_occurred": bool(result.exact_cluster_tree_built),
            "candidate_scoring_evaluated": bool(result.candidate_scoring_evaluated),
            "candidate_scoring_skip_reason": result.candidate_scoring_skip_reason,
            "tree_implementation": str(result.tree_implementation),
        }
    )
    return record


def main() -> None:
    from phospy.science.signalomes.clustering import (
        SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    )

    below_guard_matrix = _build_realistic_scoring_matrix(
        n_sites=1_800,
        n_kinases=60,
        seed=5101,
    )
    above_guard_matrix = _build_realistic_scoring_matrix(
        n_sites=2_120,
        n_kinases=60,
        seed=5102,
    )

    fixtures = (
        _Fixture(
            name="signalome_exact_below_guardrail_v1",
            scoring_matrix=below_guard_matrix,
            site_to_protein=_build_site_to_protein(
                below_guard_matrix.index, sites_per_protein=3
            ),
            candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
            max_exact_tree_sites=2_000,
            max_full_candidate_scoring_sites=2_000,
            max_clusters=10,
        ),
        _Fixture(
            name="signalome_exact_guardrail_trigger_v1",
            scoring_matrix=above_guard_matrix,
            site_to_protein=_build_site_to_protein(
                above_guard_matrix.index, sites_per_protein=3
            ),
            candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
            max_exact_tree_sites=2_000,
            max_full_candidate_scoring_sites=2_500,
            max_clusters=10,
        ),
        _Fixture(
            name="signalome_candidate_scoring_full_v1",
            scoring_matrix=below_guard_matrix,
            site_to_protein=_build_site_to_protein(
                below_guard_matrix.index, sites_per_protein=3
            ),
            candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
            max_exact_tree_sites=2_500,
            max_full_candidate_scoring_sites=2_500,
            max_clusters=10,
        ),
        _Fixture(
            name="signalome_candidate_scoring_sampled_v1",
            scoring_matrix=below_guard_matrix,
            site_to_protein=_build_site_to_protein(
                below_guard_matrix.index, sites_per_protein=3
            ),
            candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
            max_exact_tree_sites=2_500,
            max_full_candidate_scoring_sites=2_500,
            max_clusters=10,
        ),
    )

    records = [_run_fixture(fixture=fixture) for fixture in fixtures]

    print("report_format=jsonl")
    print("benchmark_suite=signalome_clustering_contracts_v3")
    for record in records:
        print(
            json.dumps(record, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        )


if __name__ == "__main__":
    main()
