from __future__ import annotations

import pytest

from tests.support.performance_contracts import (
    SPS_DISCOVERY_RUNTIME_SECONDS_MAX,
)
from tests.support.sps_discovery_performance import (
    MODERATE_SPS_BENCHMARK_CONFIG,
    build_sps_benchmark_references,
    run_sps_discovery_benchmark,
)

pytestmark = [pytest.mark.performance, pytest.mark.release_gate]


def test_sps_discovery_moderate_scale_production_path() -> None:
    config = MODERATE_SPS_BENCHMARK_CONFIG
    references = build_sps_benchmark_references(config)

    benchmark = run_sps_discovery_benchmark(references, config=config)

    boundaries = benchmark.result.provenance.selection_boundaries
    assert len(references) == config.n_references
    assert all(
        reference.intensities.shape == (config.n_sites, config.n_samples)
        for reference in references
    )
    assert boundaries.total_unique_sites == config.n_sites
    assert boundaries.sites_meeting_dataset_overlap == config.n_sites
    assert boundaries.sites_with_valid_stability == config.n_sites
    assert boundaries.sites_ranked == config.n_sites
    assert boundaries.sites_selected == config.top_n
    assert benchmark.runtime_seconds < SPS_DISCOVERY_RUNTIME_SECONDS_MAX
