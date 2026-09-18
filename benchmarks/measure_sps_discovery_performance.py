"""Measure the production SPS discovery workflow and its stability hot path."""
# ruff: noqa: E402

from __future__ import annotations

import argparse
import platform
import sys
import time
from collections.abc import Sequence
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
for _path in (ROOT, SRC):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from phospy.science.batch_correction.sps_discovery import (
    _reference_stability_evidence,
)
from tests.support.performance_contracts import SPS_DISCOVERY_RUNTIME_SECONDS_MAX
from tests.support.sps_discovery_performance import (
    MODERATE_SPS_BENCHMARK_CONFIG,
    SpsDiscoveryBenchmarkConfig,
    build_sps_benchmark_references,
    run_sps_discovery_benchmark,
)

STRESS_SPS_BENCHMARK_CONFIG = SpsDiscoveryBenchmarkConfig(
    n_references=3,
    n_sites=50_000,
    n_samples=48,
    n_conditions=4,
)


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scale",
        choices=("moderate", "stress"),
        default="moderate",
        help="moderate=10,000x24; stress=50,000x48 per reference",
    )
    parser.add_argument(
        "--skip-stability-profile",
        action="store_true",
        help="skip the separate stability-helper timing pass",
    )
    return parser.parse_args(argv)


def _selected_config(scale: str) -> SpsDiscoveryBenchmarkConfig:
    if scale == "stress":
        return STRESS_SPS_BENCHMARK_CONFIG
    return MODERATE_SPS_BENCHMARK_CONFIG


def _dependency_version(distribution: str) -> str:
    try:
        return metadata.version(distribution)
    except metadata.PackageNotFoundError:
        return "unavailable"


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    config = _selected_config(args.scale)
    references = build_sps_benchmark_references(config)
    benchmark = run_sps_discovery_benchmark(references, config=config)
    boundaries = benchmark.result.provenance.selection_boundaries

    stability_seconds: float | None = None
    stability_valid_site_entries: int | None = None
    if not args.skip_stability_profile:
        started = time.perf_counter()
        evidence = tuple(
            _reference_stability_evidence(reference) for reference in references
        )
        stability_seconds = time.perf_counter() - started
        stability_valid_site_entries = sum(
            len(item.scores_by_site) for item in evidence
        )

    print(f"benchmark=sps_discovery_{args.scale}")
    print(f"reference_count={config.n_references}")
    print(f"phosphosite_count_per_reference={config.n_sites}")
    print(f"sample_count_per_reference={config.n_samples}")
    print(f"condition_count_per_reference={config.n_conditions}")
    print(f"replicates_per_condition={config.replicates_per_condition}")
    print(f"sites_entering_stability_calculation={config.n_sites}")
    print(f"stability_site_entries_total={config.n_references * config.n_sites}")
    print(f"sites_with_valid_stability={boundaries.sites_with_valid_stability}")
    print(f"sites_reaching_consensus_ranking={boundaries.sites_ranked}")
    print(f"selected_site_count={boundaries.sites_selected}")
    print(f"total_sps_discovery_seconds={benchmark.runtime_seconds:.6f}")
    print(f"release_runtime_threshold_seconds={SPS_DISCOVERY_RUNTIME_SECONDS_MAX:.6f}")
    print(
        "release_runtime_threshold_met="
        f"{str(benchmark.runtime_seconds < SPS_DISCOVERY_RUNTIME_SECONDS_MAX).lower()}"
    )
    if stability_seconds is None:
        print("stability_calculation_seconds=skipped")
        print("stability_valid_site_entries=skipped")
    else:
        print(f"stability_calculation_seconds={stability_seconds:.6f}")
        print(f"stability_valid_site_entries={stability_valid_site_entries}")
    print(f"python_version={platform.python_version()}")
    print(f"platform={platform.platform()}")
    print(f"numpy_version={_dependency_version('numpy')}")
    print(f"pandas_version={_dependency_version('pandas')}")
    print(f"scipy_version={_dependency_version('scipy')}")
    print(f"phospy_version={_dependency_version('phospy')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
