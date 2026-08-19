#!/usr/bin/env python3
"""Measure duplicate-correlation estimator and GLS hot-path performance.

This is an explicitly invoked local benchmark. It is intentionally bounded and
machine-observational: it is suitable for detecting accidental algorithmic
regressions during development, but it is not a release-scale gate.
"""
# ruff: noqa: E402

from __future__ import annotations

import argparse
import json
import sys
import time
import tracemalloc
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from phospy.science.differential.compound_symmetry_gls import (
    fit_duplicate_correlation_gls,
)
from phospy.science.differential.duplicate_correlation import (
    estimate_duplicate_correlation_reml_consensus,
)

DOCUMENTED_MAIN_METRIC_KEYS = (
    "sites",
    "samples",
    "blocks",
    "observations_per_block",
    "missing_fraction",
    "correlation_estimation_seconds",
    "gls_seconds",
    "total_workflow_seconds",
    "estimated_feature_correlation_count",
    "consensus_correlation",
    "unique_missingness_masks",
    "factorization_cache_size",
    "factorization_cache_hit_count",
    "factorization_cache_reuse_count",
    "analysis_matrix_memory_bytes",
    "covariance_matrix_memory_bytes",
    "approximate_matrix_memory_bytes",
    "peak_tracemalloc_bytes",
    "approximate_peak_matrix_memory_bytes",
)


@dataclass(frozen=True, slots=True)
class BenchmarkConfig:
    sites: int = 1_200
    blocks: int = 8
    observations_per_block: int = 3
    missing_fraction: float = 0.08
    seed: int = 20_260_818
    block_correlation: float = 0.30
    effect_scale: float = 0.25


@dataclass(frozen=True, slots=True)
class BenchmarkInputs:
    matrix: np.ndarray
    design: np.ndarray
    block_ids: tuple[str, ...]
    feature_ids: tuple[str, ...]
    coefficient_names: tuple[str, ...]
    contrasts: np.ndarray
    contrast_names: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class BenchmarkResult:
    config: BenchmarkConfig
    metrics: dict[str, float | int]


def default_config() -> BenchmarkConfig:
    return BenchmarkConfig()


def build_inputs(config: BenchmarkConfig) -> BenchmarkInputs:
    if config.sites < 1:
        raise ValueError("sites must be positive")
    if config.blocks < 3:
        raise ValueError("blocks must be at least 3")
    if config.observations_per_block < 2:
        raise ValueError("observations_per_block must be at least 2")
    if not 0.0 <= config.missing_fraction < 0.5:
        raise ValueError("missing_fraction must be in [0.0, 0.5)")
    if not 0.0 <= config.block_correlation < 1.0:
        raise ValueError("block_correlation must be in [0.0, 1.0)")

    rng = np.random.default_rng(config.seed)
    samples = config.blocks * config.observations_per_block
    condition_count = config.observations_per_block
    condition_positions = np.tile(
        np.arange(condition_count, dtype=np.int64),
        config.blocks,
    )
    block_positions = np.repeat(
        np.arange(config.blocks, dtype=np.int64),
        config.observations_per_block,
    )

    design = np.zeros((samples, condition_count), dtype=np.float64)
    design[np.arange(samples, dtype=np.int64), condition_positions] = 1.0

    residual_variance = 1.0 - config.block_correlation
    block_variance = config.block_correlation
    feature_baseline = rng.normal(12.0, 0.8, size=(config.sites, 1))
    condition_effects = rng.normal(
        0.0,
        config.effect_scale,
        size=(config.sites, condition_count),
    )
    block_effects = rng.normal(
        0.0,
        float(np.sqrt(block_variance)),
        size=(config.sites, config.blocks),
    )
    residual = rng.normal(
        0.0,
        float(np.sqrt(residual_variance)),
        size=(config.sites, samples),
    )
    matrix = (
        feature_baseline
        + condition_effects[:, condition_positions]
        + block_effects[:, block_positions]
        + residual
    ).astype(np.float64)
    _apply_repeated_missingness_patterns(
        matrix,
        config=config,
        rng=rng,
    )

    contrast_count = max(1, condition_count - 1)
    contrasts = np.zeros((condition_count, contrast_count), dtype=np.float64)
    contrast_names: list[str] = []
    for contrast_position in range(contrast_count):
        numerator = contrast_position + 1
        contrasts[0, contrast_position] = -1.0
        contrasts[numerator, contrast_position] = 1.0
        contrast_names.append(f"C{numerator + 1}_vs_C1")

    return BenchmarkInputs(
        matrix=matrix,
        design=design,
        block_ids=tuple(
            f"block_{position + 1:03d}" for position in block_positions.tolist()
        ),
        feature_ids=tuple(
            f"site_{position + 1:06d}" for position in range(config.sites)
        ),
        coefficient_names=tuple(
            f"C{position + 1}" for position in range(condition_count)
        ),
        contrasts=contrasts,
        contrast_names=tuple(contrast_names),
    )


def _apply_repeated_missingness_patterns(
    matrix: np.ndarray,
    *,
    config: BenchmarkConfig,
    rng: np.random.Generator,
) -> None:
    samples = int(matrix.shape[1])
    if config.missing_fraction == 0.0:
        return

    pattern_count = min(16, max(2, samples))
    missing_per_pattern = max(1, int(round(config.missing_fraction * samples)))
    patterns: list[np.ndarray] = [np.ones(samples, dtype=bool)]
    for pattern_position in range(1, pattern_count):
        mask = np.ones(samples, dtype=bool)
        start = (pattern_position * missing_per_pattern) % samples
        missing_positions = (start + np.arange(missing_per_pattern)) % samples
        mask[missing_positions] = False
        patterns.append(mask)

    for row_position in range(matrix.shape[0]):
        pattern = patterns[int(row_position % len(patterns))]
        if not bool(pattern.all()):
            matrix[row_position, ~pattern] = np.nan

    random_missing = rng.random(size=matrix.shape) < (config.missing_fraction / 12.0)
    matrix[random_missing] = np.nan


def run_benchmark(config: BenchmarkConfig) -> BenchmarkResult:
    inputs = build_inputs(config)
    samples = int(inputs.matrix.shape[1])
    analysis_matrix_memory_bytes = int(inputs.matrix.nbytes)
    covariance_matrix_memory_bytes = int(
        samples * samples * np.dtype(np.float64).itemsize
    )
    approximate_matrix_memory_bytes = int(
        inputs.matrix.nbytes
        + inputs.design.nbytes
        + inputs.contrasts.nbytes
        + covariance_matrix_memory_bytes
    )

    tracemalloc.start()
    workflow_started = time.perf_counter()

    correlation_started = time.perf_counter()
    consensus = estimate_duplicate_correlation_reml_consensus(
        inputs.matrix,
        inputs.design,
        inputs.block_ids,
        feature_ids=inputs.feature_ids,
        design_column_names=inputs.coefficient_names,
    )
    correlation_estimation_seconds = time.perf_counter() - correlation_started
    if not consensus.success or consensus.consensus_correlation is None:
        raise RuntimeError(
            "duplicate-correlation benchmark consensus failed: "
            f"{consensus.failure_reason}"
        )

    gls_started = time.perf_counter()
    gls_fit = fit_duplicate_correlation_gls(
        inputs.matrix,
        inputs.design,
        inputs.block_ids,
        consensus.consensus_correlation,
        feature_ids=inputs.feature_ids,
        coefficient_names=inputs.coefficient_names,
        contrasts=inputs.contrasts,
        contrast_names=inputs.contrast_names,
    )
    gls_seconds = time.perf_counter() - gls_started
    total_workflow_seconds = time.perf_counter() - workflow_started
    _, peak_tracemalloc_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    row_masks = tuple(tuple(row) for row in np.isfinite(inputs.matrix))
    unique_missingness_masks = len(set(row_masks))
    metrics: dict[str, float | int] = {
        "sites": config.sites,
        "samples": samples,
        "blocks": config.blocks,
        "observations_per_block": config.observations_per_block,
        "missing_fraction": config.missing_fraction,
        "correlation_estimation_seconds": correlation_estimation_seconds,
        "gls_seconds": gls_seconds,
        "total_workflow_seconds": total_workflow_seconds,
        "estimated_feature_correlation_count": consensus.estimated_feature_count,
        "consensus_correlation": float(consensus.consensus_correlation),
        "unique_missingness_masks": unique_missingness_masks,
        "factorization_cache_size": int(gls_fit.factorization_cache_size),
        "factorization_cache_hit_count": int(gls_fit.factorization_cache_hit_count),
        "factorization_cache_reuse_count": int(gls_fit.factorization_cache_hit_count),
        "analysis_matrix_memory_bytes": analysis_matrix_memory_bytes,
        "covariance_matrix_memory_bytes": covariance_matrix_memory_bytes,
        "approximate_matrix_memory_bytes": approximate_matrix_memory_bytes,
        "peak_tracemalloc_bytes": int(peak_tracemalloc_bytes),
        "approximate_peak_matrix_memory_bytes": max(
            approximate_matrix_memory_bytes,
            int(peak_tracemalloc_bytes),
        ),
    }
    return BenchmarkResult(config=config, metrics=metrics)


def parse_args(argv: list[str] | None = None) -> BenchmarkConfig:
    defaults = default_config()
    parser = argparse.ArgumentParser(
        description="Measure duplicate-correlation estimator and GLS hot paths."
    )
    parser.add_argument("--sites", type=int, default=defaults.sites)
    parser.add_argument("--blocks", type=int, default=defaults.blocks)
    parser.add_argument(
        "--observations-per-block",
        type=int,
        default=defaults.observations_per_block,
    )
    parser.add_argument(
        "--missing-fraction", type=float, default=defaults.missing_fraction
    )
    parser.add_argument("--seed", type=int, default=defaults.seed)
    parser.add_argument(
        "--block-correlation", type=float, default=defaults.block_correlation
    )
    parser.add_argument("--effect-scale", type=float, default=defaults.effect_scale)
    args = parser.parse_args(argv)
    return BenchmarkConfig(
        sites=args.sites,
        blocks=args.blocks,
        observations_per_block=args.observations_per_block,
        missing_fraction=args.missing_fraction,
        seed=args.seed,
        block_correlation=args.block_correlation,
        effect_scale=args.effect_scale,
    )


def report_payload(result: BenchmarkResult) -> dict[str, Any]:
    return {
        "benchmark": "duplicate_correlation_performance",
        "observation_scope": (
            "bounded local machine observation; not a mandatory release-scale gate"
        ),
        "config": asdict(result.config),
        "metrics": dict(result.metrics),
    }


def main(argv: list[str] | None = None) -> None:
    result = run_benchmark(parse_args(argv))
    payload_json = json.dumps(report_payload(result), sort_keys=True)
    print(f"duplicate_correlation_performance_payload_json={payload_json}")
    for key in DOCUMENTED_MAIN_METRIC_KEYS:
        print(f"{key}={result.metrics[key]}")


if __name__ == "__main__":
    main()
