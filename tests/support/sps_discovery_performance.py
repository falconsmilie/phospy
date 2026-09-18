from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.advanced import (
    SpsDiscoveryConfig,
    SpsDiscoveryRequest,
    SpsDiscoveryResult,
    SpsDiscoveryWorkflow,
    SpsReferenceDataset,
)
from phospy.science.sites.site_keys import (
    ProteinScopedPhosphositeKey,
    encode_site_key,
)

SPS_BENCHMARK_SEED = 20260917


@dataclass(frozen=True, slots=True)
class SpsDiscoveryBenchmarkConfig:
    n_references: int
    n_sites: int
    n_samples: int
    n_conditions: int
    top_n: int = 100

    @property
    def replicates_per_condition(self) -> int:
        return self.n_samples // self.n_conditions


@dataclass(frozen=True, slots=True)
class SpsDiscoveryBenchmarkRun:
    result: SpsDiscoveryResult
    runtime_seconds: float


MODERATE_SPS_BENCHMARK_CONFIG = SpsDiscoveryBenchmarkConfig(
    n_references=3,
    n_sites=10_000,
    n_samples=24,
    n_conditions=3,
)


def build_sps_benchmark_references(
    config: SpsDiscoveryBenchmarkConfig,
) -> tuple[SpsReferenceDataset, ...]:
    if config.n_samples % config.n_conditions != 0:
        raise ValueError("SPS benchmark samples must divide evenly by condition")
    site_keys = _benchmark_site_keys(config.n_sites)
    site_positions = np.arange(config.n_sites, dtype=int)
    references: list[SpsReferenceDataset] = []
    for reference_number in range(config.n_references):
        rng = np.random.default_rng(SPS_BENCHMARK_SEED + reference_number)
        columns, condition_by_sample = _sample_layout(
            reference_number=reference_number,
            config=config,
        )
        stable_component = (
            ((site_positions * 17 + reference_number * 11) % 101) - 50
        ) / 500.0
        condition_effects = np.tile(
            np.linspace(
                -0.35,
                0.35,
                num=config.n_conditions,
                dtype=float,
            ),
            (config.n_sites, 1),
        )
        condition_effects *= (
            1.0 + ((site_positions + reference_number) % 17)[:, None] / 20.0
        )
        condition_effects = np.repeat(
            condition_effects,
            config.replicates_per_condition,
            axis=1,
        )
        replicate_noise = rng.normal(
            loc=0.0,
            scale=0.025,
            size=(config.n_sites, config.n_samples),
        )
        values = np.round(
            stable_component[:, None] + condition_effects + replicate_noise,
            decimals=6,
        )

        # Each reference loses a different, non-overlapping subset of sites from
        # one complete condition. Those sites remain eligible through the other
        # references and exercise partial-reference consensus contribution.
        missing_rows = (site_positions % 307) == (reference_number * 79)
        missing_condition = reference_number % config.n_conditions
        missing_start = missing_condition * config.replicates_per_condition
        missing_stop = missing_start + config.replicates_per_condition
        values[missing_rows, missing_start:missing_stop] = np.nan

        row_order = np.roll(
            site_positions,
            (reference_number * 997) % config.n_sites,
        )
        sample_order = np.roll(
            np.arange(config.n_samples, dtype=int),
            reference_number * config.replicates_per_condition,
        )
        matrix = pd.DataFrame(
            values[row_order][:, sample_order],
            index=site_keys[row_order],
            columns=pd.Index(columns)[sample_order],
            dtype=float,
        )
        references.append(
            SpsReferenceDataset.from_condition_relative_log2(
                dataset_id=f"sps_benchmark_reference_{reference_number + 1}",
                intensities=matrix,
                condition_by_sample=condition_by_sample,
                log2_scale_established_by=(
                    "deterministic SPS benchmark log2 preparation"
                ),
                baseline_centering_established_by=(
                    "deterministic SPS benchmark reference centering"
                ),
                organism="rat",
                baseline_context="benchmark control condition",
                reference_context="synthetic rat SPS performance fixture",
                source_name="deterministic-sps-performance-fixture",
                source_version="1",
                source_uri="https://example.test/sps/performance-fixture",
            )
        )
    return tuple(references)


def run_sps_discovery_benchmark(
    references: tuple[SpsReferenceDataset, ...],
    *,
    config: SpsDiscoveryBenchmarkConfig,
) -> SpsDiscoveryBenchmarkRun:
    request = SpsDiscoveryRequest(
        reference_datasets=references,
        config=SpsDiscoveryConfig(
            top_n=config.top_n,
            minimum_reference_datasets=config.n_references,
            minimum_datasets_per_site=max(2, config.n_references - 1),
            minimum_shared_sites=config.top_n,
        ),
    )
    started = time.perf_counter()
    result = SpsDiscoveryWorkflow().run(request)
    runtime_seconds = time.perf_counter() - started
    return SpsDiscoveryBenchmarkRun(
        result=result,
        runtime_seconds=runtime_seconds,
    )


def _benchmark_site_keys(n_sites: int) -> pd.Index:
    residues = ("S", "T", "Y")
    return pd.Index(
        [
            encode_site_key(
                ProteinScopedPhosphositeKey(
                    organism="human",
                    protein_namespace="uniprot",
                    protein_identifier=f"BMARK{position:06d}",
                    residue=residues[position % len(residues)],
                    position=(position % 2_000) + 1,
                )
            )
            for position in range(n_sites)
        ],
        name="site_key",
    )


def _sample_layout(
    *,
    reference_number: int,
    config: SpsDiscoveryBenchmarkConfig,
) -> tuple[tuple[str, ...], dict[str, str]]:
    columns = tuple(
        f"reference_{reference_number + 1}_condition_{condition + 1}_"
        f"replicate_{replicate + 1}"
        for condition in range(config.n_conditions)
        for replicate in range(config.replicates_per_condition)
    )
    conditions = {
        sample_id: f"condition_{condition + 1}"
        for condition in range(config.n_conditions)
        for sample_id in columns[
            condition * config.replicates_per_condition : (condition + 1)
            * config.replicates_per_condition
        ]
    }
    return columns, conditions
