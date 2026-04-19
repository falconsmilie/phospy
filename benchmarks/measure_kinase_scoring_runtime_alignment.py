#!/usr/bin/env python3
from __future__ import annotations

import sys
import time
import tracemalloc
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def _build_offlane_stress_references(
    *,
    dataset_sites: pd.Index,
    reference_map: pd.DataFrame,
    reference_sequences: pd.DataFrame,
    extra_kinases: int = 160,
    sites_per_kinase: int = 8,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    dataset_site_set = set(dataset_sites.astype(str).tolist())
    offlane_site_pool = [
        str(site_id)
        for site_id in reference_sequences.index.astype(str).tolist()
        if str(site_id) not in dataset_site_set
    ]
    if not offlane_site_pool:
        raise RuntimeError("off-lane benchmark requires reference-only sequence sites")

    expanded_map = reference_map.copy(deep=True)
    pool_size = len(offlane_site_pool)
    for kinase_index in range(extra_kinases):
        kinase = f"OFFLANE_KINASE_{kinase_index + 1:03d}"
        offsets = np.arange(sites_per_kinase, dtype=int)
        start = (kinase_index * sites_per_kinase) % pool_size
        selected_positions = ((start + offsets) % pool_size).astype(int, copy=False)
        selected_sites = [
            offlane_site_pool[position] for position in selected_positions
        ]
        expanded_map = pd.concat(
            [
                expanded_map,
                pd.DataFrame(
                    {
                        "kinase": [kinase for _ in selected_sites],
                        "substrate_site": selected_sites,
                    }
                ),
            ],
            axis=0,
            ignore_index=True,
        )
    return expanded_map, reference_sequences


def _run_once(
    *,
    include_diagnostic_scoring_tables: bool,
    emulate_unfiltered_motif_lane: bool,
) -> tuple[float, float]:
    import phospy.workflows.kinase.executor as kinase_executor
    from phospy import (
        KinasePredictionConfig,
        KinaseScoringConfig,
        KinaseWorkflow,
        KinaseWorkflowRequest,
        ReferenceBundle,
    )
    from phospy.references.resolution import BundledReferenceProvider
    from tests.support.rewrite_fixture_data import build_rat_l6_dataset

    dataset = build_rat_l6_dataset(n_sites=260)
    bundled_references = BundledReferenceProvider().run(dataset.organism)
    expanded_map, expanded_sequences = _build_offlane_stress_references(
        dataset_sites=dataset.phospho.index,
        reference_map=bundled_references.kinase_substrate_map,
        reference_sequences=bundled_references.site_sequences,
    )
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferenceBundle(
            organism=bundled_references.organism,
            kinase_substrate_map=expanded_map,
            site_sequences=expanded_sequences,
        ),
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            include_diagnostic_scoring_tables=include_diagnostic_scoring_tables,
        ),
        prediction_config=KinasePredictionConfig(top_k=6, ensemble_size=12),
        activity_config=None,
    )

    original_build_motif_library = kinase_executor.build_motif_library
    if emulate_unfiltered_motif_lane:
        full_map = request.references.kinase_substrate_map

        def _force_unfiltered_map(*, kinase_substrate_map, site_sequences, flank_size):
            _ = kinase_substrate_map
            return original_build_motif_library(
                kinase_substrate_map=full_map,
                site_sequences=site_sequences,
                flank_size=flank_size,
            )

        kinase_executor.build_motif_library = _force_unfiltered_map
    try:
        tracemalloc.start()
        start = time.perf_counter()
        KinaseWorkflow().run(request)
        elapsed_seconds = time.perf_counter() - start
        _current_bytes, peak_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()
    finally:
        kinase_executor.build_motif_library = original_build_motif_library
    return elapsed_seconds, peak_bytes / (1024.0 * 1024.0)


def _average_measurements(
    *,
    repeats: int,
    include_diagnostic_scoring_tables: bool,
    emulate_unfiltered_motif_lane: bool,
) -> tuple[float, float]:
    runtimes: list[float] = []
    peaks: list[float] = []
    for _ in range(repeats):
        runtime, peak_mib = _run_once(
            include_diagnostic_scoring_tables=include_diagnostic_scoring_tables,
            emulate_unfiltered_motif_lane=emulate_unfiltered_motif_lane,
        )
        runtimes.append(runtime)
        peaks.append(peak_mib)
    return float(sum(runtimes) / len(runtimes)), float(sum(peaks) / len(peaks))


def main() -> None:
    repeats = 5
    default_runtime, default_peak = _average_measurements(
        repeats=repeats,
        include_diagnostic_scoring_tables=False,
        emulate_unfiltered_motif_lane=False,
    )
    default_unfiltered_runtime, default_unfiltered_peak = _average_measurements(
        repeats=repeats,
        include_diagnostic_scoring_tables=False,
        emulate_unfiltered_motif_lane=True,
    )
    diagnostic_runtime, diagnostic_peak = _average_measurements(
        repeats=repeats,
        include_diagnostic_scoring_tables=True,
        emulate_unfiltered_motif_lane=False,
    )

    filtered_vs_unfiltered_speedup = (
        default_unfiltered_runtime / default_runtime
        if default_runtime > 0.0
        else float("inf")
    )
    default_vs_diagnostic_speedup = (
        diagnostic_runtime / default_runtime if default_runtime > 0.0 else float("inf")
    )
    filtered_peak_reduction_mib = default_unfiltered_peak - default_peak
    diagnostic_peak_reduction_mib = diagnostic_peak - default_peak

    print(f"repeats={repeats}")
    print(f"default_runtime_seconds={default_runtime:.6f}")
    print(f"default_peak_mib={default_peak:.3f}")
    print(f"default_unfiltered_runtime_seconds={default_unfiltered_runtime:.6f}")
    print(f"default_unfiltered_peak_mib={default_unfiltered_peak:.3f}")
    print(f"diagnostic_runtime_seconds={diagnostic_runtime:.6f}")
    print(f"diagnostic_peak_mib={diagnostic_peak:.3f}")
    print(
        "default_runtime_speedup_vs_unfiltered_motif="
        f"{filtered_vs_unfiltered_speedup:.3f}x"
    )
    print(f"default_peak_mib_reduction_vs_unfiltered={filtered_peak_reduction_mib:.3f}")
    print(
        "default_runtime_speedup_vs_diagnostic_tables="
        f"{default_vs_diagnostic_speedup:.3f}x"
    )
    print(
        "default_peak_mib_reduction_vs_diagnostic_tables="
        f"{diagnostic_peak_reduction_mib:.3f}"
    )


if __name__ == "__main__":
    main()
