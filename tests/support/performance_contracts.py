from __future__ import annotations

import time
import tracemalloc
from collections.abc import Callable

import numpy as np
import pandas as pd

DEFAULT_PERFORMANCE_SEED = 20260426
_AMINO_ACIDS = tuple("ARNDCEQGHILKMFPSTWYV")


def deterministic_site_ids(
    count: int,
    *,
    start: int = 1,
    gene_prefix: str = "GENE",
) -> pd.Index:
    return pd.Index(
        [
            f"{gene_prefix}{position:05d};S{((position - 1) % 401) + 1};"
            for position in range(start, start + count)
        ],
        name="site_id",
    )


def deterministic_sample_columns(
    count: int,
    *,
    prefix: str = "sample",
) -> pd.Index:
    return pd.Index([f"{prefix}_{index + 1:02d}" for index in range(count)])


def deterministic_matrix(
    *,
    n_sites: int,
    n_samples: int,
    seed: int = DEFAULT_PERFORMANCE_SEED,
    site_ids: pd.Index | None = None,
    sample_columns: pd.Index | None = None,
) -> pd.DataFrame:
    resolved_site_ids = (
        deterministic_site_ids(n_sites) if site_ids is None else site_ids.copy()
    )
    resolved_sample_columns = (
        deterministic_sample_columns(n_samples)
        if sample_columns is None
        else sample_columns.copy()
    )
    rng = np.random.default_rng(seed)
    values = rng.normal(
        loc=10.0,
        scale=2.5,
        size=(int(n_sites), int(n_samples)),
    ).astype(float)
    values = np.round(values, decimals=6)
    return pd.DataFrame(
        values,
        index=resolved_site_ids,
        columns=resolved_sample_columns,
        dtype=float,
    )


def with_missing_fraction(
    matrix: pd.DataFrame,
    *,
    missing_fraction: float,
    seed: int = DEFAULT_PERFORMANCE_SEED,
) -> pd.DataFrame:
    if missing_fraction <= 0.0:
        return matrix.copy(deep=True)
    if missing_fraction >= 1.0:
        return pd.DataFrame(
            np.nan,
            index=matrix.index.copy(),
            columns=matrix.columns.copy(),
            dtype=float,
        )

    values = matrix.to_numpy(dtype=float, copy=True)
    rng = np.random.default_rng(seed)
    cell_count = int(values.size)
    missing_count = int(round(cell_count * float(missing_fraction)))
    missing_count = min(max(missing_count, 0), cell_count)
    if missing_count > 0:
        missing_positions = rng.choice(cell_count, size=missing_count, replace=False)
        flat = values.reshape(-1)
        flat[missing_positions] = np.nan
    return pd.DataFrame(
        values,
        index=matrix.index.copy(),
        columns=matrix.columns.copy(),
        dtype=float,
    )


def deterministic_site_sequence_series(
    site_ids: pd.Index,
    *,
    window_width: int = 15,
) -> pd.Series:
    if window_width < 1:
        raise ValueError("window_width must be >= 1")
    center = window_width // 2
    sequences: list[str] = []
    for site_index, _site_id in enumerate(site_ids.astype(str).tolist()):
        letters = [
            _AMINO_ACIDS[(site_index * 11 + position * 7) % len(_AMINO_ACIDS)]
            for position in range(window_width)
        ]
        letters[center] = "S"
        sequences.append("".join(letters))
    return pd.Series(sequences, index=site_ids.copy(), dtype=str, name="site_sequence")


def deterministic_site_metadata(
    site_ids: pd.Index,
    *,
    include_protein_id: bool = True,
    sequence_width: int = 15,
) -> pd.DataFrame:
    site_text = site_ids.astype(str)
    genes = [token.split(";")[0] for token in site_text]
    sites = [token.split(";")[1] for token in site_text]
    metadata = pd.DataFrame(
        {
            "gene_symbol": genes,
            "site": sites,
            "site_sequence": deterministic_site_sequence_series(
                site_ids,
                window_width=sequence_width,
            ).tolist(),
        },
        index=site_ids.copy(),
    )
    if include_protein_id:
        metadata.loc[:, "protein_id"] = genes
    return metadata


def deterministic_site_sequence_frame(
    site_ids: pd.Index,
    *,
    sequence_width: int = 15,
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_sequence": deterministic_site_sequence_series(
                site_ids, window_width=sequence_width
            )
        },
        index=site_ids.copy(),
    )


def deterministic_kinase_substrate_map(
    *,
    dataset_site_ids: pd.Index,
    eligible_kinase_count: int,
    substrates_per_kinase: int,
    offlane_kinase_count: int = 0,
    offlane_sites_per_kinase: int | None = None,
    offlane_site_ids: pd.Index | None = None,
) -> pd.DataFrame:
    if eligible_kinase_count < 0:
        raise ValueError("eligible_kinase_count must be >= 0")
    if substrates_per_kinase < 1:
        raise ValueError("substrates_per_kinase must be >= 1")
    if offlane_kinase_count < 0:
        raise ValueError("offlane_kinase_count must be >= 0")

    dataset_sites = np.asarray(dataset_site_ids.astype(str), dtype=object)
    if dataset_sites.size == 0 and eligible_kinase_count > 0:
        raise ValueError("dataset_site_ids must be non-empty when eligible kinases > 0")

    rows: list[dict[str, str]] = []
    for kinase_index in range(int(eligible_kinase_count)):
        kinase = f"KINASE_{kinase_index + 1:03d}"
        start = kinase_index * substrates_per_kinase
        selected_sites = _cyclic_take(
            dataset_sites, start=start, count=substrates_per_kinase
        )
        rows.extend(
            {"kinase": kinase, "substrate_site": str(site_id)}
            for site_id in selected_sites.tolist()
        )

    resolved_offlane_sites_per_kinase = (
        substrates_per_kinase
        if offlane_sites_per_kinase is None
        else int(offlane_sites_per_kinase)
    )
    if resolved_offlane_sites_per_kinase < 1:
        raise ValueError("offlane_sites_per_kinase must be >= 1")
    if offlane_kinase_count > 0:
        resolved_offlane_site_ids = (
            deterministic_site_ids(
                int(offlane_kinase_count) * resolved_offlane_sites_per_kinase + 32,
                start=500_000,
                gene_prefix="OFFSITE",
            )
            if offlane_site_ids is None
            else offlane_site_ids.copy()
        )
        offlane_sites = np.asarray(resolved_offlane_site_ids.astype(str), dtype=object)
        for kinase_index in range(int(offlane_kinase_count)):
            kinase = f"OFFLANE_KINASE_{kinase_index + 1:03d}"
            start = kinase_index * resolved_offlane_sites_per_kinase
            selected_sites = _cyclic_take(
                offlane_sites,
                start=start,
                count=resolved_offlane_sites_per_kinase,
            )
            rows.extend(
                {"kinase": kinase, "substrate_site": str(site_id)}
                for site_id in selected_sites.tolist()
            )

    return pd.DataFrame.from_records(rows, columns=["kinase", "substrate_site"]).astype(
        {"kinase": str, "substrate_site": str}
    )


def median_runtime_seconds(
    func: Callable[[], object],
    *,
    repeats: int = 3,
    warmup: bool = True,
) -> float:
    if repeats < 1:
        raise ValueError("repeats must be >= 1")
    if warmup:
        func()
    durations: list[float] = []
    for _ in range(repeats):
        started = time.perf_counter()
        func()
        durations.append(time.perf_counter() - started)
    return float(np.median(np.asarray(durations, dtype=float)))


def measure_runtime_and_peak_mib(
    func: Callable[[], object],
    *,
    warmup: bool = True,
) -> tuple[object, float, float]:
    if warmup:
        func()
    tracemalloc.start()
    started = time.perf_counter()
    result = func()
    elapsed = time.perf_counter() - started
    _current_bytes, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    peak_mib = float(peak_bytes) / (1024.0 * 1024.0)
    return result, float(elapsed), peak_mib


def median_runtime_and_peak_mib(
    func: Callable[[], object],
    *,
    repeats: int = 3,
    warmup: bool = True,
) -> tuple[object, float, float]:
    if repeats < 1:
        raise ValueError("repeats must be >= 1")
    if warmup:
        func()
    runtimes: list[float] = []
    peaks: list[float] = []
    result: object = None
    for _ in range(repeats):
        tracemalloc.start()
        started = time.perf_counter()
        result = func()
        runtimes.append(time.perf_counter() - started)
        _current_bytes, peak_bytes = tracemalloc.get_traced_memory()
        peaks.append(float(peak_bytes) / (1024.0 * 1024.0))
        tracemalloc.stop()
    return (
        result,
        float(np.median(np.asarray(runtimes, dtype=float))),
        float(np.median(np.asarray(peaks, dtype=float))),
    )


def _cyclic_take(values: np.ndarray, *, start: int, count: int) -> np.ndarray:
    if values.size == 0:
        return np.asarray([], dtype=object)
    offsets = np.arange(count, dtype=int)
    positions = ((int(start) + offsets) % int(values.size)).astype(int, copy=False)
    return values[positions]


__all__ = [
    "DEFAULT_PERFORMANCE_SEED",
    "deterministic_kinase_substrate_map",
    "deterministic_matrix",
    "deterministic_sample_columns",
    "deterministic_site_ids",
    "deterministic_site_metadata",
    "deterministic_site_sequence_frame",
    "deterministic_site_sequence_series",
    "median_runtime_and_peak_mib",
    "median_runtime_seconds",
    "measure_runtime_and_peak_mib",
    "with_missing_fraction",
]
