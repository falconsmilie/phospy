from __future__ import annotations

import time
import tracemalloc
from collections.abc import Callable
from typing import NamedTuple

import numpy as np
import pandas as pd

DEFAULT_PERFORMANCE_SEED = 20260426
_AMINO_ACIDS = tuple("ARNDCEQGHILKMFPSTWYV")


class KnnImputationBenchmarkTier(NamedTuple):
    case_id: str
    site_count: int
    sample_count: int
    missing_target_rows: int
    missing_cells_per_target_row: int
    runtime_seconds_max: float


# Representative performance fixture dimensions (CI-safe, scientifically realistic).
PREPROCESSING_CONTRACT_N_SITES = 5_000
PREPROCESSING_CONTRACT_N_SAMPLES = 12
PREPROCESSING_CONTRACT_MISSING_FRACTION = 0.12
WORKFLOW_SMOKE_CONTRACT_N_SITES = 800
WORKFLOW_SMOKE_CONTRACT_N_SAMPLES = 8
WORKFLOW_SMOKE_CONTRACT_N_CONDITIONS = 2
WORKFLOW_SMOKE_CONTRACT_MISSING_FRACTION = 0.08
WORKFLOW_MEDIUM_CONTRACT_N_SITES = 3_000
WORKFLOW_MEDIUM_CONTRACT_N_SAMPLES = 12
WORKFLOW_MEDIUM_CONTRACT_N_CONDITIONS = 4
WORKFLOW_MEDIUM_CONTRACT_MISSING_FRACTION = 0.18
SIGNALOME_CONTRACT_N_SITES = 2_000
SIGNALOME_CONTRACT_N_KINASES = 100
ADAPTIVE_PREDICTION_CONTRACT_N_SITES = 2_000
ADAPTIVE_PREDICTION_CONTRACT_N_KINASES = 100
ADAPTIVE_PREDICTION_CONTRACT_CANDIDATE_KINASES = 12
ADAPTIVE_PREDICTION_CONTRACT_TOP_K = 48
SSGSEA_ACTIVITY_CONTRACT_N_SITES = 720
SSGSEA_ACTIVITY_CONTRACT_N_KINASES = 32
SSGSEA_ACTIVITY_CONTRACT_N_PROFILES = 6
SSGSEA_ACTIVITY_CONTRACT_SUBSTRATES_PER_KINASE = 24
SSGSEA_ACTIVITY_CONTRACT_PERMUTATIONS = 48
DATASET_VALIDATION_SMALL_N_SITES = 100
DATASET_VALIDATION_SMALL_N_SAMPLES = 6
DATASET_VALIDATION_MEDIUM_N_SITES = 10_000
DATASET_VALIDATION_MEDIUM_N_SAMPLES = 24

# CI benchmark/performance ceilings (loose enough for stability, strict enough
# to catch major accidental regressions).
SIGNALOME_BELOW_THRESHOLD_RUNTIME_SECONDS_MAX = 8.0
SIGNALOME_ABOVE_THRESHOLD_RUNTIME_SECONDS_MAX = 8.0

QUANTILE_RUNTIME_SECONDS_MAX = 12.0
QUANTILE_PEAK_MIB_MAX = 384.0
KNN_IMPUTATION_BENCHMARK_SITE_COUNTS = (10_000, 25_000, 50_000)
KNN_IMPUTATION_BENCHMARK_N_SAMPLES = 12
KNN_IMPUTATION_BENCHMARK_MISSING_TARGET_ROWS = 96
KNN_IMPUTATION_RUNTIME_SECONDS_MAX_BY_SITE_COUNT = {
    10_000: 5.0,
    25_000: 8.0,
    50_000: 12.0,
}
KNN_IMPUTATION_PEAK_MIB_MAX = 384.0
KNN_IMPUTATION_BENCHMARK_TIERS = (
    KnnImputationBenchmarkTier(
        case_id="sparse_10k_12samples_96targets",
        site_count=10_000,
        sample_count=12,
        missing_target_rows=96,
        missing_cells_per_target_row=1,
        runtime_seconds_max=5.0,
    ),
    KnnImputationBenchmarkTier(
        case_id="sparse_25k_12samples_96targets",
        site_count=25_000,
        sample_count=12,
        missing_target_rows=96,
        missing_cells_per_target_row=1,
        runtime_seconds_max=8.0,
    ),
    KnnImputationBenchmarkTier(
        case_id="sparse_50k_12samples_96targets",
        site_count=50_000,
        sample_count=12,
        missing_target_rows=96,
        missing_cells_per_target_row=1,
        runtime_seconds_max=12.0,
    ),
    KnnImputationBenchmarkTier(
        case_id="moderate_10k_24samples_256targets",
        site_count=10_000,
        sample_count=24,
        missing_target_rows=256,
        missing_cells_per_target_row=2,
        runtime_seconds_max=8.0,
    ),
    KnnImputationBenchmarkTier(
        case_id="moderate_25k_24samples_512targets",
        site_count=25_000,
        sample_count=24,
        missing_target_rows=512,
        missing_cells_per_target_row=2,
        runtime_seconds_max=15.0,
    ),
    KnnImputationBenchmarkTier(
        case_id="moderate_50k_24samples_768targets",
        site_count=50_000,
        sample_count=24,
        missing_target_rows=768,
        missing_cells_per_target_row=2,
        runtime_seconds_max=30.0,
    ),
)

DATASET_BUILD_SMOKE_RUNTIME_SECONDS_MAX = 40.0
DATASET_BUILD_SMOKE_PEAK_MIB_MAX = 320.0
DATASET_BUILD_MEDIUM_RUNTIME_SECONDS_MAX = 210.0
DATASET_BUILD_MEDIUM_PEAK_MIB_MAX = 900.0

DIFFERENTIAL_WORKFLOW_SMOKE_RUNTIME_SECONDS_MAX = 15.0
DIFFERENTIAL_WORKFLOW_SMOKE_PEAK_MIB_MAX = 320.0
DIFFERENTIAL_WORKFLOW_MEDIUM_RUNTIME_SECONDS_MAX = 140.0
DIFFERENTIAL_WORKFLOW_MEDIUM_PEAK_MIB_MAX = 640.0
EMPIRICAL_BAYES_TREND_SMALL_N_FEATURES = 1_000
EMPIRICAL_BAYES_TREND_MEDIUM_N_FEATURES = 10_000
EMPIRICAL_BAYES_TREND_LARGE_N_FEATURES = 50_000
EMPIRICAL_BAYES_TREND_SMALL_RUNTIME_SECONDS_MAX = 3.0
EMPIRICAL_BAYES_TREND_MEDIUM_RUNTIME_SECONDS_MAX = 8.0
EMPIRICAL_BAYES_TREND_LARGE_RUNTIME_SECONDS_MAX = 30.0
EMPIRICAL_BAYES_TREND_SMALL_PEAK_MIB_MAX = 64.0
EMPIRICAL_BAYES_TREND_MEDIUM_PEAK_MIB_MAX = 128.0
EMPIRICAL_BAYES_TREND_LARGE_PEAK_MIB_MAX = 256.0

MOTIF_RUNTIME_SECONDS_MAX = 18.0
MOTIF_PEAK_MIB_MAX = 640.0

# Includes workflow provenance hashing; keep broad enough for slower Windows
# runners while still catching order-of-magnitude reference-filtering regressions.
KINASE_FILTERED_REFERENCE_RUNTIME_SECONDS_MAX = 15.0
KINASE_FILTERED_REFERENCE_PEAK_MIB_MAX = 300.0

SSGSEA_ACTIVITY_RUNTIME_SECONDS_MAX = 12.0
SSGSEA_ACTIVITY_PEAK_MIB_MAX = 192.0

DIAGNOSTIC_RUNTIME_RATIO_MULTIPLIER = 5.0
DIAGNOSTIC_RUNTIME_ABSOLUTE_SECONDS = 1.0

ADAPTIVE_PREDICTION_RUNTIME_SECONDS_MAX = 35.0
ADAPTIVE_PREDICTION_PEAK_MIB_MAX = 900.0

SIGNALOME_BACKEND_SMALL_RUNTIME_SECONDS_MAX = 10.0
SIGNALOME_BACKEND_MEDIUM_RUNTIME_SECONDS_MAX = 18.0
SIGNALOME_BACKEND_SMALL_PEAK_MIB_MAX = 256.0
SIGNALOME_BACKEND_MEDIUM_PEAK_MIB_MAX = 512.0
SIGNALOME_NEAR_THRESHOLD_RUNTIME_SECONDS_MAX = 12.0
SIGNALOME_FULL_GUARD_RUNTIME_SECONDS_MAX = 3.0
SIGNALOME_WORKFLOW_RUNTIME_SECONDS_MAX = 20.0
SIGNALOME_WORKFLOW_PEAK_MIB_MAX = 700.0
SIGNALOME_WORKFLOW_PRECONDITIONED_RUNTIME_SECONDS_MAX = 20.0
SIGNALOME_WORKFLOW_PRECONDITIONED_PEAK_MIB_MAX = 700.0

SIGNALOME_CLUSTER_TREE_BENCHMARK_N_SITES = 500
SIGNALOME_CLUSTER_TREE_BENCHMARK_N_KINASES = 8
SIGNALOME_CLUSTER_TREE_RUNTIME_SECONDS_MAX = 25.0
SIGNALOME_CLUSTER_TREE_PEAK_MIB_MAX = 256.0

BUNDLE_PUBLISH_RUNTIME_SECONDS_MAX = 25.0
BUNDLE_PUBLISH_PEAK_MIB_MAX = 700.0

PROVENANCE_HASHING_RUNTIME_SECONDS_MAX = 18.0
PROVENANCE_HASHING_PEAK_MIB_MAX = 768.0

DATASET_VALIDATION_MEDIUM_CONSTRUCTION_RUNTIME_SECONDS_MAX = 420.0
DATASET_VALIDATION_BOUNDED_SITE_METADATA_RUNTIME_SECONDS_MAX = 120.0
DATASET_VALIDATION_BOUNDED_ALIGNMENT_RUNTIME_SECONDS_MAX = 45.0
DATASET_VALIDATION_BOUNDED_OBSERVATION_MASK_RUNTIME_SECONDS_MAX = 20.0
IMPORTER_MEDIUM_NORMALISATION_RUNTIME_SECONDS_MAX = 180.0
WORKFLOW_VALIDATION_BOUNDED_SEQUENCE_RUNTIME_SECONDS_MAX = 45.0

_PHOSPHO_RESIDUES = ("S", "T", "Y")


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
        values.flat[missing_positions] = np.nan
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


def deterministic_analysis_ready_site_keys(
    count: int,
    *,
    start: int = 1,
    gene_prefix: str = "PERFGENE",
    organism: str = "rat",
    protein_namespace: str = "protein_id",
) -> pd.Index:
    """Return deterministic encoded site_key labels for validation scale gates."""

    labels: list[str] = []
    for offset, position in enumerate(range(start, start + int(count))):
        residue = _PHOSPHO_RESIDUES[offset % len(_PHOSPHO_RESIDUES)]
        protein_identifier = f"{gene_prefix}{position:05d}"
        labels.append(
            "phospy:v1|"
            f"organism={organism}|"
            f"protein_namespace={protein_namespace}|"
            f"protein_identifier={protein_identifier}|"
            f"residue={residue}|"
            f"position={position}"
        )
    return pd.Index(labels, name="site_key")


def deterministic_analysis_ready_site_metadata(
    site_keys: pd.Index,
    *,
    start: int = 1,
    gene_prefix: str = "PERFGENE",
    organism: str = "rat",
    protein_namespace: str = "protein_id",
    sequence_width: int = 31,
    localisation_confidence: float = 0.95,
) -> pd.DataFrame:
    """Build strict analysis-ready site metadata aligned to encoded site keys."""

    if sequence_width < 3 or sequence_width % 2 == 0:
        raise ValueError("sequence_width must be an odd integer >= 3")
    count = int(site_keys.size)
    positions = list(range(start, start + count))
    residues = [
        _PHOSPHO_RESIDUES[index % len(_PHOSPHO_RESIDUES)] for index in range(count)
    ]
    protein_identifiers = [f"{gene_prefix}{position:05d}" for position in positions]
    sites = [
        f"{residue}{position}"
        for residue, position in zip(residues, positions, strict=True)
    ]
    center = sequence_width // 2
    sequences: list[str] = []
    for row_index, residue in enumerate(residues):
        letters = [
            _AMINO_ACIDS[(row_index * 11 + column_index * 7) % len(_AMINO_ACIDS)]
            for column_index in range(sequence_width)
        ]
        letters[center] = residue
        sequences.append("".join(letters))
    return pd.DataFrame(
        {
            "site_key": site_keys.astype(str).tolist(),
            "display_id": [
                f"{protein_identifier};{site};"
                for protein_identifier, site in zip(
                    protein_identifiers,
                    sites,
                    strict=True,
                )
            ],
            "organism": [organism] * count,
            "protein_namespace": [protein_namespace] * count,
            "protein_identifier": protein_identifiers,
            "gene_symbol": protein_identifiers,
            "site": sites,
            "protein_id": protein_identifiers,
            "site_sequence": sequences,
            "localisation_confidence": [float(localisation_confidence)] * count,
        },
        index=site_keys.copy(),
    )


def deterministic_analysis_ready_dataset_tables(
    *,
    n_sites: int,
    n_samples: int,
    seed: int = DEFAULT_PERFORMANCE_SEED,
    start: int = 1,
    gene_prefix: str = "PERFGENE",
    sequence_width: int = 31,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return aligned phospho/site-metadata frames for constructor gates."""

    site_keys = deterministic_analysis_ready_site_keys(
        n_sites,
        start=start,
        gene_prefix=gene_prefix,
    )
    phospho = deterministic_matrix(
        n_sites=n_sites,
        n_samples=n_samples,
        seed=seed,
        site_ids=site_keys,
    )
    phospho = phospho.clip(lower=0.0)
    site_metadata = deterministic_analysis_ready_site_metadata(
        site_keys,
        start=start,
        gene_prefix=gene_prefix,
        sequence_width=sequence_width,
    )
    return phospho, site_metadata


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
            "localisation_confidence": [0.95] * len(genes),
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


def deterministic_maxquant_source_table(
    *,
    n_sites: int,
    n_samples: int,
    seed: int = DEFAULT_PERFORMANCE_SEED,
    start: int = 1,
    gene_prefix: str = "MQGENE",
) -> pd.DataFrame:
    """Build a deterministic MaxQuant-like phosphosite table for importer gates."""

    positions = list(range(start, start + int(n_sites)))
    residues = [
        _PHOSPHO_RESIDUES[index % len(_PHOSPHO_RESIDUES)]
        for index in range(int(n_sites))
    ]
    proteins = [f"{gene_prefix}{position:05d}" for position in positions]
    sites = [
        f"{residue}{position}"
        for residue, position in zip(residues, positions, strict=True)
    ]
    sequence_windows: list[str] = []
    for row_index, residue in enumerate(residues):
        letters = [
            _AMINO_ACIDS[(row_index * 5 + column_index * 3) % len(_AMINO_ACIDS)]
            for column_index in range(15)
        ]
        letters[7] = residue
        sequence_windows.append("".join(letters))
    source = pd.DataFrame(
        {
            "Leading proteins": proteins,
            "Gene names": proteins,
            "Modified site": sites,
            "Localization prob": [0.95] * int(n_sites),
            "Sequence": [f"PEPTIDE{index:05d}" for index in range(int(n_sites))],
            "Modified sequence": [
                f"PEPTIDE(ph){index:05d}" for index in range(int(n_sites))
            ],
            "Sequence window": sequence_windows,
            "Potential contaminant": [""] * int(n_sites),
            "Reverse": [""] * int(n_sites),
        }
    )
    values = deterministic_matrix(
        n_sites=n_sites,
        n_samples=n_samples,
        seed=seed,
        site_ids=pd.Index([f"row_{index}" for index in range(int(n_sites))]),
        sample_columns=pd.Index(
            [f"sample_{index + 1:02d}" for index in range(int(n_samples))]
        ),
    )
    for sample_id in values.columns.astype(str).tolist():
        source[f"Intensity {sample_id}"] = values.loc[:, sample_id].astype(str).tolist()
    return source


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


def measure_wall_clock(
    func: Callable[[], object],
    *,
    warmup: bool = True,
) -> tuple[object, float]:
    """Measure ordinary wall-clock runtime without allocation tracing."""

    if warmup:
        func()
    started = time.perf_counter()
    result = func()
    elapsed = time.perf_counter() - started
    return result, float(elapsed)


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
    "ADAPTIVE_PREDICTION_CONTRACT_CANDIDATE_KINASES",
    "ADAPTIVE_PREDICTION_CONTRACT_N_KINASES",
    "ADAPTIVE_PREDICTION_CONTRACT_N_SITES",
    "ADAPTIVE_PREDICTION_CONTRACT_TOP_K",
    "ADAPTIVE_PREDICTION_PEAK_MIB_MAX",
    "ADAPTIVE_PREDICTION_RUNTIME_SECONDS_MAX",
    "BUNDLE_PUBLISH_PEAK_MIB_MAX",
    "BUNDLE_PUBLISH_RUNTIME_SECONDS_MAX",
    "DATASET_BUILD_MEDIUM_PEAK_MIB_MAX",
    "DATASET_BUILD_MEDIUM_RUNTIME_SECONDS_MAX",
    "DATASET_BUILD_SMOKE_PEAK_MIB_MAX",
    "DATASET_BUILD_SMOKE_RUNTIME_SECONDS_MAX",
    "DATASET_VALIDATION_BOUNDED_ALIGNMENT_RUNTIME_SECONDS_MAX",
    "DATASET_VALIDATION_BOUNDED_OBSERVATION_MASK_RUNTIME_SECONDS_MAX",
    "DATASET_VALIDATION_BOUNDED_SITE_METADATA_RUNTIME_SECONDS_MAX",
    "DATASET_VALIDATION_MEDIUM_CONSTRUCTION_RUNTIME_SECONDS_MAX",
    "DATASET_VALIDATION_MEDIUM_N_SAMPLES",
    "DATASET_VALIDATION_MEDIUM_N_SITES",
    "DATASET_VALIDATION_SMALL_N_SAMPLES",
    "DATASET_VALIDATION_SMALL_N_SITES",
    "DEFAULT_PERFORMANCE_SEED",
    "DIFFERENTIAL_WORKFLOW_MEDIUM_PEAK_MIB_MAX",
    "DIFFERENTIAL_WORKFLOW_MEDIUM_RUNTIME_SECONDS_MAX",
    "DIFFERENTIAL_WORKFLOW_SMOKE_PEAK_MIB_MAX",
    "DIFFERENTIAL_WORKFLOW_SMOKE_RUNTIME_SECONDS_MAX",
    "DIAGNOSTIC_RUNTIME_ABSOLUTE_SECONDS",
    "DIAGNOSTIC_RUNTIME_RATIO_MULTIPLIER",
    "EMPIRICAL_BAYES_TREND_LARGE_N_FEATURES",
    "EMPIRICAL_BAYES_TREND_LARGE_PEAK_MIB_MAX",
    "EMPIRICAL_BAYES_TREND_LARGE_RUNTIME_SECONDS_MAX",
    "EMPIRICAL_BAYES_TREND_MEDIUM_N_FEATURES",
    "EMPIRICAL_BAYES_TREND_MEDIUM_PEAK_MIB_MAX",
    "EMPIRICAL_BAYES_TREND_MEDIUM_RUNTIME_SECONDS_MAX",
    "EMPIRICAL_BAYES_TREND_SMALL_N_FEATURES",
    "EMPIRICAL_BAYES_TREND_SMALL_PEAK_MIB_MAX",
    "EMPIRICAL_BAYES_TREND_SMALL_RUNTIME_SECONDS_MAX",
    "IMPORTER_MEDIUM_NORMALISATION_RUNTIME_SECONDS_MAX",
    "KINASE_FILTERED_REFERENCE_PEAK_MIB_MAX",
    "KINASE_FILTERED_REFERENCE_RUNTIME_SECONDS_MAX",
    "KNN_IMPUTATION_BENCHMARK_TIERS",
    "KNN_IMPUTATION_BENCHMARK_MISSING_TARGET_ROWS",
    "KNN_IMPUTATION_BENCHMARK_N_SAMPLES",
    "KNN_IMPUTATION_BENCHMARK_SITE_COUNTS",
    "KNN_IMPUTATION_PEAK_MIB_MAX",
    "KNN_IMPUTATION_RUNTIME_SECONDS_MAX_BY_SITE_COUNT",
    "KnnImputationBenchmarkTier",
    "MOTIF_PEAK_MIB_MAX",
    "MOTIF_RUNTIME_SECONDS_MAX",
    "PREPROCESSING_CONTRACT_MISSING_FRACTION",
    "PREPROCESSING_CONTRACT_N_SAMPLES",
    "PREPROCESSING_CONTRACT_N_SITES",
    "PROVENANCE_HASHING_PEAK_MIB_MAX",
    "PROVENANCE_HASHING_RUNTIME_SECONDS_MAX",
    "QUANTILE_PEAK_MIB_MAX",
    "QUANTILE_RUNTIME_SECONDS_MAX",
    "SIGNALOME_ABOVE_THRESHOLD_RUNTIME_SECONDS_MAX",
    "SIGNALOME_BACKEND_MEDIUM_PEAK_MIB_MAX",
    "SIGNALOME_BACKEND_MEDIUM_RUNTIME_SECONDS_MAX",
    "SIGNALOME_BACKEND_SMALL_PEAK_MIB_MAX",
    "SIGNALOME_BACKEND_SMALL_RUNTIME_SECONDS_MAX",
    "SIGNALOME_BELOW_THRESHOLD_RUNTIME_SECONDS_MAX",
    "SIGNALOME_CLUSTER_TREE_BENCHMARK_N_KINASES",
    "SIGNALOME_CLUSTER_TREE_BENCHMARK_N_SITES",
    "SIGNALOME_CLUSTER_TREE_PEAK_MIB_MAX",
    "SIGNALOME_CLUSTER_TREE_RUNTIME_SECONDS_MAX",
    "SIGNALOME_CONTRACT_N_KINASES",
    "SIGNALOME_CONTRACT_N_SITES",
    "SIGNALOME_FULL_GUARD_RUNTIME_SECONDS_MAX",
    "SIGNALOME_NEAR_THRESHOLD_RUNTIME_SECONDS_MAX",
    "SIGNALOME_WORKFLOW_PEAK_MIB_MAX",
    "SIGNALOME_WORKFLOW_PRECONDITIONED_PEAK_MIB_MAX",
    "SIGNALOME_WORKFLOW_PRECONDITIONED_RUNTIME_SECONDS_MAX",
    "SIGNALOME_WORKFLOW_RUNTIME_SECONDS_MAX",
    "SSGSEA_ACTIVITY_CONTRACT_N_KINASES",
    "SSGSEA_ACTIVITY_CONTRACT_N_PROFILES",
    "SSGSEA_ACTIVITY_CONTRACT_N_SITES",
    "SSGSEA_ACTIVITY_CONTRACT_PERMUTATIONS",
    "SSGSEA_ACTIVITY_CONTRACT_SUBSTRATES_PER_KINASE",
    "SSGSEA_ACTIVITY_PEAK_MIB_MAX",
    "SSGSEA_ACTIVITY_RUNTIME_SECONDS_MAX",
    "WORKFLOW_VALIDATION_BOUNDED_SEQUENCE_RUNTIME_SECONDS_MAX",
    "deterministic_analysis_ready_dataset_tables",
    "deterministic_analysis_ready_site_keys",
    "deterministic_analysis_ready_site_metadata",
    "deterministic_kinase_substrate_map",
    "deterministic_maxquant_source_table",
    "deterministic_matrix",
    "deterministic_sample_columns",
    "deterministic_site_ids",
    "deterministic_site_metadata",
    "deterministic_site_sequence_frame",
    "deterministic_site_sequence_series",
    "median_runtime_and_peak_mib",
    "median_runtime_seconds",
    "measure_wall_clock",
    "measure_runtime_and_peak_mib",
    "with_missing_fraction",
    "WORKFLOW_MEDIUM_CONTRACT_MISSING_FRACTION",
    "WORKFLOW_MEDIUM_CONTRACT_N_CONDITIONS",
    "WORKFLOW_MEDIUM_CONTRACT_N_SAMPLES",
    "WORKFLOW_MEDIUM_CONTRACT_N_SITES",
    "WORKFLOW_SMOKE_CONTRACT_MISSING_FRACTION",
    "WORKFLOW_SMOKE_CONTRACT_N_CONDITIONS",
    "WORKFLOW_SMOKE_CONTRACT_N_SAMPLES",
    "WORKFLOW_SMOKE_CONTRACT_N_SITES",
]
