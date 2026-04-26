"""Motif scoring kernels used by kinase workflow scoring."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

AMINO_ACIDS: tuple[str, ...] = (
    "A",
    "R",
    "N",
    "D",
    "C",
    "E",
    "Q",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
)

DEFAULT_MOTIF_FLANK_SIZE = 7
_ASCII_LOOKUP_SIZE = 256
_INVALID_AMINO_ACID_INDEX = -1
_AMINO_ACID_INDEX_LOOKUP = np.full(
    _ASCII_LOOKUP_SIZE,
    _INVALID_AMINO_ACID_INDEX,
    dtype=np.int16,
)
for _row_index, _amino_acid in enumerate(AMINO_ACIDS):
    _AMINO_ACID_INDEX_LOOKUP[ord(_amino_acid)] = _row_index


@dataclass(frozen=True, slots=True)
class MotifScoringResult:
    """Motif score matrices and window metadata for one scoring run."""

    motif_scores: pd.DataFrame
    motif_sizes: pd.Series
    sequence_windows: pd.Series


def build_motif_library(
    *,
    kinase_substrate_map: pd.DataFrame,
    site_sequences: pd.Series,
    flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
) -> tuple[dict[str, pd.DataFrame], pd.Series]:
    """Build per-kinase motif frequency matrices from reference sequences."""

    if flank_size < 0:
        raise ValueError("flank_size must be >= 0")

    if site_sequences.empty:
        return {}, pd.Series(dtype=float, name="motif_size")

    sequence_index = set(site_sequences.index.astype(str))
    frequency_matrices: dict[str, pd.DataFrame] = {}
    motif_sizes: dict[str, float] = {}
    for kinase, grouped in kinase_substrate_map.groupby("kinase", sort=False):
        sites = list(dict.fromkeys(grouped.loc[:, "substrate_site"].astype(str)))
        sequence_sites = [site for site in sites if site in sequence_index]
        if not sequence_sites:
            continue
        windows = _coerce_sequence_series(
            site_sequences.loc[sequence_sites],
            flank_size=flank_size,
        )
        windows = windows.dropna()
        if windows.empty:
            continue
        frequency_matrices[str(kinase)] = _build_frequency_matrix_from_windows(
            windows,
            context=f"kinase={kinase}",
        )
        motif_sizes[str(kinase)] = float(len(windows))

    size_series = pd.Series(motif_sizes, dtype=float, name="motif_size")
    size_series.index.name = "kinase"
    return frequency_matrices, size_series


def build_motif_library_from_sequences(
    *,
    motif_sequences: Mapping[str, Sequence[str]],
    flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
) -> tuple[dict[str, pd.DataFrame], pd.Series]:
    """Build motif frequency matrices from explicit per-kinase sequences."""

    if flank_size < 0:
        raise ValueError("flank_size must be >= 0")
    frequency_matrices: dict[str, pd.DataFrame] = {}
    motif_sizes: dict[str, float] = {}
    for kinase, sequences in motif_sequences.items():
        windows = _coerce_sequence_series(sequences, flank_size=flank_size).dropna()
        if windows.empty:
            continue
        frequency_matrices[str(kinase)] = _build_frequency_matrix_from_windows(
            windows,
            context=f"kinase={kinase}",
        )
        motif_sizes[str(kinase)] = float(len(windows))
    size_series = pd.Series(motif_sizes, dtype=float, name="motif_size")
    size_series.index.name = "kinase"
    return frequency_matrices, size_series


def score_phosphosite_motifs(
    *,
    site_sequences: Mapping[str, str] | Sequence[str] | pd.Series,
    motif_frequency_matrices: Mapping[str, pd.DataFrame],
    motif_sizes: pd.Series,
    site_index: Sequence[str] | None = None,
    min_motif_size: int = 1,
    flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
) -> MotifScoringResult:
    """Score phosphosite sequences against per-kinase motif frequency matrices.

    Runtime scales with scored sites, eligible kinases, and motif window width.
    """

    if min_motif_size < 1:
        raise ValueError("min_motif_size must be >= 1")

    windows = _coerce_sequence_series(
        site_sequences,
        site_index=site_index,
        flank_size=flank_size,
    )
    kinases = [
        kinase
        for kinase in motif_frequency_matrices
        if kinase in motif_sizes.index
        and float(motif_sizes.loc[kinase]) >= float(min_motif_size)
    ]
    motif_scores = pd.DataFrame(
        np.nan,
        index=windows.index.copy(),
        columns=kinases,
        dtype=float,
    )
    if kinases:
        first_matrix = motif_frequency_matrices[kinases[0]]
        encoded_windows = _encode_sequence_positions(
            windows,
            width=int(first_matrix.shape[1]),
        )
        for kinase in kinases:
            frequency_values = _coerce_frequency_matrix(
                motif_frequency_matrices[kinase]
            ).to_numpy(dtype=float, copy=False)
            motif_scores.loc[:, kinase] = _score_encoded_sequences(
                encoded_sequences=encoded_windows,
                frequency_values=frequency_values,
            )
    motif_scores = minmax_scale_columns(motif_scores)

    selected_sizes = motif_sizes.loc[kinases].astype(float).copy()
    selected_sizes.index.name = "kinase"
    return MotifScoringResult(
        motif_scores=motif_scores,
        motif_sizes=selected_sizes,
        sequence_windows=windows,
    )


def minmax_scale_columns(mat: pd.DataFrame) -> pd.DataFrame:
    """Apply column-wise min-max scaling used by baseline motif scoring."""

    scaled = mat.astype(float).copy()
    for column in scaled.columns:
        values = scaled.loc[:, column]
        min_value = float(values.min())
        max_value = float(values.max())
        denominator = max_value - min_value
        if denominator == 0.0:
            scaled.loc[:, column] = np.nan
        else:
            scaled.loc[:, column] = (values - min_value) / denominator
    return scaled


def _build_frequency_matrix_from_windows(
    windows: pd.Series,
    *,
    context: str,
) -> pd.DataFrame:
    if windows.empty:
        raise ValueError(f"{context} must contain at least one sequence")
    width = len(str(windows.iloc[0]))
    if any(len(str(window)) != width for window in windows):
        raise ValueError(f"{context} must contain same-length sequence windows")

    frequency_values = np.zeros((len(AMINO_ACIDS), width), dtype=float)
    encoded_windows = _encode_sequence_positions(windows, width=width)
    valid_rows, valid_cols = np.nonzero(encoded_windows >= 0)
    if valid_rows.size > 0:
        np.add.at(
            frequency_values,
            (encoded_windows[valid_rows, valid_cols], valid_cols),
            1.0,
        )
    frequency_values /= float(len(windows))
    return pd.DataFrame(
        frequency_values,
        index=list(AMINO_ACIDS),
        columns=[f"p{i}" for i in range(1, width + 1)],
        dtype=float,
    )


def _coerce_frequency_matrix(frequency_mat: pd.DataFrame) -> pd.DataFrame:
    matrix = frequency_mat.astype(float).copy()
    matrix = matrix.reindex(index=list(AMINO_ACIDS), fill_value=0.0)
    matrix.columns = [f"p{i}" for i in range(1, matrix.shape[1] + 1)]
    return matrix


def _coerce_sequence_series(
    seqs: Mapping[str, str] | Sequence[str] | pd.Series,
    *,
    site_index: Sequence[str] | None = None,
    flank_size: int | None = DEFAULT_MOTIF_FLANK_SIZE,
) -> pd.Series:
    if isinstance(seqs, pd.Series):
        series = seqs.copy()
    elif isinstance(seqs, Mapping):
        series = pd.Series(dict(seqs), dtype=object)
    else:
        seq_list = list(seqs)
        if site_index is None:
            series = pd.Series(seq_list, dtype=object)
        else:
            if len(seq_list) != len(site_index):
                raise ValueError("site_index must have same length as seqs")
            series = pd.Series(seq_list, index=list(site_index), dtype=object)

    if site_index is not None:
        missing = [site for site in site_index if site not in series.index]
        if missing:
            raise ValueError(
                f"site_sequences missing entries for {', '.join(map(str, missing[:5]))}"
            )
        series = series.loc[list(site_index)]

    return series.map(lambda value: _extract_sequence_window(value, flank_size))


def _extract_sequence_window(value: object, flank_size: int | None) -> object:
    if value is None or pd.isna(value):
        return np.nan
    sequence = str(value).upper()
    if sequence == "":
        return np.nan
    if flank_size is None:
        return sequence
    window_size = (2 * flank_size) + 1
    if len(sequence) <= window_size:
        return sequence
    mid = len(sequence) // 2
    start = mid - flank_size
    stop = mid + flank_size + 1
    return sequence[start:stop]


def _encode_sequence_positions(
    sequences: Sequence[object] | pd.Series,
    *,
    width: int,
) -> np.ndarray:
    encoded = np.full(
        (len(sequences), width),
        _INVALID_AMINO_ACID_INDEX,
        dtype=np.int16,
    )
    if width == 0:
        return encoded
    for row_index, sequence in enumerate(sequences):
        if sequence is None or pd.isna(sequence):
            continue
        text = str(sequence)[:width]
        if text == "":
            continue
        code_points = np.fromiter(
            (ord(character) for character in text),
            dtype=np.int32,
            count=len(text),
        )
        valid_ascii = code_points < _ASCII_LOOKUP_SIZE
        if not np.any(valid_ascii):
            continue
        valid_positions = np.flatnonzero(valid_ascii)
        encoded[row_index, valid_positions] = _AMINO_ACID_INDEX_LOOKUP[
            code_points[valid_positions]
        ]
    return encoded


def _score_encoded_sequences(
    *,
    encoded_sequences: np.ndarray,
    frequency_values: np.ndarray,
) -> np.ndarray:
    if encoded_sequences.size == 0:
        return np.zeros(encoded_sequences.shape[0], dtype=float)
    valid_mask = encoded_sequences >= 0
    safe_indices = np.where(valid_mask, encoded_sequences, 0)
    position_indices = np.arange(encoded_sequences.shape[1])
    position_scores = frequency_values[safe_indices, position_indices]
    position_scores = np.where(valid_mask, position_scores, 0.0)
    return position_scores.sum(axis=1, dtype=float)


__all__ = [
    "DEFAULT_MOTIF_FLANK_SIZE",
    "MotifScoringResult",
    "build_motif_library",
    "build_motif_library_from_sequences",
    "minmax_scale_columns",
    "score_phosphosite_motifs",
]
