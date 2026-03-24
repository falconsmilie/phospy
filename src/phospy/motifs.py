from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .validation.errors import (
    InputCompatibilityError,
    PhospyValidationError,
    TableSchemaError,
)

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


@dataclass(slots=True)
class MotifScoringResult:
    motif_scores: pd.DataFrame
    motif_sizes: pd.Series
    sequence_windows: pd.Series


class KinaseMotifScorer:
    """Score phosphosite sequences against per-kinase motif frequency matrices."""

    def __init__(
        self,
        motif_frequency_matrices: Mapping[str, pd.DataFrame],
        motif_sizes: pd.Series,
        flank_size: int = 7,
    ) -> None:
        _validate_non_negative_int(flank_size, name="flank_size")
        if not motif_frequency_matrices:
            msg = "motif_frequency_matrices must not be empty"
            raise InputCompatibilityError(msg)

        self.motif_frequency_matrices = {
            kinase: _coerce_frequency_matrix(matrix)
            for kinase, matrix in motif_frequency_matrices.items()
        }
        self.motif_sizes = motif_sizes.astype(float).copy()
        self.flank_size = flank_size

        missing = [
            kinase
            for kinase in self.motif_frequency_matrices
            if kinase not in self.motif_sizes.index
        ]
        if missing:
            msg = f"motif_sizes is missing entries for: {', '.join(missing)}"
            raise InputCompatibilityError(msg)

    @classmethod
    def from_substrate_sequences(
        cls,
        motif_sequences: Mapping[str, Sequence[str]],
        flank_size: int = 7,
    ) -> KinaseMotifScorer:
        matrices: dict[str, pd.DataFrame] = {}
        sizes: dict[str, float] = {}
        for kinase, sequences in motif_sequences.items():
            matrices[kinase] = create_frequency_matrix(
                sequences,
                flank_size=flank_size,
            )
            sizes[kinase] = float(len(sequences))
        return cls(
            motif_frequency_matrices=matrices,
            motif_sizes=pd.Series(sizes, dtype=float),
            flank_size=flank_size,
        )

    def score_sequences(
        self,
        seqs: Mapping[str, str] | Sequence[str] | pd.Series,
        site_index: Sequence[str] | None = None,
        min_motif_size: int = 1,
    ) -> MotifScoringResult:
        _validate_positive_int(min_motif_size, name="min_motif_size")
        windows = _coerce_sequence_series(
            seqs=seqs,
            site_index=site_index,
            flank_size=self.flank_size,
        )

        kinases = [
            kinase
            for kinase in self.motif_frequency_matrices
            if float(self.motif_sizes.loc[kinase]) >= float(min_motif_size)
        ]
        motif_scores = pd.DataFrame(
            np.nan,
            index=windows.index.copy(),
            columns=kinases,
            dtype=float,
        )
        for kinase in kinases:
            motif_scores.loc[:, kinase] = frequency_scoring(
                sequence_list=windows,
                frequency_mat=self.motif_frequency_matrices[kinase],
            )
        motif_scores = minmax_scale_columns(motif_scores)

        motif_sizes = self.motif_sizes.loc[kinases].astype(float).copy()
        motif_sizes.index.name = "kinase"
        return MotifScoringResult(
            motif_scores=motif_scores,
            motif_sizes=motif_sizes,
            sequence_windows=windows,
        )


def create_frequency_matrix(
    substrates_seq: Sequence[str] | pd.Series,
    flank_size: int = 7,
) -> pd.DataFrame:
    """Create an amino-acid frequency matrix from substrate sequences."""

    _validate_non_negative_int(flank_size, name="flank_size")
    windows = _coerce_sequence_series(substrates_seq, flank_size=flank_size)
    if windows.empty:
        msg = "substrates_seq must contain at least one sequence"
        raise TableSchemaError(msg)

    width = len(str(windows.iloc[0]))
    if any(len(str(window)) != width for window in windows):
        msg = "All sequences must have the same window length"
        raise TableSchemaError(msg)

    frequency_mat = pd.DataFrame(
        0.0,
        index=list(AMINO_ACIDS),
        columns=[f"p{i}" for i in range(1, width + 1)],
        dtype=float,
    )

    for window in windows:
        for position, aa in enumerate(str(window), start=1):
            if aa == "_":
                continue
            if aa in frequency_mat.index:
                frequency_mat.loc[aa, f"p{position}"] += 1.0

    frequency_mat /= float(len(windows))
    return frequency_mat


def frequency_scoring(
    sequence_list: Sequence[str] | pd.Series,
    frequency_mat: pd.DataFrame,
) -> pd.Series:
    """Score phosphosite windows against a single motif frequency matrix."""

    frequency_mat = _coerce_frequency_matrix(frequency_mat)
    sequences = _coerce_sequence_series(sequence_list, flank_size=None)
    max_width = frequency_mat.shape[1]

    scores = pd.Series(0.0, index=sequences.index.copy(), dtype=float)
    for site, sequence in sequences.items():
        if pd.isna(sequence):
            scores.loc[site] = 0.0
            continue

        score = 0.0
        for position, aa in enumerate(str(sequence), start=1):
            if position > max_width:
                break
            if aa in frequency_mat.index:
                score += float(frequency_mat.loc[aa, f"p{position}"])
        scores.loc[site] = score
    return scores


def score_phosphosite_motifs(
    seqs: Mapping[str, str] | Sequence[str] | pd.Series,
    motif_frequency_matrices: Mapping[str, pd.DataFrame],
    motif_sizes: pd.Series,
    site_index: Sequence[str] | None = None,
    min_motif_size: int = 1,
    flank_size: int = 7,
) -> MotifScoringResult:
    scorer = KinaseMotifScorer(
        motif_frequency_matrices=motif_frequency_matrices,
        motif_sizes=motif_sizes,
        flank_size=flank_size,
    )
    return scorer.score_sequences(
        seqs=seqs,
        site_index=site_index,
        min_motif_size=min_motif_size,
    )


def minmax_scale_columns(mat: pd.DataFrame) -> pd.DataFrame:
    """Apply PhosR-style column-wise min-max scaling."""

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


def _coerce_frequency_matrix(frequency_mat: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(frequency_mat, pd.DataFrame):
        msg = "frequency_mat must be a pandas DataFrame"
        raise TypeError(msg)

    matrix = frequency_mat.astype(float).copy()
    matrix = matrix.reindex(index=list(AMINO_ACIDS), fill_value=0.0)
    matrix.columns = [f"p{i}" for i in range(1, matrix.shape[1] + 1)]
    return matrix


def _coerce_sequence_series(
    seqs: Mapping[str, str] | Sequence[str] | pd.Series,
    site_index: Sequence[str] | None = None,
    flank_size: int | None = 7,
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
                msg = "site_index must have the same length as seqs"
                raise TableSchemaError(msg)
            series = pd.Series(seq_list, index=list(site_index), dtype=object)

    if site_index is not None:
        missing = [site for site in site_index if site not in series.index]
        if missing:
            msg = f"seqs is missing entries for: {', '.join(missing)}"
            raise TableSchemaError(msg)
        series = series.loc[list(site_index)].copy()

    return series.map(lambda value: _extract_sequence_window(value, flank_size))


def _extract_sequence_window(value: object, flank_size: int | None) -> str:
    if value is None or pd.isna(value):
        return np.nan  # type: ignore[return-value]

    sequence = str(value)
    if sequence == "":
        return "_"

    if flank_size is None:
        return sequence

    window_size = (2 * flank_size) + 1
    if len(sequence) <= window_size:
        return sequence

    mid = len(sequence) // 2
    start = mid - flank_size
    stop = mid + flank_size + 1
    return sequence[start:stop]


def _validate_non_negative_int(value: int, name: str) -> None:
    if value < 0:
        raise PhospyValidationError(f"{name} must be at least 0")


def _validate_positive_int(value: int, name: str) -> None:
    if value < 1:
        raise PhospyValidationError(f"{name} must be at least 1")
