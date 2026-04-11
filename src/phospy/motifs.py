from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .errors import (
    InputCompatibilityError,
    TableSchemaError,
)
from .validation.values.numeric import validate_non_negative_int, validate_positive_int

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

_ASCII_LOOKUP_SIZE = 256
_INVALID_AMINO_ACID_INDEX = -1
_AMINO_ACID_INDEX_LOOKUP = np.full(
    _ASCII_LOOKUP_SIZE,
    _INVALID_AMINO_ACID_INDEX,
    dtype=np.int16,
)
for _row_index, _amino_acid in enumerate(AMINO_ACIDS):
    _AMINO_ACID_INDEX_LOOKUP[ord(_amino_acid)] = _row_index

_ALLOWED_SEQUENCE_CHARACTERS: frozenset[str] = frozenset(AMINO_ACIDS) | frozenset({"_"})


@dataclass(slots=True)
class MotifScoringResult:
    """Motif scoring tables and sequence windows for one scoring run."""

    motif_scores: pd.DataFrame
    motif_sizes: pd.Series
    sequence_windows: pd.Series


@dataclass(slots=True)
class ValidatedMotifLibrary:
    """Trusted validated motif library bundle used by motif scoring setup.

    The contained frequency matrices and size series remain mutable pandas data,
    so this wrapper is not presented as an immutable value object.
    """

    motif_frequency_matrices: dict[str, pd.DataFrame]
    motif_sizes: pd.Series


class KinaseMotifScorer:
    """Score phosphosite sequences against per-kinase motif frequency matrices."""

    def __init__(
        self,
        motif_frequency_matrices: Mapping[str, pd.DataFrame],
        motif_sizes: pd.Series,
        flank_size: int = 7,
    ) -> None:
        validate_non_negative_int(flank_size, name="flank_size")
        if not motif_frequency_matrices:
            msg = "motif_frequency_matrices must not be empty"
            raise InputCompatibilityError(msg)

        self.motif_frequency_matrices = {
            kinase: _coerce_frequency_matrix(matrix)
            for kinase, matrix in motif_frequency_matrices.items()
        }
        motif_widths = {
            _motif_matrix_width(matrix)
            for matrix in self.motif_frequency_matrices.values()
        }
        if len(motif_widths) != 1:
            msg = "All motif frequency matrices must have the same window width"
            raise InputCompatibilityError(msg)
        self._motif_width = motif_widths.pop()
        self._motif_frequency_values = {
            kinase: matrix.to_numpy(dtype=float, copy=False)
            for kinase, matrix in self.motif_frequency_matrices.items()
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
        validated_library = build_validated_motif_library(
            motif_sequences=motif_sequences,
            flank_size=flank_size,
        )
        return cls(
            motif_frequency_matrices=validated_library.motif_frequency_matrices,
            motif_sizes=validated_library.motif_sizes,
            flank_size=flank_size,
        )

    def score_sequences(
        self,
        seqs: Mapping[str, str] | Sequence[str] | pd.Series,
        site_index: Sequence[str] | None = None,
        min_motif_size: int = 1,
    ) -> MotifScoringResult:
        validate_positive_int(min_motif_size, name="min_motif_size")
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
        if kinases:
            encoded_windows = _encode_sequence_positions(
                windows,
                width=self._motif_width,
            )
            for kinase in kinases:
                motif_scores.loc[:, kinase] = _score_encoded_sequences(
                    encoded_sequences=encoded_windows,
                    frequency_values=self._motif_frequency_values[kinase],
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

    validate_non_negative_int(flank_size, name="flank_size")
    windows = _coerce_sequence_series(substrates_seq, flank_size=flank_size)
    return _build_frequency_matrix_from_windows(
        windows,
        context="substrates_seq",
    )


def build_validated_motif_library(
    motif_sequences: Mapping[str, Sequence[str]],
    *,
    flank_size: int = 7,
    context: str = "motif_sequences",
) -> ValidatedMotifLibrary:
    validate_non_negative_int(flank_size, name="flank_size")

    matrices: dict[str, pd.DataFrame] = {}
    sizes: dict[str, float] = {}
    widths: set[int] = set()

    for kinase, sequences in motif_sequences.items():
        windows = _coerce_sequence_series(sequences, flank_size=flank_size)
        try:
            matrices[kinase] = _build_frequency_matrix_from_windows(
                windows,
                context=f"{context} for kinase {kinase}",
            )
        except TableSchemaError as error:
            if "same window length" in str(error):
                msg = f"{context} for kinase {kinase} must use a consistent sequence width"
                raise InputCompatibilityError(msg) from error
            raise
        sizes[kinase] = float(len(windows))
        widths.add(matrices[kinase].shape[1])

    if len(widths) > 1:
        msg = f"{context} must use the same sequence width across kinases"
        raise InputCompatibilityError(msg)

    return ValidatedMotifLibrary(
        motif_frequency_matrices=matrices,
        motif_sizes=pd.Series(sizes, dtype=float),
    )


def frequency_scoring(
    sequence_list: Sequence[str] | pd.Series,
    frequency_mat: pd.DataFrame,
) -> pd.Series:
    """Score phosphosite windows against a single motif frequency matrix."""

    frequency_mat = _coerce_frequency_matrix(frequency_mat)
    sequences = _coerce_sequence_series(sequence_list, flank_size=None)
    frequency_values = frequency_mat.to_numpy(dtype=float, copy=False)
    encoded_sequences = _encode_sequence_positions(
        sequences,
        width=frequency_values.shape[1],
    )
    score_values = _score_encoded_sequences(
        encoded_sequences=encoded_sequences,
        frequency_values=frequency_values,
    )
    return pd.Series(score_values, index=sequences.index.copy(), dtype=float)


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


def _build_frequency_matrix_from_windows(
    windows: pd.Series,
    *,
    context: str,
) -> pd.DataFrame:
    if windows.empty:
        msg = f"{context} must contain at least one sequence"
        raise TableSchemaError(msg)

    width = len(str(windows.iloc[0]))
    if any(len(str(window)) != width for window in windows):
        msg = "All sequences must have the same window length"
        raise TableSchemaError(msg)

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


def _encode_sequence_positions(
    sequences: Sequence[object] | pd.Series,
    width: int,
) -> np.ndarray:
    """Encode sequence characters into amino-acid row indices for vectorised scoring."""
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
    encoded_sequences: np.ndarray,
    frequency_values: np.ndarray,
) -> np.ndarray:
    """Score encoded sequence windows against a frequency matrix using NumPy indexing."""
    if encoded_sequences.size == 0:
        return np.zeros(encoded_sequences.shape[0], dtype=float)

    valid_mask = encoded_sequences >= 0
    safe_indices = np.where(valid_mask, encoded_sequences, 0)
    position_indices = np.arange(encoded_sequences.shape[1])
    position_scores = frequency_values[safe_indices, position_indices]
    position_scores = np.where(valid_mask, position_scores, 0.0)
    return position_scores.sum(axis=1, dtype=float)


def _motif_matrix_width(frequency_mat: pd.DataFrame) -> int:
    return int(frequency_mat.shape[1])


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
        series = seqs
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
        series = series.loc[list(site_index)]

    return series.map(lambda value: _extract_sequence_window(value, flank_size))


def _extract_sequence_window(value: object, flank_size: int | None) -> str:
    if value is None or pd.isna(value):
        return np.nan  # type: ignore[return-value]

    sequence = str(value).upper()
    if sequence == "":
        return "_"

    if flank_size is None:
        return _validate_sequence_window(sequence)

    window_size = (2 * flank_size) + 1
    if len(sequence) <= window_size:
        return _validate_sequence_window(sequence)

    mid = len(sequence) // 2
    start = mid - flank_size
    stop = mid + flank_size + 1
    return _validate_sequence_window(sequence[start:stop])


def _validate_sequence_window(sequence: str) -> str:
    invalid_characters = sorted(
        {
            character
            for character in sequence
            if character not in _ALLOWED_SEQUENCE_CHARACTERS
        }
    )
    if invalid_characters:
        invalid_text = ", ".join(repr(character) for character in invalid_characters)
        msg = (
            "sequence contains invalid amino-acid characters: "
            f"{invalid_text}; allowed characters are the 20 standard amino acids and '_'"
        )
        raise TableSchemaError(msg)
    return sequence
