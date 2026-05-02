from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from phospy.prediction.motif_scoring.models import (
    _AMINO_ACID_INDEX_LOOKUP,
    _ASCII_LOOKUP_SIZE,
    _INVALID_AMINO_ACID_INDEX,
    AMINO_ACIDS,
    _is_missing_scalar,
)


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
        index=pd.Index(AMINO_ACIDS),
        columns=pd.Index([f"p{i}" for i in range(1, width + 1)]),
        dtype=float,
    )


def _coerce_frequency_matrix(frequency_mat: pd.DataFrame) -> pd.DataFrame:
    matrix = frequency_mat.astype(float).copy()
    matrix = matrix.reindex(index=list(AMINO_ACIDS), fill_value=0.0)
    matrix.columns = [f"p{i}" for i in range(1, matrix.shape[1] + 1)]
    return matrix


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
        if _is_missing_scalar(sequence):
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


def _require_fully_supported_encoded_sequences(
    encoded_sequences: np.ndarray,
    *,
    context: str,
) -> None:
    if (encoded_sequences < 0).any():
        raise ValueError(f"{context} contains unsupported residues")
