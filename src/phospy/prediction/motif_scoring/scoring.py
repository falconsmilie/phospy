from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from phospy.prediction.motif_scoring.frequency_matrices import (
    _coerce_frequency_matrix,
    _encode_sequence_positions,
    _require_fully_supported_encoded_sequences,
    _score_encoded_sequences,
)
from phospy.prediction.motif_scoring.models import (
    DEFAULT_MOTIF_FLANK_SIZE,
    SEQUENCE_SEMANTICS_CENTRED_SEQUENCE,
    SEQUENCE_SEMANTICS_CENTRED_WINDOW,
    MotifLibraryValidationResult,
    MotifScoringResult,
    SequenceSemantics,
    _extract_sequence_window,
    _is_missing_scalar,
    _normalize_sequence_value,
)
from phospy.prediction.motif_scoring.scaling import minmax_scale_columns
from phospy.prediction.sequence_validation import (
    SEQUENCE_VALIDATION_STATUS_VALID,
    MotifSequenceValidator,
    SequenceValidationInput,
    SequenceValidationResult,
)


def score_phosphosite_motifs(
    *,
    site_sequences: Mapping[str, str] | Sequence[str] | pd.Series,
    motif_frequency_matrices: Mapping[str, pd.DataFrame],
    motif_sizes: pd.Series,
    site_index: Sequence[str] | None = None,
    min_motif_size: int = 1,
    flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
    sequence_semantics: SequenceSemantics = SEQUENCE_SEMANTICS_CENTRED_WINDOW,
    library_validation: MotifLibraryValidationResult | None = None,
) -> MotifScoringResult:
    """Score phosphosite sequences against per-kinase motif frequency matrices.

    `sequence_semantics="centred_window"` (default) requires each sequence to be a
    centred phosphosite window matching the motif width exactly.
    `sequence_semantics="centred_sequence"` accepts centred odd-length sequences
    with length >= motif width and centre-extracts the scoring window.

    Runtime scales with scored sites, eligible kinases, and motif window width.
    """

    if min_motif_size < 1:
        raise ValueError("min_motif_size must be >= 1")

    raw_sequences = _coerce_sequence_series(
        site_sequences,
        site_index=site_index,
    )
    kinases = [
        kinase
        for kinase in motif_frequency_matrices
        if kinase in motif_sizes.index
        and float(motif_sizes.loc[kinase]) >= float(min_motif_size)
    ]
    expected_window_size = _resolve_expected_window_size(
        motif_frequency_matrices=motif_frequency_matrices,
        kinases=kinases,
        flank_size=flank_size,
    )
    validation = _validate_sequence_windows(
        raw_sequences,
        expected_window_size=expected_window_size,
        sequence_semantics=sequence_semantics,
    )
    windows = _materialize_scoring_windows(
        raw_sequences=raw_sequences,
        validation=validation,
        expected_window_size=expected_window_size,
        flank_size=flank_size,
        sequence_semantics=sequence_semantics,
    )

    motif_scores = pd.DataFrame(
        np.nan,
        index=raw_sequences.index.copy(),
        columns=pd.Index(kinases),
        dtype=float,
    )
    if kinases and validation.valid_sequences > 0:
        valid_site_ids = [
            row.site_id
            for row in validation.rows
            if row.status == SEQUENCE_VALIDATION_STATUS_VALID
        ]
        valid_windows = windows.loc[valid_site_ids]
        first_matrix = motif_frequency_matrices[kinases[0]]
        encoded_windows = _encode_sequence_positions(
            valid_windows,
            width=int(first_matrix.shape[1]),
        )
        _require_fully_supported_encoded_sequences(
            encoded_windows,
            context="motif scoring",
        )
        for kinase in kinases:
            frequency_values = _coerce_frequency_matrix(
                motif_frequency_matrices[kinase]
            ).to_numpy(dtype=float, copy=False)
            motif_scores.loc[valid_windows.index, kinase] = _score_encoded_sequences(
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
        sequence_validation=validation,
        library_validation=library_validation,
    )


def _coerce_sequence_series(
    seqs: Mapping[str, str] | Sequence[str] | pd.Series,
    *,
    site_index: Sequence[str] | None = None,
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

    return series.map(_normalize_sequence_value)


def _resolve_expected_window_size(
    *,
    motif_frequency_matrices: Mapping[str, pd.DataFrame],
    kinases: Sequence[str],
    flank_size: int,
) -> int:
    if flank_size < 0:
        raise ValueError("flank_size must be >= 0")
    fallback = (2 * flank_size) + 1
    if not kinases:
        return fallback
    widths = {
        int(_coerce_frequency_matrix(motif_frequency_matrices[kinase]).shape[1])
        for kinase in kinases
    }
    if len(widths) != 1:
        raise ValueError("motif_frequency_matrices must share one sequence width")
    width = next(iter(widths))
    if width <= 0:
        raise ValueError("motif_frequency_matrices sequence width must be > 0")
    return width


def _validate_sequence_windows(
    windows: pd.Series,
    *,
    expected_window_size: int,
    sequence_semantics: SequenceSemantics,
) -> SequenceValidationResult:
    _require_supported_sequence_semantics(sequence_semantics)
    validator = MotifSequenceValidator(
        expected_window_size=expected_window_size,
        require_phospho_centre_residue=True,
        enforce_site_identity_format=True,
        window_length_policy=(
            "exact"
            if sequence_semantics == SEQUENCE_SEMANTICS_CENTRED_WINDOW
            else "centred_superset"
        ),
        residue_validation_scope=(
            "full_sequence"
            if sequence_semantics == SEQUENCE_SEMANTICS_CENTRED_WINDOW
            else "centre_window"
        ),
    )
    rows = [
        SequenceValidationInput(
            site_id=str(site_id),
            site_sequence=sequence,
            site_identity=str(site_id),
        )
        for site_id, sequence in windows.items()
    ]
    return validator.run(rows=rows)


def _materialize_scoring_windows(
    *,
    raw_sequences: pd.Series,
    validation: SequenceValidationResult,
    expected_window_size: int,
    flank_size: int,
    sequence_semantics: SequenceSemantics,
) -> pd.Series:
    windows = pd.Series(
        [None] * len(raw_sequences),
        index=raw_sequences.index.copy(),
        dtype=object,
    )
    for row in validation.rows:
        if row.sequence is None:
            continue
        if row.status != SEQUENCE_VALIDATION_STATUS_VALID:
            windows.loc[row.site_id] = row.sequence
            continue
        if sequence_semantics == SEQUENCE_SEMANTICS_CENTRED_WINDOW:
            windows.loc[row.site_id] = row.sequence
            continue
        extracted = _extract_sequence_window(row.sequence, flank_size)
        if _is_missing_scalar(extracted):
            windows.loc[row.site_id] = row.sequence
            continue
        if len(str(extracted)) != expected_window_size:
            windows.loc[row.site_id] = row.sequence
            continue
        windows.loc[row.site_id] = str(extracted)
    return windows


def _require_supported_sequence_semantics(
    sequence_semantics: SequenceSemantics,
) -> None:
    if sequence_semantics not in {
        SEQUENCE_SEMANTICS_CENTRED_WINDOW,
        SEQUENCE_SEMANTICS_CENTRED_SEQUENCE,
    }:
        raise ValueError(
            "sequence_semantics must be 'centred_window' or 'centred_sequence'"
        )
