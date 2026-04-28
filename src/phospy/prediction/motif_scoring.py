"""Motif scoring kernels used by kinase workflow scoring."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

from phospy.prediction.sequence_validation import (
    SEQUENCE_VALIDATION_STATUS_INVALID_SITE_ID,
    SEQUENCE_VALIDATION_STATUS_MISSING_SEQUENCE,
    SEQUENCE_VALIDATION_STATUS_NON_PHOSPHO_CENTRE_RESIDUE,
    SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE,
    SEQUENCE_VALIDATION_STATUS_SHORT_SEQUENCE,
    SEQUENCE_VALIDATION_STATUS_SITE_RESIDUE_MISMATCH,
    SEQUENCE_VALIDATION_STATUS_UNSUPPORTED_RESIDUE_CHARACTER,
    SEQUENCE_VALIDATION_STATUS_VALID,
    SUPPORTED_AMINO_ACIDS,
    MotifSequenceValidator,
    SequenceValidationInput,
    SequenceValidationResult,
    SequenceValidationStatus,
)

AMINO_ACIDS: tuple[str, ...] = SUPPORTED_AMINO_ACIDS
DEFAULT_MOTIF_FLANK_SIZE = 7
SEQUENCE_SEMANTICS_CENTRED_WINDOW = "centred_window"
SEQUENCE_SEMANTICS_CENTRED_SEQUENCE = "centred_sequence"
SequenceSemantics = Literal["centred_window", "centred_sequence"]
_MOTIF_LIBRARY_VALIDATION_ATTR = "motif_library_validation"
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
    sequence_validation: SequenceValidationResult
    library_validation: MotifLibraryValidationResult | None = None


@dataclass(frozen=True, slots=True)
class MotifLibraryValidationRow:
    """Per-reference motif-library validation outcome."""

    reference_id: str
    site_id: str | None
    kinase: str
    sequence: str | None
    status: SequenceValidationStatus
    reason: str | None
    expected_centre_residue: str | None
    observed_centre_residue: str | None
    sequence_length: int | None


@dataclass(frozen=True, slots=True)
class MotifLibraryValidationResult:
    """Structured diagnostics for motif-library/reference sequence validation."""

    total_reference_sequences: int
    accepted_reference_sequences: int
    excluded_reference_sequences: int
    missing_sequences: int
    short_sequences: int
    off_centre_sequences: int
    site_residue_mismatches: int
    invalid_site_ids: int
    non_phospho_centre_residues: int
    unsupported_residue_characters: int
    sequences_excluded_from_motif_profile_construction: int
    expected_window_size: int
    supported_amino_acids: tuple[str, ...]
    accepted_window_length_policy: str
    unsupported_residue_policy: str
    excluded_reference_ids: tuple[str, ...]
    rows: tuple[MotifLibraryValidationRow, ...]

    def summary(self) -> dict[str, object]:
        """Return compact diagnostics for run-level reporting/provenance."""

        return {
            "reference_sequences_provided": self.total_reference_sequences,
            "reference_sequences_accepted": self.accepted_reference_sequences,
            "reference_sequences_excluded": self.excluded_reference_sequences,
            "excluded_missing_sequence": self.missing_sequences,
            "excluded_short_window": self.short_sequences,
            "excluded_unsupported_residue": self.unsupported_residue_characters,
            "excluded_off_centre_residue": self.off_centre_sequences,
            "excluded_site_residue_mismatch": self.site_residue_mismatches,
            "excluded_invalid_site_id": self.invalid_site_ids,
            "excluded_non_phospho_centre_residue": (self.non_phospho_centre_residues),
            "sequences_excluded_from_motif_profile_construction": (
                self.sequences_excluded_from_motif_profile_construction
            ),
            "expected_window_size": self.expected_window_size,
            "accepted_window_length_policy": self.accepted_window_length_policy,
            "unsupported_residue_policy": self.unsupported_residue_policy,
            "supported_amino_acids": list(self.supported_amino_acids),
        }

    def to_frame(self) -> pd.DataFrame:
        """Return row-level library validation provenance as a compact table."""

        columns = [
            "reference_id",
            "site_id",
            "kinase",
            "sequence",
            "status",
            "reason",
            "observed_centre_residue",
            "expected_centre_residue",
            "sequence_length",
        ]
        return pd.DataFrame(
            [
                {
                    "reference_id": row.reference_id,
                    "site_id": row.site_id,
                    "kinase": row.kinase,
                    "sequence": row.sequence,
                    "status": row.status,
                    "reason": row.reason,
                    "observed_centre_residue": row.observed_centre_residue,
                    "expected_centre_residue": row.expected_centre_residue,
                    "sequence_length": row.sequence_length,
                }
                for row in self.rows
            ],
            columns=columns,
        )


@dataclass(frozen=True, slots=True)
class _LibraryCandidate:
    reference_id: str
    site_id: str | None
    kinase: str
    sequence_input: object


@dataclass(frozen=True, slots=True)
class ExplicitMotifSequence:
    """Structured explicit motif-library sequence metadata."""

    reference_id: str
    site_id: str | None
    kinase: str
    sequence: object


def build_motif_library(
    *,
    kinase_substrate_map: pd.DataFrame,
    site_sequences: pd.Series,
    flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
) -> tuple[dict[str, pd.DataFrame], pd.Series]:
    """Build per-kinase motif frequency matrices from reference sequences."""

    if flank_size < 0:
        raise ValueError("flank_size must be >= 0")

    sequence_lookup = {
        str(site_id): sequence for site_id, sequence in site_sequences.items()
    }
    candidates: list[_LibraryCandidate] = []
    for kinase, grouped in kinase_substrate_map.groupby("kinase", sort=False):
        sites = list(dict.fromkeys(grouped.loc[:, "substrate_site"].astype(str)))
        for reference_id in sites:
            candidates.append(
                _LibraryCandidate(
                    reference_id=reference_id,
                    site_id=reference_id,
                    kinase=str(kinase),
                    sequence_input=_normalize_sequence_value(
                        sequence_lookup.get(reference_id, np.nan)
                    ),
                )
            )

    frequency_matrices, size_series, validation = _build_motif_library_from_candidates(
        candidates,
        flank_size=flank_size,
    )
    size_series.attrs[_MOTIF_LIBRARY_VALIDATION_ATTR] = validation
    return frequency_matrices, size_series


def build_motif_library_from_sequences(
    *,
    motif_sequences: Mapping[
        str, Sequence[str | Mapping[str, object] | ExplicitMotifSequence]
    ],
    flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
) -> tuple[dict[str, pd.DataFrame], pd.Series]:
    """Build motif frequency matrices from explicit per-kinase sequences.

    Each per-kinase entry can be either:
    - a bare sequence string (legacy/less-informative mode), or
    - structured metadata carrying `reference_id`, optional `site_id`,
      optional `kinase`, and `sequence`.
    """

    if flank_size < 0:
        raise ValueError("flank_size must be >= 0")

    candidates: list[_LibraryCandidate] = []
    for kinase, sequences in motif_sequences.items():
        kinase_name = str(kinase)
        for index, entry in enumerate(sequences):
            explicit_entry = _normalize_explicit_motif_sequence(
                entry=entry,
                default_kinase=kinase_name,
                default_reference_id=f"{kinase_name}#{index}",
            )
            candidates.append(
                _LibraryCandidate(
                    reference_id=explicit_entry.reference_id,
                    site_id=explicit_entry.site_id,
                    kinase=explicit_entry.kinase,
                    sequence_input=_normalize_sequence_value(explicit_entry.sequence),
                )
            )
    frequency_matrices, size_series, validation = _build_motif_library_from_candidates(
        candidates,
        flank_size=flank_size,
    )
    size_series.attrs[_MOTIF_LIBRARY_VALIDATION_ATTR] = validation
    return frequency_matrices, size_series


def get_motif_library_validation(
    motif_sizes: pd.Series,
) -> MotifLibraryValidationResult | None:
    """Extract attached motif-library validation diagnostics from motif sizes."""

    validation = motif_sizes.attrs.get(_MOTIF_LIBRARY_VALIDATION_ATTR)
    if isinstance(validation, MotifLibraryValidationResult):
        return validation
    return None


def _build_motif_library_from_candidates(
    candidates: Sequence[_LibraryCandidate],
    *,
    flank_size: int,
) -> tuple[dict[str, pd.DataFrame], pd.Series, MotifLibraryValidationResult]:
    expected_window_size = (2 * flank_size) + 1
    validator = MotifSequenceValidator(
        expected_window_size=expected_window_size,
        require_phospho_centre_residue=True,
        enforce_site_identity_format=True,
        window_length_policy="centred_superset",
        residue_validation_scope="centre_window",
    )
    validation_inputs = [
        SequenceValidationInput(
            site_id=candidate.reference_id,
            site_sequence=candidate.sequence_input,
            site_identity=candidate.site_id,
        )
        for candidate in candidates
    ]
    validation_result = validator.run(rows=validation_inputs)
    status_counts = Counter(row.status for row in validation_result.rows)

    accepted_windows_by_kinase: dict[str, list[str]] = {}
    validation_rows: list[MotifLibraryValidationRow] = []
    excluded_reference_ids: list[str] = []
    for candidate, row in zip(candidates, validation_result.rows, strict=True):
        validation_rows.append(
            MotifLibraryValidationRow(
                reference_id=candidate.reference_id,
                site_id=candidate.site_id,
                kinase=candidate.kinase,
                sequence=row.sequence,
                status=row.status,
                reason=row.reason,
                expected_centre_residue=row.expected_centre_residue,
                observed_centre_residue=row.observed_centre_residue,
                sequence_length=row.sequence_length,
            )
        )
        if row.status != SEQUENCE_VALIDATION_STATUS_VALID:
            excluded_reference_ids.append(candidate.reference_id)
            continue
        if row.sequence is None:
            continue
        accepted_window = _extract_sequence_window(row.sequence, flank_size)
        if _is_missing_scalar(accepted_window):
            excluded_reference_ids.append(candidate.reference_id)
            continue
        accepted_windows_by_kinase.setdefault(candidate.kinase, []).append(
            str(accepted_window)
        )

    frequency_matrices: dict[str, pd.DataFrame] = {}
    motif_sizes: dict[str, float] = {}
    for kinase, windows in accepted_windows_by_kinase.items():
        if not windows:
            continue
        windows_series = pd.Series(windows, dtype=object)
        frequency_matrices[kinase] = _build_frequency_matrix_from_windows(
            windows_series,
            context=f"kinase={kinase}",
        )
        motif_sizes[kinase] = float(len(windows_series))

    size_series = pd.Series(motif_sizes, dtype=float, name="motif_size")
    size_series.index.name = "kinase"
    validation = MotifLibraryValidationResult(
        total_reference_sequences=validation_result.total_sequences,
        accepted_reference_sequences=validation_result.valid_sequences,
        excluded_reference_sequences=validation_result.invalid_sequences,
        missing_sequences=status_counts[SEQUENCE_VALIDATION_STATUS_MISSING_SEQUENCE],
        short_sequences=status_counts[SEQUENCE_VALIDATION_STATUS_SHORT_SEQUENCE],
        off_centre_sequences=status_counts[
            SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE
        ],
        site_residue_mismatches=status_counts[
            SEQUENCE_VALIDATION_STATUS_SITE_RESIDUE_MISMATCH
        ],
        invalid_site_ids=status_counts[SEQUENCE_VALIDATION_STATUS_INVALID_SITE_ID],
        non_phospho_centre_residues=status_counts[
            SEQUENCE_VALIDATION_STATUS_NON_PHOSPHO_CENTRE_RESIDUE
        ],
        unsupported_residue_characters=status_counts[
            SEQUENCE_VALIDATION_STATUS_UNSUPPORTED_RESIDUE_CHARACTER
        ],
        sequences_excluded_from_motif_profile_construction=(
            validation_result.invalid_sequences
        ),
        expected_window_size=expected_window_size,
        supported_amino_acids=tuple(AMINO_ACIDS),
        accepted_window_length_policy=(
            f"input sequences must be centred and odd-length with minimum length "
            f"{expected_window_size}; scoring windows are centre-extracted to exactly "
            f"{expected_window_size} residues (2 * flank_size + 1)"
        ),
        unsupported_residue_policy=(
            "exclude any window containing non-canonical amino acids; "
            f"supported residues: {', '.join(AMINO_ACIDS)}"
        ),
        excluded_reference_ids=tuple(excluded_reference_ids),
        rows=tuple(validation_rows),
    )
    return frequency_matrices, size_series, validation


def _normalize_explicit_motif_sequence(
    *,
    entry: str | Mapping[str, object] | ExplicitMotifSequence,
    default_kinase: str,
    default_reference_id: str,
) -> ExplicitMotifSequence:
    if isinstance(entry, ExplicitMotifSequence):
        reference_id = _coerce_identifier(
            entry.reference_id, fallback=default_reference_id
        )
        kinase = _coerce_identifier(entry.kinase, fallback=default_kinase)
        if kinase != default_kinase:
            raise ValueError(
                "explicit motif sequence kinase metadata must match the parent "
                f"mapping key '{default_kinase}'"
            )
        return ExplicitMotifSequence(
            reference_id=reference_id,
            site_id=_coerce_optional_identifier(entry.site_id),
            kinase=kinase,
            sequence=entry.sequence,
        )
    if isinstance(entry, Mapping):
        reference_id = _coerce_identifier(
            entry.get("reference_id"),
            fallback=default_reference_id,
        )
        entry_kinase = _coerce_identifier(entry.get("kinase"), fallback=default_kinase)
        if entry_kinase != default_kinase:
            raise ValueError(
                "explicit motif sequence kinase metadata must match the parent "
                f"mapping key '{default_kinase}'"
            )
        return ExplicitMotifSequence(
            reference_id=reference_id,
            site_id=_coerce_optional_identifier(entry.get("site_id")),
            kinase=entry_kinase,
            sequence=entry.get("sequence"),
        )
    return ExplicitMotifSequence(
        reference_id=default_reference_id,
        site_id=None,
        kinase=default_kinase,
        sequence=entry,
    )


def _coerce_identifier(value: object, *, fallback: str) -> str:
    identifier = _coerce_optional_identifier(value)
    if identifier is None:
        return str(fallback)
    return identifier


def _coerce_optional_identifier(value: object) -> str | None:
    if value is None:
        return None
    try:
        if bool(pd.isna(value)):
            return None
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    if text == "":
        return None
    return text


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
        columns=kinases,
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


def _normalize_sequence_value(value: object) -> object:
    if value is None:
        return np.nan
    try:
        if bool(pd.isna(value)):
            return np.nan
    except (TypeError, ValueError):
        pass
    if not isinstance(value, str):
        return value
    sequence = value.strip().upper()
    if sequence == "":
        return np.nan
    return sequence


def _is_missing_scalar(value: object) -> bool:
    if value is None:
        return True
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _extract_sequence_window(value: object, flank_size: int | None) -> object:
    if _is_missing_scalar(value):
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
    windows = pd.Series(np.nan, index=raw_sequences.index.copy(), dtype=object)
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


def _require_fully_supported_encoded_sequences(
    encoded_sequences: np.ndarray,
    *,
    context: str,
) -> None:
    if (encoded_sequences < 0).any():
        raise ValueError(f"{context} contains unsupported residues")


__all__ = [
    "DEFAULT_MOTIF_FLANK_SIZE",
    "ExplicitMotifSequence",
    "MotifLibraryValidationResult",
    "MotifLibraryValidationRow",
    "MotifScoringResult",
    "SEQUENCE_SEMANTICS_CENTRED_SEQUENCE",
    "SEQUENCE_SEMANTICS_CENTRED_WINDOW",
    "build_motif_library",
    "build_motif_library_from_sequences",
    "get_motif_library_validation",
    "minmax_scale_columns",
    "score_phosphosite_motifs",
]
