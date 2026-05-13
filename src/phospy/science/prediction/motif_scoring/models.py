from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

from phospy.science.prediction.sequence_validation import (
    SUPPORTED_AMINO_ACIDS,
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
            columns=pd.Index(columns),
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


def _normalize_sequence_value(value: object) -> object:
    if value is None:
        return np.nan
    if bool(pd.Series((value,), dtype="object").isna().iat[0]):
        return np.nan
    if not isinstance(value, str):
        return value
    sequence = value.strip().upper()
    if sequence == "":
        return np.nan
    return sequence


def _is_missing_scalar(value: object) -> bool:
    if value is None:
        return True
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


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
