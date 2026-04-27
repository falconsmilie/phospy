"""Sequence-context validation for motif scoring."""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

import pandas as pd

SUPPORTED_AMINO_ACIDS: tuple[str, ...] = (
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
_SUPPORTED_AMINO_ACID_SET = frozenset(SUPPORTED_AMINO_ACIDS)
_SITE_IDENTITY_PATTERN = re.compile(r"^\s*[^;]+?\s*;\s*(?P<site>[^;]+?)\s*;\s*$")
PHOSPHO_COMPATIBLE_RESIDUES: tuple[str, ...] = ("S", "T", "Y")
_PHOSPHO_COMPATIBLE_RESIDUE_SET = frozenset(PHOSPHO_COMPATIBLE_RESIDUES)

SEQUENCE_VALIDATION_STATUS_VALID = "valid"
SEQUENCE_VALIDATION_STATUS_MISSING_SEQUENCE = "missing_sequence"
SEQUENCE_VALIDATION_STATUS_SHORT_SEQUENCE = "short_sequence"
SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE = "off_centre_sequence"
SEQUENCE_VALIDATION_STATUS_SITE_RESIDUE_MISMATCH = "site_residue_mismatch"
SEQUENCE_VALIDATION_STATUS_INVALID_SITE_ID = "invalid_site_id"
SEQUENCE_VALIDATION_STATUS_NON_PHOSPHO_CENTRE_RESIDUE = "non_phospho_centre_residue"
SEQUENCE_VALIDATION_STATUS_UNSUPPORTED_RESIDUE_CHARACTER = (
    "unsupported_residue_character"
)

SequenceValidationStatus = Literal[
    "valid",
    "missing_sequence",
    "short_sequence",
    "off_centre_sequence",
    "site_residue_mismatch",
    "invalid_site_id",
    "non_phospho_centre_residue",
    "unsupported_residue_character",
]


@dataclass(frozen=True, slots=True)
class SequenceValidationInput:
    """Input row for motif sequence-context validation."""

    site_id: str
    site_sequence: object
    site_identity: str | None = None


@dataclass(frozen=True, slots=True)
class SequenceValidationRow:
    """Per-site validation outcome."""

    site_id: str
    sequence: str | None
    status: SequenceValidationStatus
    reason: str | None
    expected_centre_residue: str | None
    observed_centre_residue: str | None
    sequence_length: int | None


@dataclass(frozen=True, slots=True)
class SequenceValidationResult:
    """Structured sequence-validation summary and per-row outcomes."""

    total_sequences: int
    valid_sequences: int
    invalid_sequences: int
    short_sequences: int
    off_centre_sequences: int
    site_residue_mismatches: int
    invalid_site_ids: int
    non_phospho_centre_residues: int
    unsupported_residue_characters: int
    sequences_excluded_from_motif_scoring: int
    excluded_site_ids: tuple[str, ...]
    rows: tuple[SequenceValidationRow, ...]

    def summary(self) -> dict[str, int]:
        """Return workflow-level summary diagnostics."""

        return {
            "total_sequences": self.total_sequences,
            "valid_sequences": self.valid_sequences,
            "invalid_sequences": self.invalid_sequences,
            "short_sequences": self.short_sequences,
            "off_centre_sequences": self.off_centre_sequences,
            "site_residue_mismatches": self.site_residue_mismatches,
            "invalid_site_ids": self.invalid_site_ids,
            "non_phospho_centre_residues": self.non_phospho_centre_residues,
            "unsupported_residue_characters": self.unsupported_residue_characters,
            "sequences_excluded_from_motif_scoring": (
                self.sequences_excluded_from_motif_scoring
            ),
        }


class MotifSequenceValidator:
    """Validate motif-scoring sequence windows for one scoring run."""

    def __init__(
        self,
        *,
        expected_window_size: int,
        supported_amino_acids: Sequence[str] = SUPPORTED_AMINO_ACIDS,
        require_phospho_centre_residue: bool = False,
        enforce_site_identity_format: bool = False,
    ) -> None:
        if expected_window_size <= 0:
            raise ValueError("expected_window_size must be > 0")
        self.expected_window_size = int(expected_window_size)
        self.expected_centre_index = self.expected_window_size // 2
        self._supported_amino_acids = frozenset(supported_amino_acids)
        self._require_phospho_centre_residue = bool(require_phospho_centre_residue)
        self._enforce_site_identity_format = bool(enforce_site_identity_format)

    def run(
        self, *, rows: Sequence[SequenceValidationInput]
    ) -> SequenceValidationResult:
        validation_rows: list[SequenceValidationRow] = []
        excluded_site_ids: list[str] = []
        status_counts: Counter[str] = Counter()

        for row in rows:
            validation_row = self._validate_row(row)
            validation_rows.append(validation_row)
            status_counts[validation_row.status] += 1
            if validation_row.status != SEQUENCE_VALIDATION_STATUS_VALID:
                excluded_site_ids.append(validation_row.site_id)

        total_sequences = len(validation_rows)
        valid_sequences = status_counts[SEQUENCE_VALIDATION_STATUS_VALID]
        invalid_sequences = total_sequences - valid_sequences
        return SequenceValidationResult(
            total_sequences=total_sequences,
            valid_sequences=valid_sequences,
            invalid_sequences=invalid_sequences,
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
            sequences_excluded_from_motif_scoring=invalid_sequences,
            excluded_site_ids=tuple(excluded_site_ids),
            rows=tuple(validation_rows),
        )

    def _validate_row(self, row: SequenceValidationInput) -> SequenceValidationRow:
        site_id = str(row.site_id)
        sequence = _coerce_sequence(row.site_sequence)
        site_identity = _coerce_site_identity(row.site_identity)
        site_identity_provided = row.site_identity is not None
        expected_residue = _parse_expected_site_residue(
            site_identity if site_identity is not None else site_id
        )
        if site_identity_provided and self._enforce_site_identity_format:
            if site_identity is None or expected_residue is None:
                return SequenceValidationRow(
                    site_id=site_id,
                    sequence=sequence,
                    status=SEQUENCE_VALIDATION_STATUS_INVALID_SITE_ID,
                    reason=(
                        "site_id must follow the '<protein>;<residue><position>;' "
                        "shape (for example 'MAPK1;S202;')"
                    ),
                    expected_centre_residue=None,
                    observed_centre_residue=(
                        None
                        if sequence is None
                        else sequence[self.expected_centre_index]
                        if len(sequence) > self.expected_centre_index
                        else None
                    ),
                    sequence_length=None if sequence is None else len(sequence),
                )
        if sequence is None:
            return SequenceValidationRow(
                site_id=site_id,
                sequence=None,
                status=SEQUENCE_VALIDATION_STATUS_MISSING_SEQUENCE,
                reason="site_sequence is missing or blank",
                expected_centre_residue=expected_residue,
                observed_centre_residue=None,
                sequence_length=None,
            )

        sequence_length = len(sequence)
        if sequence_length < self.expected_window_size:
            return SequenceValidationRow(
                site_id=site_id,
                sequence=sequence,
                status=SEQUENCE_VALIDATION_STATUS_SHORT_SEQUENCE,
                reason=(
                    f"sequence length {sequence_length} is shorter than expected "
                    f"window length {self.expected_window_size}"
                ),
                expected_centre_residue=expected_residue,
                observed_centre_residue=None,
                sequence_length=sequence_length,
            )

        unsupported_characters = sorted(
            {
                character
                for character in sequence
                if character not in self._supported_amino_acids
            }
        )
        if unsupported_characters:
            joined = ", ".join(unsupported_characters)
            return SequenceValidationRow(
                site_id=site_id,
                sequence=sequence,
                status=SEQUENCE_VALIDATION_STATUS_UNSUPPORTED_RESIDUE_CHARACTER,
                reason=f"unsupported residue character(s): {joined}",
                expected_centre_residue=expected_residue,
                observed_centre_residue=sequence[self.expected_centre_index],
                sequence_length=sequence_length,
            )

        observed_residue = sequence[self.expected_centre_index]
        if (
            self._require_phospho_centre_residue
            and observed_residue not in _PHOSPHO_COMPATIBLE_RESIDUE_SET
        ):
            return SequenceValidationRow(
                site_id=site_id,
                sequence=sequence,
                status=SEQUENCE_VALIDATION_STATUS_NON_PHOSPHO_CENTRE_RESIDUE,
                reason=(
                    f"centre residue '{observed_residue}' is not phospho-compatible; "
                    f"expected one of {', '.join(PHOSPHO_COMPATIBLE_RESIDUES)}"
                ),
                expected_centre_residue=expected_residue,
                observed_centre_residue=observed_residue,
                sequence_length=sequence_length,
            )
        if expected_residue is not None and observed_residue != expected_residue:
            if expected_residue in sequence:
                status = SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE
            else:
                status = SEQUENCE_VALIDATION_STATUS_SITE_RESIDUE_MISMATCH
            return SequenceValidationRow(
                site_id=site_id,
                sequence=sequence,
                status=status,
                reason=(
                    f"expected centre residue '{expected_residue}' at index "
                    f"{self.expected_centre_index}; observed '{observed_residue}'"
                ),
                expected_centre_residue=expected_residue,
                observed_centre_residue=observed_residue,
                sequence_length=sequence_length,
            )

        return SequenceValidationRow(
            site_id=site_id,
            sequence=sequence,
            status=SEQUENCE_VALIDATION_STATUS_VALID,
            reason=None,
            expected_centre_residue=expected_residue,
            observed_centre_residue=observed_residue,
            sequence_length=sequence_length,
        )


def _coerce_sequence(value: object) -> str | None:
    if value is None:
        return None
    try:
        if bool(pd.isna(value)):
            return None
    except (TypeError, ValueError):
        pass
    if not isinstance(value, str):
        return None
    sequence = value.strip().upper()
    if sequence == "":
        return None
    return sequence


def _coerce_site_identity(value: object) -> str | None:
    if value is None:
        return None
    try:
        if bool(pd.isna(value)):
            return None
    except (TypeError, ValueError):
        pass
    site_identity = str(value).strip()
    if site_identity == "":
        return None
    return site_identity


def _parse_expected_site_residue(site_id: str) -> str | None:
    match = _SITE_IDENTITY_PATTERN.fullmatch(site_id)
    if match is None:
        return None
    site_token = match.group("site").strip().upper()
    if site_token == "":
        return None
    residue = site_token[0]
    if residue not in _SUPPORTED_AMINO_ACID_SET:
        return None
    return residue


__all__ = [
    "MotifSequenceValidator",
    "PHOSPHO_COMPATIBLE_RESIDUES",
    "SEQUENCE_VALIDATION_STATUS_INVALID_SITE_ID",
    "SEQUENCE_VALIDATION_STATUS_MISSING_SEQUENCE",
    "SEQUENCE_VALIDATION_STATUS_NON_PHOSPHO_CENTRE_RESIDUE",
    "SEQUENCE_VALIDATION_STATUS_OFF_CENTRE_SEQUENCE",
    "SEQUENCE_VALIDATION_STATUS_SHORT_SEQUENCE",
    "SEQUENCE_VALIDATION_STATUS_SITE_RESIDUE_MISMATCH",
    "SEQUENCE_VALIDATION_STATUS_UNSUPPORTED_RESIDUE_CHARACTER",
    "SEQUENCE_VALIDATION_STATUS_VALID",
    "SUPPORTED_AMINO_ACIDS",
    "SequenceValidationInput",
    "SequenceValidationResult",
    "SequenceValidationRow",
    "SequenceValidationStatus",
]
