"""Validation helpers for local FASTA sequence repositories."""

from __future__ import annotations

from pathlib import Path

from phospy.errors.input import PhosPyInputError

LOOKUP_STATUS_FOUND = "found"
LOOKUP_STATUS_MISSING_ACCESSION = "missing_accession"
LOOKUP_STATUS_ACCESSION_NOT_FOUND = "accession_not_found"
LOOKUP_STATUS_AMBIGUOUS_ACCESSION = "ambiguous_accession"

_ALLOWED_SEQUENCE_CHARACTERS = frozenset("ACDEFGHIKLMNPQRSTVWYXBZUOJ*")


def ensure_local_path(path: str | Path) -> Path:
    """Normalize and validate that an input points to a local filesystem path."""

    if isinstance(path, str) and "://" in path:
        raise PhosPyInputError(
            "local FASTA repository only supports filesystem paths; "
            "remote URLs are not supported"
        )
    return Path(path)


def normalize_lookup_accession(accession: object) -> str | None:
    """Normalize lookup accession input for strict repository matching."""

    if not isinstance(accession, str):
        return None
    normalized = accession.strip()
    if normalized == "":
        return None
    return normalized


def validate_sequence_characters(
    sequence: str,
    *,
    accession: str,
    line_number: int,
) -> None:
    """Validate that a parsed FASTA sequence uses supported residue symbols."""

    invalid = sorted(
        {
            character
            for character in sequence
            if character not in _ALLOWED_SEQUENCE_CHARACTERS
        }
    )
    if invalid:
        joined = ", ".join(invalid)
        raise PhosPyInputError(
            "invalid sequence character(s) in FASTA record "
            f"'{accession}' at line {line_number}: {joined}; "
            "allowed characters are uppercase amino-acid letters plus X, B, Z, U, O, J, and *"
        )
