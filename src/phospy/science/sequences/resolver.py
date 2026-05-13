"""Phosphosite-to-sequence resolution over local FASTA repositories."""

from __future__ import annotations

import re
from dataclasses import dataclass

from phospy.errors.input import PhosPyInputError
from phospy.science.sequences.repository import FastaProteinSequenceRepository
from phospy.science.sequences.validation import (
    LOOKUP_STATUS_ACCESSION_NOT_FOUND,
    LOOKUP_STATUS_AMBIGUOUS_ACCESSION,
    LOOKUP_STATUS_MISSING_ACCESSION,
    normalize_lookup_accession,
)

RESOLUTION_STATUS_RESOLVED = "resolved"
RESOLUTION_STATUS_MISSING_ACCESSION = "missing_accession"
RESOLUTION_STATUS_ACCESSION_NOT_FOUND = "accession_not_found"
RESOLUTION_STATUS_INVALID_SITE_TOKEN = "invalid_site_token"
RESOLUTION_STATUS_SITE_OUT_OF_BOUNDS = "site_out_of_bounds"
RESOLUTION_STATUS_RESIDUE_MISMATCH = "residue_mismatch"
RESOLUTION_STATUS_INSUFFICIENT_FLANKING_SEQUENCE = "insufficient_flanking_sequence"
RESOLUTION_STATUS_AMBIGUOUS_ACCESSION = "ambiguous_accession"

_SITE_TOKEN_PATTERN = re.compile(r"^\s*([A-Za-z])([1-9][0-9]*)\s*$")


@dataclass(frozen=True, slots=True)
class PhosphositeSequenceResolutionRequest:
    accession: object
    site_token: object
    flank_size: int


@dataclass(frozen=True, slots=True)
class PhosphositeSequenceResolutionResult:
    accession: str | None
    site_token: str | None
    status: str
    site_sequence: str | None
    reason: str | None
    flank_size: int
    protein_length: int | None
    site_position: int | None
    expected_residue: str | None
    observed_residue: str | None


class PhosphositeSequenceResolver:
    """Resolve one phosphosite token to a centred sequence window."""

    def run(
        self,
        request: PhosphositeSequenceResolutionRequest,
        repository: FastaProteinSequenceRepository,
    ) -> PhosphositeSequenceResolutionResult:
        flank_size = _validate_flank_size(request.flank_size)
        normalized_accession = normalize_lookup_accession(request.accession)
        if normalized_accession is None:
            return PhosphositeSequenceResolutionResult(
                accession=None,
                site_token=_normalize_site_token_text(request.site_token),
                status=RESOLUTION_STATUS_MISSING_ACCESSION,
                site_sequence=None,
                reason="accession must be a non-empty string",
                flank_size=flank_size,
                protein_length=None,
                site_position=None,
                expected_residue=None,
                observed_residue=None,
            )

        lookup = repository.lookup(normalized_accession)
        if lookup.status in {
            LOOKUP_STATUS_MISSING_ACCESSION,
            RESOLUTION_STATUS_MISSING_ACCESSION,
        }:
            return PhosphositeSequenceResolutionResult(
                accession=None,
                site_token=_normalize_site_token_text(request.site_token),
                status=RESOLUTION_STATUS_MISSING_ACCESSION,
                site_sequence=None,
                reason=lookup.reason or "accession must be a non-empty string",
                flank_size=flank_size,
                protein_length=None,
                site_position=None,
                expected_residue=None,
                observed_residue=None,
            )
        if lookup.status in {
            LOOKUP_STATUS_ACCESSION_NOT_FOUND,
            RESOLUTION_STATUS_ACCESSION_NOT_FOUND,
        }:
            return PhosphositeSequenceResolutionResult(
                accession=normalized_accession,
                site_token=_normalize_site_token_text(request.site_token),
                status=RESOLUTION_STATUS_ACCESSION_NOT_FOUND,
                site_sequence=None,
                reason=lookup.reason or "accession not found in FASTA source",
                flank_size=flank_size,
                protein_length=None,
                site_position=None,
                expected_residue=None,
                observed_residue=None,
            )
        if lookup.status in {
            LOOKUP_STATUS_AMBIGUOUS_ACCESSION,
            RESOLUTION_STATUS_AMBIGUOUS_ACCESSION,
        }:
            return PhosphositeSequenceResolutionResult(
                accession=normalized_accession,
                site_token=_normalize_site_token_text(request.site_token),
                status=RESOLUTION_STATUS_AMBIGUOUS_ACCESSION,
                site_sequence=None,
                reason=lookup.reason or "accession resolves to multiple FASTA records",
                flank_size=flank_size,
                protein_length=None,
                site_position=None,
                expected_residue=None,
                observed_residue=None,
            )

        token_parsed = _parse_site_token(request.site_token)
        if token_parsed is None:
            return PhosphositeSequenceResolutionResult(
                accession=normalized_accession,
                site_token=None,
                status=RESOLUTION_STATUS_INVALID_SITE_TOKEN,
                site_sequence=None,
                reason=(
                    "site token must follow '<residue><position>' with one residue letter "
                    "and a positive integer position (for example 'S123')"
                ),
                flank_size=flank_size,
                protein_length=None,
                site_position=None,
                expected_residue=None,
                observed_residue=None,
            )

        site_token, expected_residue, site_position = token_parsed
        record = lookup.record
        if record is None:
            return PhosphositeSequenceResolutionResult(
                accession=normalized_accession,
                site_token=site_token,
                status=RESOLUTION_STATUS_ACCESSION_NOT_FOUND,
                site_sequence=None,
                reason="accession lookup did not return a sequence record",
                flank_size=flank_size,
                protein_length=None,
                site_position=site_position,
                expected_residue=expected_residue,
                observed_residue=None,
            )

        sequence = record.sequence
        protein_length = len(sequence)
        if site_position > protein_length:
            return PhosphositeSequenceResolutionResult(
                accession=normalized_accession,
                site_token=site_token,
                status=RESOLUTION_STATUS_SITE_OUT_OF_BOUNDS,
                site_sequence=None,
                reason=(
                    f"site position {site_position} exceeds protein length "
                    f"{protein_length}"
                ),
                flank_size=flank_size,
                protein_length=protein_length,
                site_position=site_position,
                expected_residue=expected_residue,
                observed_residue=None,
            )

        observed_residue = sequence[site_position - 1]
        if observed_residue != expected_residue:
            return PhosphositeSequenceResolutionResult(
                accession=normalized_accession,
                site_token=site_token,
                status=RESOLUTION_STATUS_RESIDUE_MISMATCH,
                site_sequence=None,
                reason=(
                    f"expected residue '{expected_residue}' at position {site_position}, "
                    f"observed '{observed_residue}'"
                ),
                flank_size=flank_size,
                protein_length=protein_length,
                site_position=site_position,
                expected_residue=expected_residue,
                observed_residue=observed_residue,
            )

        n_terminal_position = site_position - flank_size
        c_terminal_position = site_position + flank_size
        if n_terminal_position < 1 or c_terminal_position > protein_length:
            return PhosphositeSequenceResolutionResult(
                accession=normalized_accession,
                site_token=site_token,
                status=RESOLUTION_STATUS_INSUFFICIENT_FLANKING_SEQUENCE,
                site_sequence=None,
                reason=(
                    f"site position {site_position} with flank size {flank_size} "
                    "does not have sufficient flanking sequence"
                ),
                flank_size=flank_size,
                protein_length=protein_length,
                site_position=site_position,
                expected_residue=expected_residue,
                observed_residue=observed_residue,
            )

        start_index = n_terminal_position - 1
        stop_index = c_terminal_position
        site_sequence = sequence[start_index:stop_index]
        return PhosphositeSequenceResolutionResult(
            accession=normalized_accession,
            site_token=site_token,
            status=RESOLUTION_STATUS_RESOLVED,
            site_sequence=site_sequence,
            reason=None,
            flank_size=flank_size,
            protein_length=protein_length,
            site_position=site_position,
            expected_residue=expected_residue,
            observed_residue=observed_residue,
        )


def _validate_flank_size(flank_size: int) -> int:
    if isinstance(flank_size, bool) or not isinstance(flank_size, int):
        raise PhosPyInputError(
            "flank_size must be an integer greater than or equal to zero"
        )
    if flank_size < 0:
        raise PhosPyInputError(
            "flank_size must be an integer greater than or equal to zero"
        )
    return flank_size


def _normalize_site_token_text(site_token: object) -> str | None:
    if not isinstance(site_token, str):
        return None
    normalized = site_token.strip()
    if normalized == "":
        return None
    return normalized


def _parse_site_token(site_token: object) -> tuple[str, str, int] | None:
    if not isinstance(site_token, str):
        return None
    match = _SITE_TOKEN_PATTERN.fullmatch(site_token)
    if match is None:
        return None
    residue = match.group(1).upper()
    position = int(match.group(2))
    return f"{residue}{position}", residue, position


__all__ = [
    "PhosphositeSequenceResolutionRequest",
    "PhosphositeSequenceResolutionResult",
    "PhosphositeSequenceResolver",
    "RESOLUTION_STATUS_ACCESSION_NOT_FOUND",
    "RESOLUTION_STATUS_AMBIGUOUS_ACCESSION",
    "RESOLUTION_STATUS_INSUFFICIENT_FLANKING_SEQUENCE",
    "RESOLUTION_STATUS_INVALID_SITE_TOKEN",
    "RESOLUTION_STATUS_MISSING_ACCESSION",
    "RESOLUTION_STATUS_RESIDUE_MISMATCH",
    "RESOLUTION_STATUS_RESOLVED",
    "RESOLUTION_STATUS_SITE_OUT_OF_BOUNDS",
]
