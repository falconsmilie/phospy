"""Local FASTA parsing helpers for protein sequence repositories."""

from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
from pathlib import Path

from phospy.errors.input import PhosPyInputError
from phospy.sequences.models import ProteinSequenceRecord
from phospy.sequences.validation import ensure_local_path, validate_sequence_characters


def load_fasta_text_and_digest(path: str | Path) -> tuple[Path, str, str]:
    """Read local FASTA file bytes and return decoded text plus SHA-256 digest."""

    normalized_path = ensure_local_path(path)
    try:
        raw_bytes = normalized_path.read_bytes()
    except FileNotFoundError as exc:
        raise PhosPyInputError(f"FASTA file does not exist: {normalized_path}") from exc
    except PermissionError as exc:
        raise PhosPyInputError(
            f"permission denied while reading FASTA file: {normalized_path}"
        ) from exc
    except OSError as exc:
        raise PhosPyInputError(
            f"failed to read FASTA file '{normalized_path}': {exc}"
        ) from exc

    digest = sha256(raw_bytes).hexdigest()
    try:
        text = raw_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise PhosPyInputError(
            f"failed to decode FASTA file '{normalized_path}' as UTF-8: {exc}"
        ) from exc
    return normalized_path, text, digest


def parse_fasta_records(text: str) -> tuple[ProteinSequenceRecord, ...]:
    """Parse FASTA text into immutable protein sequence records."""

    records: list[ProteinSequenceRecord] = []
    current_header: str | None = None
    current_sequence_parts: list[str] = []
    current_header_line = 0

    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if line == "":
            continue
        if line.startswith(">"):
            if current_header is not None:
                records.append(
                    _build_record(
                        header=current_header,
                        sequence_parts=current_sequence_parts,
                        header_line=current_header_line,
                    )
                )
            current_header = line
            current_sequence_parts = []
            current_header_line = line_number
            continue

        if current_header is None:
            raise PhosPyInputError(
                "invalid FASTA input: encountered sequence content before first header "
                f"at line {line_number}"
            )
        sanitized = "".join(raw_line.split()).upper()
        if sanitized:
            accession = _extract_accession_from_header(current_header)
            validate_sequence_characters(
                sanitized,
                accession=accession,
                line_number=line_number,
            )
            current_sequence_parts.append(sanitized)

    if current_header is not None:
        records.append(
            _build_record(
                header=current_header,
                sequence_parts=current_sequence_parts,
                header_line=current_header_line,
            )
        )
    return tuple(records)


def collect_duplicate_accessions(
    records: tuple[ProteinSequenceRecord, ...],
) -> tuple[str, ...]:
    """Collect sorted duplicate accession identifiers from parsed records."""

    counts: dict[str, int] = defaultdict(int)
    for record in records:
        counts[record.accession] += 1
    duplicates = sorted(accession for accession, count in counts.items() if count > 1)
    return tuple(duplicates)


def _build_record(
    *,
    header: str,
    sequence_parts: list[str],
    header_line: int,
) -> ProteinSequenceRecord:
    header_content = header[1:].strip()
    if header_content == "":
        raise PhosPyInputError(
            f"invalid FASTA header at line {header_line}: header is empty"
        )

    source_record_id = header_content.split(maxsplit=1)[0]
    accession = _extract_accession_from_header(header)
    sequence = "".join(sequence_parts).upper()
    validate_sequence_characters(sequence, accession=accession, line_number=header_line)
    return ProteinSequenceRecord(
        accession=accession,
        sequence=sequence,
        description=header_content,
        source_record_id=source_record_id,
    )


def _extract_accession_from_header(header: str) -> str:
    header_content = header[1:].strip()
    if header_content == "":
        raise PhosPyInputError(
            "invalid FASTA header: missing record identifier after '>'"
        )
    header_token = header_content.split(maxsplit=1)[0]
    if header_token.count("|") >= 2:
        first = header_token.find("|")
        second = header_token.find("|", first + 1)
        accession = header_token[first + 1 : second].strip()
        if accession == "":
            raise PhosPyInputError(
                f"invalid UniProt-style FASTA header '{header_content}': accession is empty"
            )
        return accession
    accession = header_token.strip()
    if accession == "":
        raise PhosPyInputError(
            f"invalid FASTA header '{header_content}': accession token is empty"
        )
    return accession
