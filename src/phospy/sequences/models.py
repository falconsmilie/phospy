"""Sequence-domain models for local protein FASTA repositories."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class ProteinSequenceRecord:
    """Single parsed FASTA protein sequence record."""

    accession: str
    sequence: str
    description: str
    source_record_id: str


@dataclass(frozen=True, slots=True)
class FastaSourceMetadata:
    """Metadata captured from one parsed FASTA source file."""

    source_label: str
    source_path: str
    sha256: str
    record_count: int
    duplicate_accessions: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ProteinSequenceLookupResult:
    """Repository lookup outcome for one requested accession."""

    accession: str
    status: str
    record: ProteinSequenceRecord | None
    reason: str | None
