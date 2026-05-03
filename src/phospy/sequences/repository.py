"""Local FASTA-backed protein sequence repository."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from phospy.sequences.fasta import (
    collect_duplicate_accessions,
    load_fasta_text_and_digest,
    parse_fasta_records,
)
from phospy.sequences.models import (
    FastaSourceMetadata,
    ProteinSequenceLookupResult,
    ProteinSequenceRecord,
)
from phospy.sequences.validation import (
    LOOKUP_STATUS_ACCESSION_NOT_FOUND,
    LOOKUP_STATUS_AMBIGUOUS_ACCESSION,
    LOOKUP_STATUS_FOUND,
    LOOKUP_STATUS_MISSING_ACCESSION,
    normalize_lookup_accession,
)


@dataclass(frozen=True, slots=True)
class FastaProteinSequenceRepository:
    """Strict local FASTA sequence repository keyed by accessions."""

    _records_by_accession: dict[str, ProteinSequenceRecord]
    _ambiguous_accessions: frozenset[str]
    metadata: FastaSourceMetadata

    @classmethod
    def from_path(
        cls,
        path: str | Path,
        *,
        source_label: str | None = None,
    ) -> FastaProteinSequenceRepository:
        normalized_path, text, digest = load_fasta_text_and_digest(path)
        records = parse_fasta_records(text)
        duplicates = collect_duplicate_accessions(records)
        duplicate_set = frozenset(duplicates)

        records_by_accession: dict[str, ProteinSequenceRecord] = {}
        for record in records:
            if record.accession in duplicate_set:
                continue
            records_by_accession[record.accession] = record

        label = source_label if source_label is not None else normalized_path.name
        if label.strip() == "":
            label = str(normalized_path)

        metadata = FastaSourceMetadata(
            source_label=label,
            source_path=str(normalized_path.resolve()),
            sha256=digest,
            record_count=len(records),
            duplicate_accessions=duplicates,
        )
        return cls(
            _records_by_accession=records_by_accession,
            _ambiguous_accessions=duplicate_set,
            metadata=metadata,
        )

    def lookup(self, accession: object) -> ProteinSequenceLookupResult:
        normalized_accession = normalize_lookup_accession(accession)
        if normalized_accession is None:
            return ProteinSequenceLookupResult(
                accession="" if accession is None else str(accession),
                status=LOOKUP_STATUS_MISSING_ACCESSION,
                record=None,
                reason="accession must be a non-empty string",
            )
        if normalized_accession in self._ambiguous_accessions:
            return ProteinSequenceLookupResult(
                accession=normalized_accession,
                status=LOOKUP_STATUS_AMBIGUOUS_ACCESSION,
                record=None,
                reason=(
                    "accession resolves to multiple FASTA records in this source; "
                    "lookup is ambiguous"
                ),
            )
        record = self._records_by_accession.get(normalized_accession)
        if record is None:
            return ProteinSequenceLookupResult(
                accession=normalized_accession,
                status=LOOKUP_STATUS_ACCESSION_NOT_FOUND,
                record=None,
                reason="accession not found in FASTA source",
            )
        return ProteinSequenceLookupResult(
            accession=normalized_accession,
            status=LOOKUP_STATUS_FOUND,
            record=record,
            reason=None,
        )
