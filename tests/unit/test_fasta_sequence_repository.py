from __future__ import annotations

from hashlib import sha256
from pathlib import Path

import pytest

from phospy.errors.input import PhosPyInputError
from phospy.sequences.repository import FastaProteinSequenceRepository


def test_uniprot_header_parsing(tmp_path: Path) -> None:
    fasta_path = tmp_path / "proteins.fasta"
    fasta_path.write_text(
        ">sp|P28482|MK01_HUMAN MAPK1_HUMAN\nMTEY\n",
        encoding="utf-8",
    )

    repository = FastaProteinSequenceRepository.from_path(fasta_path)
    lookup = repository.lookup("P28482")

    assert lookup.status == "found"
    assert lookup.record is not None
    assert lookup.record.accession == "P28482"
    assert lookup.record.source_record_id == "sp|P28482|MK01_HUMAN"


def test_simple_header_parsing(tmp_path: Path) -> None:
    fasta_path = tmp_path / "proteins.fasta"
    fasta_path.write_text(
        ">P99999 simple-protein\nACDE\n",
        encoding="utf-8",
    )

    repository = FastaProteinSequenceRepository.from_path(fasta_path)
    lookup = repository.lookup("P99999")

    assert lookup.status == "found"
    assert lookup.record is not None
    assert lookup.record.accession == "P99999"
    assert lookup.record.source_record_id == "P99999"


def test_multiline_sequence_concatenation(tmp_path: Path) -> None:
    fasta_path = tmp_path / "proteins.fasta"
    fasta_path.write_text(
        ">P1 protein\nac d\n\nE\tF\ng\n",
        encoding="utf-8",
    )

    repository = FastaProteinSequenceRepository.from_path(fasta_path)
    lookup = repository.lookup("P1")

    assert lookup.status == "found"
    assert lookup.record is not None
    assert lookup.record.sequence == "ACDEFG"


def test_file_digest_stability(tmp_path: Path) -> None:
    fasta_bytes = b">P1 protein\nACDE\n"
    fasta_path = tmp_path / "proteins.fasta"
    fasta_path.write_bytes(fasta_bytes)

    repository_a = FastaProteinSequenceRepository.from_path(fasta_path)
    repository_b = FastaProteinSequenceRepository.from_path(fasta_path)
    expected_digest = sha256(fasta_bytes).hexdigest()

    assert repository_a.metadata.sha256 == expected_digest
    assert repository_b.metadata.sha256 == expected_digest


def test_duplicate_accession_handling(tmp_path: Path) -> None:
    fasta_path = tmp_path / "proteins.fasta"
    fasta_path.write_text(
        ">sp|P11111|REC1\nAAAA\n>sp|P11111|REC2\nCCCC\n",
        encoding="utf-8",
    )

    repository = FastaProteinSequenceRepository.from_path(fasta_path)
    lookup = repository.lookup("P11111")

    assert repository.metadata.record_count == 2
    assert repository.metadata.duplicate_accessions == ("P11111",)
    assert lookup.status == "ambiguous_accession"
    assert lookup.record is None


def test_missing_accession_lookup(tmp_path: Path) -> None:
    fasta_path = tmp_path / "proteins.fasta"
    fasta_path.write_text(
        ">P1 protein\nACDE\n",
        encoding="utf-8",
    )
    repository = FastaProteinSequenceRepository.from_path(fasta_path)

    lookup = repository.lookup(None)

    assert lookup.status == "missing_accession"
    assert lookup.record is None
    assert lookup.reason is not None


def test_accession_not_found_lookup(tmp_path: Path) -> None:
    fasta_path = tmp_path / "proteins.fasta"
    fasta_path.write_text(
        ">P1 protein\nACDE\n",
        encoding="utf-8",
    )
    repository = FastaProteinSequenceRepository.from_path(fasta_path)

    lookup = repository.lookup("P40404")

    assert lookup.status == "accession_not_found"
    assert lookup.record is None


def test_invalid_sequence_character_rejection(tmp_path: Path) -> None:
    fasta_path = tmp_path / "proteins.fasta"
    fasta_path.write_text(
        ">P1 protein\nACD@\n",
        encoding="utf-8",
    )

    with pytest.raises(PhosPyInputError, match="invalid sequence character"):
        FastaProteinSequenceRepository.from_path(fasta_path)
