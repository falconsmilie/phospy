from __future__ import annotations

from pathlib import Path

import pytest

from phospy.errors.input import PhosPyInputError
from phospy.science.sequences.models import FastaSourceMetadata, ProteinSequenceRecord
from phospy.science.sequences.repository import FastaProteinSequenceRepository
from phospy.science.sequences.resolver import (
    RESOLUTION_STATUS_ACCESSION_NOT_FOUND,
    RESOLUTION_STATUS_AMBIGUOUS_ACCESSION,
    RESOLUTION_STATUS_INSUFFICIENT_FLANKING_SEQUENCE,
    RESOLUTION_STATUS_INVALID_SITE_TOKEN,
    RESOLUTION_STATUS_MISSING_ACCESSION,
    RESOLUTION_STATUS_RESIDUE_MISMATCH,
    RESOLUTION_STATUS_RESOLVED,
    RESOLUTION_STATUS_SITE_OUT_OF_BOUNDS,
    PhosphositeSequenceResolutionRequest,
    PhosphositeSequenceResolver,
)


def test_successful_resolution_for_p28482_t185(tmp_path: Path) -> None:
    sequence = ("A" * 184) + "T" + ("C" * 40)
    repository = _build_repository(
        tmp_path=tmp_path,
        fasta_text=f">sp|P28482|MAPK1_HUMAN\n{sequence}\n",
    )

    result = PhosphositeSequenceResolver().run(
        PhosphositeSequenceResolutionRequest(
            accession="P28482",
            site_token="T185",
            flank_size=7,
        ),
        repository,
    )

    assert result.status == RESOLUTION_STATUS_RESOLVED
    assert result.reason is None
    assert result.site_sequence == ("A" * 7) + "T" + ("C" * 7)
    assert result.protein_length == len(sequence)
    assert result.site_position == 185
    assert result.expected_residue == "T"
    assert result.observed_residue == "T"


def test_missing_accession_returns_status(tmp_path: Path) -> None:
    repository = _build_repository(tmp_path=tmp_path, fasta_text=">P1\nACDE\n")

    result = PhosphositeSequenceResolver().run(
        PhosphositeSequenceResolutionRequest(
            accession=" ",
            site_token="S2",
            flank_size=1,
        ),
        repository,
    )

    assert result.status == RESOLUTION_STATUS_MISSING_ACCESSION
    assert result.reason


def test_accession_not_found_returns_status(tmp_path: Path) -> None:
    repository = _build_repository(tmp_path=tmp_path, fasta_text=">P1\nACDE\n")

    result = PhosphositeSequenceResolver().run(
        PhosphositeSequenceResolutionRequest(
            accession="P40404",
            site_token="S2",
            flank_size=1,
        ),
        repository,
    )

    assert result.status == RESOLUTION_STATUS_ACCESSION_NOT_FOUND
    assert result.reason


def test_ambiguous_accession_returns_status(tmp_path: Path) -> None:
    repository = _build_repository(
        tmp_path=tmp_path,
        fasta_text=">sp|P11111|REC1\nAAAA\n>sp|P11111|REC2\nCCCC\n",
    )

    result = PhosphositeSequenceResolver().run(
        PhosphositeSequenceResolutionRequest(
            accession="P11111",
            site_token="S2",
            flank_size=1,
        ),
        repository,
    )

    assert result.status == RESOLUTION_STATUS_AMBIGUOUS_ACCESSION
    assert result.reason


def test_invalid_site_token_returns_status(tmp_path: Path) -> None:
    repository = _build_repository(tmp_path=tmp_path, fasta_text=">P1\nASAA\n")

    result = PhosphositeSequenceResolver().run(
        PhosphositeSequenceResolutionRequest(
            accession="P1",
            site_token="S0",
            flank_size=1,
        ),
        repository,
    )

    assert result.status == RESOLUTION_STATUS_INVALID_SITE_TOKEN
    assert result.reason


def test_out_of_bounds_site_returns_status(tmp_path: Path) -> None:
    repository = _build_repository(tmp_path=tmp_path, fasta_text=">P1\nASAA\n")

    result = PhosphositeSequenceResolver().run(
        PhosphositeSequenceResolutionRequest(
            accession="P1",
            site_token="S99",
            flank_size=1,
        ),
        repository,
    )

    assert result.status == RESOLUTION_STATUS_SITE_OUT_OF_BOUNDS
    assert result.reason


def test_residue_mismatch_returns_status(tmp_path: Path) -> None:
    repository = _build_repository(tmp_path=tmp_path, fasta_text=">P1\nAAAAA\n")

    result = PhosphositeSequenceResolver().run(
        PhosphositeSequenceResolutionRequest(
            accession="P1",
            site_token="S3",
            flank_size=1,
        ),
        repository,
    )

    assert result.status == RESOLUTION_STATUS_RESIDUE_MISMATCH
    assert result.expected_residue == "S"
    assert result.observed_residue == "A"
    assert result.reason


def test_insufficient_n_terminal_flank_returns_status(tmp_path: Path) -> None:
    repository = _build_repository(tmp_path=tmp_path, fasta_text=">P1\nASAAA\n")

    result = PhosphositeSequenceResolver().run(
        PhosphositeSequenceResolutionRequest(
            accession="P1",
            site_token="S2",
            flank_size=2,
        ),
        repository,
    )

    assert result.status == RESOLUTION_STATUS_INSUFFICIENT_FLANKING_SEQUENCE
    assert result.reason


def test_insufficient_c_terminal_flank_returns_status(tmp_path: Path) -> None:
    repository = _build_repository(tmp_path=tmp_path, fasta_text=">P1\nAAASA\n")

    result = PhosphositeSequenceResolver().run(
        PhosphositeSequenceResolutionRequest(
            accession="P1",
            site_token="S4",
            flank_size=2,
        ),
        repository,
    )

    assert result.status == RESOLUTION_STATUS_INSUFFICIENT_FLANKING_SEQUENCE
    assert result.reason


def test_central_residue_correctness() -> None:
    sequence = "QQQSQQQ"
    result = PhosphositeSequenceResolver().run(
        PhosphositeSequenceResolutionRequest(
            accession="P1",
            site_token="s4",
            flank_size=3,
        ),
        _repository_from_sequence(sequence),
    )

    assert result.status == RESOLUTION_STATUS_RESOLVED
    assert result.site_sequence is not None
    centre_index = result.flank_size
    assert result.site_sequence[centre_index] == "S"
    assert result.site_sequence[centre_index] == result.expected_residue


def test_output_sequence_length() -> None:
    result = PhosphositeSequenceResolver().run(
        PhosphositeSequenceResolutionRequest(
            accession="P1",
            site_token="T5",
            flank_size=2,
        ),
        _repository_from_sequence("AAAATAAAA"),
    )

    assert result.status == RESOLUTION_STATUS_RESOLVED
    assert result.site_sequence is not None
    assert len(result.site_sequence) == (2 * result.flank_size) + 1


def test_negative_flank_size_raises_input_error() -> None:
    with pytest.raises(PhosPyInputError, match="flank_size"):
        PhosphositeSequenceResolver().run(
            PhosphositeSequenceResolutionRequest(
                accession="P1",
                site_token="S1",
                flank_size=-1,
            ),
            _repository_from_sequence("SAAAA"),
        )


def _build_repository(
    *,
    tmp_path: Path,
    fasta_text: str,
) -> FastaProteinSequenceRepository:
    fasta_path = tmp_path / "proteins.fasta"
    fasta_path.write_text(fasta_text, encoding="utf-8")
    return FastaProteinSequenceRepository.from_path(fasta_path)


def _repository_from_sequence(sequence: str) -> FastaProteinSequenceRepository:
    return FastaProteinSequenceRepository(
        _records_by_accession={
            "P1": ProteinSequenceRecord(
                accession="P1",
                sequence=sequence,
                description="P1",
                source_record_id="P1",
            )
        },
        _ambiguous_accessions=frozenset(),
        metadata=FastaSourceMetadata(
            source_label="in-memory",
            source_path="in-memory",
            sha256="",
            record_count=1,
            duplicate_accessions=(),
        ),
    )
