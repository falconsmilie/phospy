from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy.errors.validation import ReferenceValidationError
from phospy.io.bundles.reference_sources import ReferenceSourceTableReader
from phospy.science.references import (
    KinaseLibraryResidueClass,
    KinaseLibraryResourceLoader,
)


def test_kinase_library_loader_loads_ser_thr_matrix(tmp_path: Path) -> None:
    path = _write_kinase_library_fixture(
        tmp_path, kinase="AKT1", residue_class="ser_thr"
    )

    resource = KinaseLibraryResourceLoader(
        source_reader=ReferenceSourceTableReader()
    ).run(path)

    assert resource.source_version == "synthetic-kl-v1"
    assert resource.score_scale == "raw_log2_enrichment"
    assert resource.organisms == ("human", "mouse")
    assert resource.sequence_window.upstream_residues == 1
    assert resource.sequence_window.downstream_residues == 1
    matrix = resource.matrix_for("AKT1", KinaseLibraryResidueClass.SER_THR)
    assert matrix.kinase_family == "AGC"
    assert matrix.kinase_group == "AGC"
    assert matrix.score_table.loc["A", -1] == pytest.approx(12.5)
    assert matrix.score_table.loc["S", 0] == pytest.approx(3.25)


def test_kinase_library_loader_loads_tyr_matrix(tmp_path: Path) -> None:
    path = _write_kinase_library_fixture(
        tmp_path,
        kinase="ABL1",
        residue_class="tyr",
        family="ABL",
        group="TK",
        amino_acids=("A", "Y"),
    )

    resource = KinaseLibraryResourceLoader(
        source_reader=ReferenceSourceTableReader()
    ).run(path)

    matrix = resource.matrix_for("ABL1", "tyr")
    assert matrix.residue_class is KinaseLibraryResidueClass.TYR
    assert matrix.kinase_group == "TK"
    assert list(matrix.score_table.columns) == [-1, 0, 1]
    assert matrix.score_table.loc["Y", 1] == pytest.approx(4.75)


def test_kinase_library_loader_rejects_invalid_residue_class(
    tmp_path: Path,
) -> None:
    path = _write_kinase_library_fixture(tmp_path, residue_class="ser_tyr")

    with pytest.raises(ReferenceValidationError, match="residue_class"):
        KinaseLibraryResourceLoader(source_reader=ReferenceSourceTableReader()).run(
            path
        )


def test_kinase_library_loader_rejects_missing_positions(
    tmp_path: Path,
) -> None:
    path = _write_kinase_library_fixture(tmp_path, positions=(-1, 1))

    with pytest.raises(ReferenceValidationError, match="missing required positions"):
        KinaseLibraryResourceLoader(source_reader=ReferenceSourceTableReader()).run(
            path
        )


def test_kinase_library_loader_rejects_bad_numeric_scores(
    tmp_path: Path,
) -> None:
    path = _write_kinase_library_fixture(tmp_path, bad_score=True)

    with pytest.raises(ReferenceValidationError, match="score values must be numeric"):
        KinaseLibraryResourceLoader(source_reader=ReferenceSourceTableReader()).run(
            path
        )


def test_kinase_library_loader_preserves_provenance(tmp_path: Path) -> None:
    path = _write_kinase_library_fixture(tmp_path)

    resource = KinaseLibraryResourceLoader(
        source_reader=ReferenceSourceTableReader()
    ).run(path)

    provenance = resource.provenance
    assert provenance.source_type == "local"
    assert provenance.source_name == "Synthetic Kinase Library"
    assert provenance.source_version == "synthetic-kl-v1"
    assert provenance.license == "CC0 synthetic fixture"
    assert provenance.score_scale == "raw_log2_enrichment"
    assert provenance.sequence_window == {
        "upstream_residues": 1,
        "downstream_residues": 1,
        "central_residue_required": True,
    }
    source_file = provenance.source_files["kinase_library"]
    assert isinstance(source_file, dict)
    assert source_file["role"] == "kinase_library"
    assert source_file["sha256"]
    assert provenance.manifest is not None
    assert provenance.manifest["license"] == "CC0 synthetic fixture"
    fingerprint_names = {item.name for item in provenance.table_fingerprints}
    assert "references.kinase_library.matrix_index" in fingerprint_names
    assert "references.kinase_library.score_table.akt1.ser_thr" in fingerprint_names


def _write_kinase_library_fixture(
    tmp_path: Path,
    *,
    kinase: str = "AKT1",
    residue_class: str = "ser_thr",
    family: str = "AGC",
    group: str = "AGC",
    amino_acids: tuple[str, ...] = ("A", "S"),
    positions: tuple[int, ...] = (-1, 0, 1),
    bad_score: bool = False,
) -> Path:
    path = tmp_path / "kinase_library.csv"
    rows: list[dict[str, object]] = []
    for amino_acid in amino_acids:
        for position in positions:
            score: object = _score_for(amino_acid, position)
            if bad_score and amino_acid == amino_acids[0] and position == positions[0]:
                score = "not_numeric"
            rows.append(
                {
                    "kinase": kinase,
                    "kinase_family": family,
                    "kinase_group": group,
                    "residue_class": residue_class,
                    "position": position,
                    "amino_acid": amino_acid,
                    "score": score,
                    "source_name": "Synthetic Kinase Library",
                    "source_version": "synthetic-kl-v1",
                    "retrieved_at": "2026-06-11",
                    "license": "CC0 synthetic fixture",
                    "score_scale": "raw_log2_enrichment",
                    "organisms": "human|mouse",
                    "upstream_residues": 1,
                    "downstream_residues": 1,
                    "central_residue_required": "true",
                }
            )
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def _score_for(amino_acid: str, position: int) -> float:
    scores = {
        ("A", -1): 12.5,
        ("A", 0): -0.5,
        ("A", 1): 1.5,
        ("S", -1): 0.25,
        ("S", 0): 3.25,
        ("S", 1): 2.25,
        ("Y", -1): 0.75,
        ("Y", 0): 5.5,
        ("Y", 1): 4.75,
    }
    return scores[(amino_acid, position)]
