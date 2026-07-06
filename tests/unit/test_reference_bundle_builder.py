from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy.api import (
    Organism,
    ReferenceBundle,
    ReferenceBundleBuilder,
    ReferenceBundleBuildRequest,
)
from phospy.errors.references import ReferenceResolutionError
from phospy.errors.validation import ReferenceValidationError


def test_reference_bundle_builder_builds_valid_human_bundle_from_local_files(
    tmp_path: Path,
) -> None:
    kinase_path, sequence_path = _write_reference_sources(
        tmp_path,
        organism_label="human",
    )

    bundle = ReferenceBundleBuilder().run(
        ReferenceBundleBuildRequest(
            organism=Organism.HUMAN,
            kinase_substrate_path=kinase_path,
            site_sequence_path=sequence_path,
            source_name="synthetic local reference",
            source_version="fixture-2026-06-11",
            retrieved_at="2026-06-11",
            license="synthetic test license",
            redistribution_status="redistributable synthetic fixture",
            identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
        )
    )

    assert isinstance(bundle, ReferenceBundle)
    assert bundle.organism is Organism.HUMAN
    assert bundle.kinase_substrate_map.loc[:, "kinase"].tolist() == [
        "AKT1",
        "MAP2K1",
    ]
    assert bundle.kinase_substrate_map.loc[:, "substrate_site"].tolist() == [
        "MAPK1;S123;",
        "MAPK1;Y185;",
    ]
    assert bundle.site_sequences.index.tolist() == ["MAPK1;S123;", "MAPK1;Y185;"]
    assert bundle.site_sequences.loc["MAPK1;S123;", "site_sequence"] == _window("S")
    assert bundle.site_sequences.loc["MAPK1;S123;", "display_id"] == "MAPK1;S123;"
    assert bundle.site_sequences.loc["MAPK1;S123;", "gene_symbol"] == "MAPK1"
    assert bundle.site_sequences.loc["MAPK1;S123;", "protein_accession"] == "P28482"


def test_reference_bundle_builder_builds_valid_mouse_bundle_from_local_files(
    tmp_path: Path,
) -> None:
    kinase_path, sequence_path = _write_reference_sources(
        tmp_path,
        organism_label="Mus musculus",
    )

    bundle = ReferenceBundleBuilder().run(
        ReferenceBundleBuildRequest(
            organism=Organism.MOUSE,
            kinase_substrate_path=kinase_path,
            site_sequence_path=sequence_path,
            source_name="synthetic local mouse reference",
            source_version="fixture-2026-06-11",
            retrieved_at="2026-06-11",
            license="synthetic test license",
            redistribution_status="redistributable synthetic fixture",
            identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
        )
    )

    assert bundle.organism is Organism.MOUSE
    assert bundle.manifest is not None
    assert bundle.manifest.organism_common_name == "mouse"
    assert bundle.provenance is not None
    assert bundle.provenance.source_type == "local"


def test_reference_bundle_builder_provenance_records_source_and_redistribution_metadata(
    tmp_path: Path,
) -> None:
    kinase_path, sequence_path = _write_reference_sources(
        tmp_path,
        organism_label="mouse",
    )

    bundle = ReferenceBundleBuilder().run(
        ReferenceBundleBuildRequest(
            organism=Organism.MOUSE,
            kinase_substrate_path=kinase_path,
            site_sequence_path=sequence_path,
            source_name="synthetic provenance source",
            source_version="v1",
            retrieved_at="2026-06-11",
            license="CC0 synthetic",
            redistribution_status="redistributable",
            identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
            bundle_id="synthetic_mouse_local_v1",
        )
    )

    assert bundle.manifest is not None
    assert bundle.manifest.bundle_id == "synthetic_mouse_local_v1"
    assert bundle.manifest.license == "CC0 synthetic"
    assert bundle.manifest.redistribution_status == "unresolved"
    assert bundle.manifest.redistribution_notes == "redistributable"
    assert bundle.manifest.source_files is not None
    assert set(bundle.manifest.source_files) == {
        "kinase_substrate",
        "site_sequences",
    }
    assert bundle.provenance is not None
    assert bundle.provenance.source_type == "local"
    assert bundle.provenance.source_name == "synthetic provenance source"
    assert bundle.provenance.source_version == "v1"
    assert bundle.provenance.retrieved_at == "2026-06-11"
    assert bundle.provenance.identifier_namespace == (
        "display_id (GENE_SYMBOL;RESIDUE;)"
    )
    assert bundle.provenance.sequence_window == {
        "upstream_residues": 15,
        "downstream_residues": 15,
        "central_residue_required": True,
    }
    assert bundle.provenance.manifest is not None
    assert bundle.provenance.manifest["license"] == "CC0 synthetic"
    assert bundle.provenance.manifest["redistribution_status"] == "unresolved"
    assert bundle.provenance.manifest["redistribution_notes"] == "redistributable"
    source_files = bundle.provenance.manifest["source_files"]
    assert isinstance(source_files, dict)
    assert source_files["kinase_substrate"]["sha256"]
    assert source_files["site_sequences"]["sha256"]
    assert {item.name for item in bundle.provenance.table_fingerprints} == {
        "references.kinase_substrate_map",
        "references.site_sequences",
    }
    assert bundle.provenance.identifier_normalisation is not None


def test_reference_bundle_builder_does_not_mark_approval_without_license_url(
    tmp_path: Path,
) -> None:
    kinase_path, sequence_path = _write_reference_sources(
        tmp_path,
        organism_label="mouse",
    )

    bundle = ReferenceBundleBuilder().run(
        ReferenceBundleBuildRequest(
            organism=Organism.MOUSE,
            kinase_substrate_path=kinase_path,
            site_sequence_path=sequence_path,
            source_name="synthetic approval source",
            source_version="v1",
            retrieved_at="2026-06-11",
            license="CC0 synthetic",
            redistribution_status="approved",
            identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
        )
    )

    assert bundle.manifest is not None
    assert bundle.manifest.redistribution_status == "unresolved"
    assert bundle.manifest.redistribution_notes == "approved"
    assert bundle.manifest.license_url is None


def test_reference_bundle_builder_rejects_missing_sequence_context(
    tmp_path: Path,
) -> None:
    kinase_path = tmp_path / "kinase.csv"
    sequence_path = tmp_path / "sequences.csv"
    pd.DataFrame(
        {
            "kinase": ["AKT1", "MAP2K1"],
            "site_id": ["MAPK1;S123;", "MAPK1;Y185;"],
        }
    ).to_csv(kinase_path, index=False)
    pd.DataFrame(
        {
            "site_id": ["MAPK1;S123;"],
            "site_sequence": [_window("S")],
        }
    ).to_csv(sequence_path, index=False)

    with pytest.raises(ReferenceValidationError, match="missing sequence entries"):
        _build_mouse_bundle(kinase_path, sequence_path)


def test_reference_bundle_builder_rejects_malformed_site_identifiers(
    tmp_path: Path,
) -> None:
    kinase_path = tmp_path / "kinase.csv"
    sequence_path = tmp_path / "sequences.csv"
    pd.DataFrame(
        {
            "kinase": ["AKT1"],
            "site_id": ["MAPK1-S123"],
        }
    ).to_csv(kinase_path, index=False)
    pd.DataFrame(
        {
            "site_id": ["MAPK1;S123;"],
            "site_sequence": [_window("S")],
        }
    ).to_csv(sequence_path, index=False)

    with pytest.raises(ReferenceValidationError, match="GENE;SITE"):
        _build_mouse_bundle(kinase_path, sequence_path)


def test_reference_bundle_builder_rejects_duplicate_kinase_substrate_edges(
    tmp_path: Path,
) -> None:
    kinase_path = tmp_path / "kinase.csv"
    sequence_path = tmp_path / "sequences.csv"
    pd.DataFrame(
        {
            "kinase": ["akt1", "AKT1"],
            "site_id": ["mapk1;s123", "MAPK1;S123;"],
        }
    ).to_csv(kinase_path, index=False)
    pd.DataFrame(
        {
            "site_id": ["MAPK1;S123;"],
            "site_sequence": [_window("S")],
        }
    ).to_csv(sequence_path, index=False)

    with pytest.raises(
        ReferenceValidationError,
        match="duplicate \\(kinase, substrate_site\\) pairs",
    ):
        _build_mouse_bundle(kinase_path, sequence_path)


def test_reference_bundle_builder_rejects_source_organism_mismatch(
    tmp_path: Path,
) -> None:
    kinase_path, sequence_path = _write_reference_sources(
        tmp_path,
        organism_label="mouse",
    )

    with pytest.raises(ReferenceResolutionError, match="does not match requested"):
        ReferenceBundleBuilder().run(
            ReferenceBundleBuildRequest(
                organism=Organism.HUMAN,
                kinase_substrate_path=kinase_path,
                site_sequence_path=sequence_path,
                source_name="synthetic local reference",
                source_version="fixture-2026-06-11",
                retrieved_at="2026-06-11",
                license="synthetic test license",
                redistribution_status="redistributable synthetic fixture",
                identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
            )
        )


def test_reference_bundle_builder_rejects_incomplete_kinase_maps_with_diagnostics(
    tmp_path: Path,
) -> None:
    kinase_path = tmp_path / "kinase.csv"
    sequence_path = tmp_path / "sequences.csv"
    pd.DataFrame({"site_id": ["MAPK1;S123;"]}).to_csv(kinase_path, index=False)
    pd.DataFrame(
        {
            "site_id": ["MAPK1;S123;"],
            "site_sequence": [_window("S")],
        }
    ).to_csv(sequence_path, index=False)

    with pytest.raises(ReferenceResolutionError) as exc_info:
        _build_mouse_bundle(kinase_path, sequence_path)
    message = str(exc_info.value)
    assert "missing required kinase column" in message
    assert "accepted aliases" in message
    assert "available columns: site_id" in message


def _write_reference_sources(
    tmp_path: Path,
    *,
    organism_label: str,
) -> tuple[Path, Path]:
    kinase_path = tmp_path / "kinase.csv"
    sequence_path = tmp_path / "sequences.csv"
    pd.DataFrame(
        {
            "kinase": [" akt1 ", "Map2k1"],
            "site_id": [" mapk1 ; s123 ", "Mapk1;Y185;"],
            "organism": [organism_label, organism_label],
        }
    ).to_csv(kinase_path, index=False)
    pd.DataFrame(
        {
            "site_id": [" mapk1 ; s123 ", "Mapk1;Y185;"],
            "site_sequence": [_window("S").lower(), _window("Y")],
            "display_id": ["mapk1;s123", "MAPK1;Y185;"],
            "gene": [" mapk1 ", "Mapk1"],
            "protein_accession": [" P28482 ", "Q63844"],
            "organism": [organism_label, organism_label],
        }
    ).to_csv(sequence_path, index=False)
    return kinase_path, sequence_path


def _build_mouse_bundle(kinase_path: Path, sequence_path: Path) -> ReferenceBundle:
    return ReferenceBundleBuilder().run(
        ReferenceBundleBuildRequest(
            organism=Organism.MOUSE,
            kinase_substrate_path=kinase_path,
            site_sequence_path=sequence_path,
            source_name="synthetic local reference",
            source_version="fixture-2026-06-11",
            retrieved_at="2026-06-11",
            license="synthetic test license",
            redistribution_status="redistributable synthetic fixture",
            identifier_namespace="display_id (GENE_SYMBOL;RESIDUE;)",
        )
    )


def _window(center: str) -> str:
    return f"{'A' * 15}{center}{'A' * 15}"
