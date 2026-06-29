from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path

import pytest

from phospy.science.references.errors import ReferenceManifestError
from phospy.science.references.validation import load_reference_manifest


def test_valid_manifest_loads_successfully(tmp_path: Path) -> None:
    manifest_path = _write_manifest_bundle(tmp_path)

    manifest = load_reference_manifest(manifest_path)

    assert manifest.reference_id == "unit_reference"
    assert manifest.reference_version == "v1"
    assert manifest.redistribution_allowed is True
    assert manifest.files[0].relative_path == "reference.csv"
    assert manifest.sequence_window.upstream_residues == 1
    assert manifest.sequence_window.downstream_residues == 1


def test_missing_required_field_raises_reference_manifest_error(
    tmp_path: Path,
) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path, remove_manifest_field="source_name"
    )

    with pytest.raises(ReferenceManifestError, match="source_name"):
        load_reference_manifest(manifest_path)


def test_missing_manifest_file_raises_reference_manifest_error(tmp_path: Path) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        file_overrides={"relative_path": "missing.csv"},
    )

    with pytest.raises(ReferenceManifestError, match="does not exist"):
        load_reference_manifest(manifest_path)


def test_sha256_mismatch_raises_reference_manifest_error(tmp_path: Path) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        file_overrides={"sha256": "0" * 64},
    )

    with pytest.raises(ReferenceManifestError, match="hash mismatch"):
        load_reference_manifest(manifest_path)


def test_sequence_window_length_without_center_index_raises(
    tmp_path: Path,
) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"sequence_window_length": 3, "sequence_center_index": None},
    )

    with pytest.raises(ReferenceManifestError, match="without sequence_center_index"):
        load_reference_manifest(manifest_path)


def test_sequence_center_index_outside_window_length_raises(tmp_path: Path) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"sequence_window_length": 3, "sequence_center_index": 3},
    )

    with pytest.raises(ReferenceManifestError, match="within sequence_window_length"):
        load_reference_manifest(manifest_path)


def test_bundled_manifest_with_disallowed_redistribution_raises(
    tmp_path: Path,
) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "redistribution_allowed": False,
            "redistribution_notes": "local/private use only",
        },
    )

    with pytest.raises(ReferenceManifestError, match="package release is blocked"):
        load_reference_manifest(manifest_path, bundled=True)


def _write_manifest_bundle(
    tmp_path: Path,
    *,
    manifest_overrides: dict[str, object] | None = None,
    file_overrides: dict[str, object] | None = None,
    remove_manifest_field: str | None = None,
) -> Path:
    data_path = tmp_path / "reference.csv"
    data = "kinase,site_id\nAKT1,MAPK1;S123;\n"
    data_path.write_text(data, encoding="utf-8")
    file_hash = sha256(data_path.read_bytes()).hexdigest()
    file_payload: dict[str, object] = {
        "relative_path": "reference.csv",
        "role": "kinase_substrate",
        "format": "csv",
        "sha256": file_hash,
        "row_count": 1,
        "column_names": ["kinase", "site_id"],
    }
    if file_overrides is not None:
        file_payload.update(file_overrides)
    manifest_payload: dict[str, object] = {
        "reference_id": "unit_reference",
        "display_name": "Unit reference",
        "organism": "Rattus norvegicus",
        "taxonomy_id": 10116,
        "protein_namespace": "display_site_id",
        "reference_version": "v1",
        "source_name": "unit source",
        "source_version": "source-v1",
        "source_url": None,
        "source_publication": None,
        "source_license": "CC0 synthetic",
        "source_license_url": None,
        "redistribution_allowed": True,
        "redistribution_notes": "synthetic test fixture",
        "derived_from": ["unit test"],
        "generated_by": "unit test",
        "generated_at_utc": "2026-06-29T00:00:00Z",
        "manifest_schema_version": "1.0",
        "files": [file_payload],
        "sequence_context_policy": "centered phosphosite sequence window",
        "sequence_window_length": 3,
        "sequence_center_index": 1,
        "allowed_sequence_alphabet": "ACDEFGHIKLMNPQRSTVWY",
    }
    if manifest_overrides is not None:
        manifest_payload.update(manifest_overrides)
    if remove_manifest_field is not None:
        manifest_payload.pop(remove_manifest_field)
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest_payload, indent=2),
        encoding="utf-8",
    )
    return manifest_path
