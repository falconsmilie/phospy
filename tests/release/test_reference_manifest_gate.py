from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path

import pytest

from phospy.science.references.errors import ReferenceManifestError
from phospy.science.references.validation import validate_bundled_reference_manifests

pytestmark = pytest.mark.release_gate


def test_current_rat_manifest_fails_release_gate_until_approval_evidence_is_recorded() -> (
    None
):
    with pytest.raises(ReferenceManifestError) as exc_info:
        validate_bundled_reference_manifests(_reference_bundles_root())

    message = str(exc_info.value)
    assert "reference_id='l6_native'" in message
    assert "field='redistribution_notes'" in message
    assert "redistribution_status='approved'" in message
    assert "not independently verified" in message


def test_reference_manifest_release_gate_fails_for_invalid_fixture(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"license_name": ""},
    )

    with pytest.raises(ReferenceManifestError, match="license_name"):
        validate_bundled_reference_manifests(root)


def test_release_gate_rejects_approved_manifest_without_redistribution_evidence(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(tmp_path, remove_redistribution_evidence=True)

    with pytest.raises(ReferenceManifestError) as exc_info:
        validate_bundled_reference_manifests(root)

    message = str(exc_info.value)
    assert "field='redistribution_evidence'" in message
    assert "requires structured exact-file redistribution evidence" in message


def test_release_gate_rejects_approved_manifest_with_unverified_notes(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "redistribution_notes": (
                "Approval for this exact bundle is not independently verified."
            ),
        },
    )

    with pytest.raises(ReferenceManifestError) as exc_info:
        validate_bundled_reference_manifests(root)

    message = str(exc_info.value)
    assert "field='redistribution_notes'" in message
    assert "contradictory approval text: 'not independently verified'" in message


def test_release_gate_error_identifies_reference_organism_namespace_field_and_status(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "redistribution_notes": (
                "Approval for this exact bundle is not independently verified."
            ),
        },
    )

    with pytest.raises(ReferenceManifestError) as exc_info:
        validate_bundled_reference_manifests(root)

    message = str(exc_info.value)
    assert "reference release gate failed:" in message
    assert "reference_id='unit_reference'" in message
    assert "display_name='Unit reference'" in message
    assert "organism='Rattus norvegicus'" in message
    assert "protein_namespace='display_site_id'" in message
    assert "field='redistribution_notes'" in message
    assert "redistribution_status='approved'" in message


def test_release_gate_accepts_approved_manifest_with_structured_exact_file_evidence(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(tmp_path)

    manifests = validate_bundled_reference_manifests(root)

    assert len(manifests) == 1
    manifest = manifests[0]
    assert manifest.reference_id == "unit_reference"
    assert manifest.redistribution_evidence is not None
    assert manifest.redistribution_evidence.applies_to_exact_packaged_files is True
    assert manifest.redistribution_evidence.verified_at.isoformat() == "2026-06-29"
    assert manifest.to_payload()["redistribution_evidence"] == {
        "evidence_type": "synthetic_fixture",
        "applies_to_exact_packaged_files": True,
        "evidence_url": "https://example.test/approval-record",
        "evidence_reference": (
            "unit test approval record for exact packaged synthetic fixture"
        ),
        "verified_at": "2026-06-29",
    }


def test_release_gate_rejects_external_only_bundled_manifest(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "license_name": None,
            "license_url": None,
            "redistribution_status": "external_only",
            "redistribution_notes": "source must be supplied externally",
        },
        remove_redistribution_evidence=True,
    )

    with pytest.raises(ReferenceManifestError) as exc_info:
        validate_bundled_reference_manifests(root)

    message = str(exc_info.value)
    assert "field='redistribution_status'" in message
    assert "redistribution_status='external_only'" in message
    assert "requires redistribution_status 'approved'" in message


def test_release_gate_rejects_unresolved_bundled_manifest(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "license_name": None,
            "license_url": None,
            "redistribution_status": "unresolved",
            "redistribution_notes": "redistribution review has not completed",
        },
        remove_redistribution_evidence=True,
    )

    with pytest.raises(ReferenceManifestError) as exc_info:
        validate_bundled_reference_manifests(root)

    message = str(exc_info.value)
    assert "field='redistribution_status'" in message
    assert "redistribution_status='unresolved'" in message
    assert "requires redistribution_status 'approved'" in message


def _reference_bundles_root() -> Path:
    return (
        Path(__file__).resolve().parents[2]
        / "src"
        / "phospy"
        / "data"
        / "reference_bundles"
    )


def _write_manifest_bundle(
    tmp_path: Path,
    *,
    manifest_overrides: dict[str, object] | None = None,
    remove_redistribution_evidence: bool = False,
) -> Path:
    bundle_root = tmp_path / "reference_bundles" / "rat" / "unit_reference"
    bundle_root.mkdir(parents=True)
    data = "kinase,site_id\nAKT1,MAPK1;S123;\n"
    data_path = bundle_root / "substrate_map.csv"
    data_path.write_text(data, encoding="utf-8")
    payload = _valid_manifest_payload(sha256(data_path.read_bytes()).hexdigest())
    if manifest_overrides is not None:
        payload.update(manifest_overrides)
    if remove_redistribution_evidence:
        payload.pop("redistribution_evidence", None)
    (bundle_root / "manifest.json").write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )
    return tmp_path / "reference_bundles"


def _valid_manifest_payload(file_hash: str) -> dict[str, object]:
    return {
        "reference_id": "unit_reference",
        "display_name": "Unit reference",
        "organism": "Rattus norvegicus",
        "taxonomy_id": 10116,
        "protein_namespace": "display_site_id",
        "reference_version": "v1",
        "source_name": "unit source",
        "source_version": None,
        "source_url": "https://example.test/reference",
        "retrieved_at": "2026-06-29",
        "table_sha256": file_hash,
        "source_publication": None,
        "license_name": "CC0 synthetic",
        "license_url": "https://example.test/license",
        "redistribution_status": "approved",
        "redistribution_notes": "redistribution approved for exact synthetic fixture",
        "redistribution_evidence": {
            "evidence_type": "synthetic_fixture",
            "applies_to_exact_packaged_files": True,
            "evidence_url": "https://example.test/approval-record",
            "evidence_reference": (
                "unit test approval record for exact packaged synthetic fixture"
            ),
            "verified_at": "2026-06-29",
        },
        "derived_from": ["unit test"],
        "generated_by": "unit test",
        "generated_at_utc": "2026-06-29T00:00:00Z",
        "manifest_schema_version": "1.1",
        "files": [
            {
                "relative_path": "substrate_map.csv",
                "role": "kinase_substrate",
                "format": "csv",
                "sha256": file_hash,
                "row_count": 1,
                "column_names": ["kinase", "site_id"],
            }
        ],
        "sequence_context_policy": "centered phosphosite sequence window",
        "sequence_window_length": 3,
        "sequence_center_index": 1,
        "allowed_sequence_alphabet": "ACDEFGHIKLMNPQRSTVWY",
    }
