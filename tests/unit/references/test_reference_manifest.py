from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path

import pytest

from phospy.science.references.errors import ReferenceManifestError
from phospy.science.references.models import RedistributionStatus
from phospy.science.references.validation import load_reference_manifest


def test_valid_approved_manifest_loads_successfully(tmp_path: Path) -> None:
    manifest_path = _write_manifest_bundle(tmp_path)

    manifest = load_reference_manifest(manifest_path)

    assert manifest.reference_id == "unit_reference"
    assert manifest.reference_version == "v1"
    assert manifest.table_sha256 == manifest.files[0].sha256
    assert manifest.license_name == "CC0 synthetic"
    assert manifest.license_url == "https://example.test/license"
    assert manifest.redistribution_status is RedistributionStatus.APPROVED
    assert manifest.redistribution_allowed is True
    assert manifest.files[0].relative_path == "reference.csv"
    assert manifest.sequence_window.upstream_residues == 1
    assert manifest.sequence_window.downstream_residues == 1


def test_approved_manifest_parses_structured_redistribution_evidence(
    tmp_path: Path,
) -> None:
    evidence = {
        "evidence_type": "upstream_package_license",
        "upstream_package": {
            "package_name": "UnitPackage",
            "package_version": "1.0.0",
            "license_name": "CC0-1.0",
            "license_url": "https://example.test/license",
        },
        "scope": {
            "reference_id": "unit_reference",
            "reference_version": "v1",
            "applies_to_exact_packaged_files": True,
            "packaged_files": ["reference.csv"],
            "applies_to_future_bundles": False,
        },
        "attribution": {
            "repository_notice_path": "NOTICE.md",
            "bundle_attribution_path": "ATTRIBUTION.md",
        },
        "independent_database_permission_claimed": False,
        "evidence_url": "https://example.test/approval-record",
        "verified_at": "2026-06-29",
        "notes": "unit approval record for exact packaged fixture",
    }
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "license_name": "CC0-1.0",
            "redistribution_evidence": evidence,
        },
    )

    manifest = load_reference_manifest(manifest_path)

    assert manifest.redistribution_evidence is not None
    assert (
        manifest.redistribution_evidence.evidence_type.value
        == "upstream_package_license"
    )
    assert manifest.redistribution_evidence.upstream_package.package_name == (
        "UnitPackage"
    )
    assert manifest.redistribution_evidence.scope.applies_to_exact_packaged_files is (
        True
    )
    assert manifest.redistribution_evidence.scope.packaged_files == ("reference.csv",)
    assert manifest.redistribution_evidence.evidence_url == (
        "https://example.test/approval-record"
    )
    assert manifest.redistribution_evidence.notes == (
        "unit approval record for exact packaged fixture"
    )
    assert (
        manifest.redistribution_evidence.independent_database_permission_claimed
        is False
    )
    assert manifest.redistribution_evidence.verified_at.isoformat() == "2026-06-29"
    assert manifest.to_payload()["redistribution_evidence"] == evidence


def test_raw_redistribution_allowed_absence_remains_supported(
    tmp_path: Path,
) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        remove_manifest_field="redistribution_allowed",
    )

    manifest = load_reference_manifest(manifest_path)

    assert manifest.raw_redistribution_allowed is None
    assert manifest.redistribution_allowed is True


def test_raw_redistribution_allowed_null_is_rejected(tmp_path: Path) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"redistribution_allowed": None},
    )

    with pytest.raises(ReferenceManifestError) as exc_info:
        load_reference_manifest(manifest_path)

    message = str(exc_info.value)
    assert "redistribution_allowed" in message
    assert "JSON Boolean" in message
    assert "None" in message


@pytest.mark.parametrize(
    "raw_value",
    [
        pytest.param(None, id="null"),
        pytest.param("true", id="string-true"),
        pytest.param("false", id="string-false"),
        pytest.param(0, id="integer-zero"),
        pytest.param(1, id="integer-one"),
        pytest.param(0.0, id="float-zero"),
        pytest.param([], id="array"),
        pytest.param({}, id="object"),
    ],
)
def test_raw_redistribution_allowed_non_boolean_values_are_rejected(
    tmp_path: Path,
    raw_value: object,
) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"redistribution_allowed": raw_value},
    )

    with pytest.raises(ReferenceManifestError) as exc_info:
        load_reference_manifest(manifest_path)

    message = str(exc_info.value)
    assert "redistribution_allowed" in message
    assert "JSON Boolean" in message
    assert type(raw_value).__name__ in message


@pytest.mark.parametrize(
    ("redistribution_status", "raw_value", "expected_allowed"),
    [
        pytest.param("approved", True, True, id="approved-true"),
        pytest.param("unresolved", False, False, id="unresolved-false"),
        pytest.param("external_only", False, False, id="external-only-false"),
    ],
)
def test_raw_redistribution_allowed_status_agreement(
    tmp_path: Path,
    redistribution_status: str,
    raw_value: bool,
    expected_allowed: bool,
) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "license_name": None if not expected_allowed else "CC0 synthetic",
            "license_url": (
                None if not expected_allowed else "https://example.test/license"
            ),
            "redistribution_status": redistribution_status,
            "redistribution_allowed": raw_value,
            "redistribution_notes": "synthetic test fixture",
        },
    )

    manifest = load_reference_manifest(manifest_path)

    assert manifest.raw_redistribution_allowed is raw_value
    assert manifest.redistribution_allowed is expected_allowed


@pytest.mark.parametrize(
    "redistribution_status",
    [
        pytest.param("unresolved", id="unresolved"),
        pytest.param("external_only", id="external-only"),
    ],
)
def test_non_releasable_status_with_true_raw_flag_is_rejected(
    tmp_path: Path,
    redistribution_status: str,
) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "license_name": None,
            "license_url": None,
            "redistribution_status": redistribution_status,
            "redistribution_allowed": True,
            "redistribution_notes": "redistribution review has not completed",
        },
    )

    with pytest.raises(ReferenceManifestError, match="contradicts"):
        load_reference_manifest(manifest_path)


def test_raw_redistribution_allowed_contradiction_raises(tmp_path: Path) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"redistribution_allowed": False},
    )

    with pytest.raises(ReferenceManifestError, match="contradicts"):
        load_reference_manifest(manifest_path)


def test_unrecognized_redistribution_evidence_field_raises(tmp_path: Path) -> None:
    evidence = _structured_evidence()
    evidence["free_text_approval"] = "not authoritative"
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"redistribution_evidence": evidence},
    )

    with pytest.raises(ReferenceManifestError) as exc_info:
        load_reference_manifest(manifest_path)

    message = str(exc_info.value)
    assert "unrecognized field" in message
    assert "redistribution_evidence.free_text_approval" in message


def test_unrecognized_top_level_manifest_field_raises(tmp_path: Path) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"extra_release_note": "reviewed by unit test"},
    )

    with pytest.raises(ReferenceManifestError) as exc_info:
        load_reference_manifest(manifest_path)

    message = str(exc_info.value)
    assert "unrecognized field" in message
    assert "extra_release_note" in message


def test_unrecognized_file_manifest_field_raises(tmp_path: Path) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        file_overrides={"extra_release_note": "reviewed by unit test"},
    )

    with pytest.raises(ReferenceManifestError) as exc_info:
        load_reference_manifest(manifest_path)

    message = str(exc_info.value)
    assert "unrecognized field" in message
    assert "files[0].extra_release_note" in message


def test_valid_external_only_manifest_loads_successfully(tmp_path: Path) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "license_name": None,
            "license_url": None,
            "redistribution_status": "external_only",
            "redistribution_allowed": False,
            "redistribution_notes": "source must be supplied externally",
        },
    )

    manifest = load_reference_manifest(manifest_path)

    assert manifest.redistribution_status is RedistributionStatus.EXTERNAL_ONLY
    assert manifest.redistribution_allowed is False
    assert manifest.license_name is None
    assert manifest.license_url is None


def test_unresolved_manifest_loads_during_ordinary_parsing(tmp_path: Path) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "license_name": None,
            "license_url": None,
            "redistribution_status": "unresolved",
            "redistribution_allowed": False,
            "redistribution_notes": "redistribution review has not completed",
        },
    )

    manifest = load_reference_manifest(manifest_path)

    assert manifest.redistribution_status is RedistributionStatus.UNRESOLVED
    assert manifest.redistribution_allowed is False


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


def test_missing_table_sha256_raises_reference_manifest_error(tmp_path: Path) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"table_sha256": ""},
    )

    with pytest.raises(ReferenceManifestError, match="table_sha256"):
        load_reference_manifest(manifest_path)


def test_missing_license_for_approved_manifest_raises(tmp_path: Path) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"license_name": "", "license_url": None},
    )

    with pytest.raises(ReferenceManifestError, match="license_name and license_url"):
        load_reference_manifest(manifest_path)


def test_invalid_redistribution_status_raises(tmp_path: Path) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"redistribution_status": "redistributable"},
    )

    with pytest.raises(ReferenceManifestError, match="redistribution_status"):
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


def test_unresolved_bundled_manifest_is_rejected(
    tmp_path: Path,
) -> None:
    manifest_path = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "license_name": None,
            "license_url": None,
            "redistribution_status": "unresolved",
            "redistribution_allowed": False,
            "redistribution_notes": "redistribution review has not completed",
        },
    )

    with pytest.raises(ReferenceManifestError, match="redistribution_status"):
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
        "source_url": "https://example.test/reference",
        "retrieved_at": "2026-06-29",
        "table_sha256": file_hash,
        "source_publication": None,
        "license_name": "CC0 synthetic",
        "license_url": "https://example.test/license",
        "redistribution_status": "approved",
        "redistribution_allowed": True,
        "redistribution_notes": "synthetic test fixture",
        "derived_from": ["unit test"],
        "generated_by": "unit test",
        "generated_at_utc": "2026-06-29T00:00:00Z",
        "manifest_schema_version": "1.1",
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


def _structured_evidence() -> dict[str, object]:
    return {
        "evidence_type": "upstream_package_license",
        "upstream_package": {
            "package_name": "UnitPackage",
            "package_version": "1.0.0",
            "license_name": "CC0 synthetic",
            "license_url": "https://example.test/license",
        },
        "scope": {
            "reference_id": "unit_reference",
            "reference_version": "v1",
            "applies_to_exact_packaged_files": True,
            "packaged_files": ["reference.csv"],
            "applies_to_future_bundles": False,
        },
        "attribution": {
            "repository_notice_path": "NOTICE.md",
            "bundle_attribution_path": "ATTRIBUTION.md",
        },
        "independent_database_permission_claimed": False,
        "evidence_url": "https://example.test/approval-record",
        "verified_at": "2026-06-29",
        "notes": "unit approval record for exact packaged fixture",
    }
