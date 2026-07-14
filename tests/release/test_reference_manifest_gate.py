from __future__ import annotations

import json
from copy import deepcopy
from datetime import date
from hashlib import sha256
from pathlib import Path

import pytest

from phospy.science.references.errors import ReferenceManifestError
from phospy.science.references.manifest import RedistributionEvidenceType
from phospy.science.references.validation import validate_bundled_reference_manifests

pytestmark = pytest.mark.release_gate


def test_bundled_rat_l6_native_manifest_is_release_eligible_with_phosr_1_20_0_license_evidence() -> (
    None
):
    manifests = validate_bundled_reference_manifests(_reference_bundles_root())

    rat_manifest = next(
        manifest for manifest in manifests if manifest.reference_id == "l6_native"
    )
    assert rat_manifest.source_version == "PhosR 1.20.0"
    assert rat_manifest.raw_redistribution_allowed is True
    assert rat_manifest.redistribution_evidence is not None
    evidence = rat_manifest.redistribution_evidence
    assert evidence.evidence_type is (
        RedistributionEvidenceType.UPSTREAM_PACKAGE_LICENSE
    )
    assert evidence.upstream_package.package_name == "PhosR"
    assert evidence.upstream_package.package_version == "1.20.0"
    assert evidence.upstream_package.license_name == "GPL-3 + file LICENSE"
    assert evidence.scope.reference_id == "l6_native"
    assert evidence.scope.reference_version == "bundled-snapshot-2026-04-16"
    assert evidence.scope.applies_to_exact_packaged_files is True
    assert evidence.scope.applies_to_future_bundles is False
    assert set(evidence.scope.packaged_files) == {
        item.relative_path for item in rat_manifest.files
    }
    assert evidence.attribution.repository_notice_path == "NOTICE.md"
    assert evidence.attribution.bundle_attribution_path == "ATTRIBUTION.md"
    assert evidence.independent_database_permission_claimed is False


def test_real_rat_manifest_keeps_expected_verification_date() -> None:
    manifests = validate_bundled_reference_manifests(_reference_bundles_root())

    rat_manifest = next(
        manifest for manifest in manifests if manifest.reference_id == "l6_native"
    )

    assert rat_manifest.redistribution_evidence is not None
    assert rat_manifest.redistribution_evidence.verified_at == date(2026, 4, 16)


def test_release_gate_accepts_approved_manifest_with_structured_exact_file_evidence(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(tmp_path)

    manifests = validate_bundled_reference_manifests(root)

    assert len(manifests) == 1
    manifest = manifests[0]
    assert manifest.reference_id == "unit_reference"
    assert manifest.redistribution_evidence is not None
    evidence = manifest.redistribution_evidence
    assert evidence.evidence_type is (
        RedistributionEvidenceType.UPSTREAM_PACKAGE_LICENSE
    )
    assert evidence.scope.applies_to_exact_packaged_files is True
    assert evidence.scope.packaged_files == ("substrate_map.csv", "ATTRIBUTION.md")
    assert evidence.verified_at is not None
    assert evidence.verified_at.isoformat() == "2026-06-29"
    assert manifest.to_payload()["redistribution_evidence"] == _valid_evidence_payload(
        reference_id="unit_reference",
        reference_version="v1",
    )


def test_release_gate_rejects_approved_manifest_without_redistribution_evidence(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(tmp_path, remove_redistribution_evidence=True)

    message = _release_error(root)

    assert "field='redistribution_evidence'" in message
    assert "actual_value=None" in message
    assert "requires structured exact-file redistribution evidence" in message


def test_release_gate_requires_verified_at_for_approved_evidence(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(tmp_path, remove_evidence_field="verified_at")

    message = _release_error(root)

    assert "field='redistribution_evidence.verified_at'" in message
    assert "actual_value=None" in message
    assert "missing or null verified_at" in message
    assert "requires an explicit verification date" in message


def test_release_gate_rejects_null_verified_at(tmp_path: Path) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        evidence_overrides={"verified_at": None},
    )

    message = _release_error(root)

    assert "field='redistribution_evidence.verified_at'" in message
    assert "actual_value=None" in message
    assert "missing or null verified_at" in message
    assert "requires an explicit verification date" in message


@pytest.mark.parametrize("source_version", [None, " "])
def test_release_gate_rejects_missing_or_blank_source_version(
    tmp_path: Path,
    source_version: object,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"source_version": source_version},
    )

    message = _release_error(root)

    assert "field='source_version'" in message
    assert "non-empty string" in message


@pytest.mark.parametrize(
    "source_version",
    ["unknown", "Unspecified", "n/a", "na", "none", "null", "tbd", "not specified"],
)
def test_release_gate_rejects_placeholder_source_version(
    tmp_path: Path,
    source_version: str,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"source_version": source_version},
    )

    message = _release_error(root)

    assert "field='source_version'" in message
    assert "must not be a placeholder" in message
    assert f"actual_value={source_version!r}" in message


def test_rat_policy_rejects_wrong_source_version(tmp_path: Path) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        rat_policy=True,
        manifest_overrides={"source_version": "PhosR 1.22.0"},
    )

    message = _release_error(root)

    assert "field='source_version'" in message
    assert "expected 'PhosR 1.20.0'; got 'PhosR 1.22.0'" in message


def test_rat_policy_rejects_wrong_upstream_package(tmp_path: Path) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        rat_policy=True,
        upstream_overrides={"package_name": "OtherPackage"},
    )

    message = _release_error(root)

    assert "field='redistribution_evidence.upstream_package.package_name'" in message
    assert "expected 'PhosR'; got 'OtherPackage'" in message


def test_rat_policy_rejects_wrong_upstream_package_version(tmp_path: Path) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        rat_policy=True,
        upstream_overrides={"package_version": "1.22.0"},
    )

    message = _release_error(root)

    assert "field='redistribution_evidence.upstream_package.package_version'" in message
    assert "expected '1.20.0'; got '1.22.0'" in message


def test_release_gate_rejects_wrong_evidence_license(tmp_path: Path) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        upstream_overrides={"license_name": "see redistribution notes"},
    )

    message = _release_error(root)

    assert "field='redistribution_evidence.upstream_package.license_name'" in message
    assert "machine-readable" in message


def test_release_gate_rejects_evidence_manifest_license_mismatch(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        upstream_overrides={"license_name": "Apache-2.0"},
    )

    message = _release_error(root)

    assert "field='redistribution_evidence.upstream_package.license_name'" in message
    assert "expected manifest license_name 'MIT'; got 'Apache-2.0'" in message


def test_release_gate_rejects_scope_reference_id_mismatch(tmp_path: Path) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        scope_overrides={"reference_id": "other_reference"},
    )

    message = _release_error(root)

    assert "field='redistribution_evidence.scope.reference_id'" in message
    assert "expected 'unit_reference'; got 'other_reference'" in message


def test_release_gate_rejects_scope_reference_version_mismatch(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        scope_overrides={"reference_version": "v2"},
    )

    message = _release_error(root)

    assert "field='redistribution_evidence.scope.reference_version'" in message
    assert "expected 'v1'; got 'v2'" in message


def test_release_gate_rejects_missing_exact_file_scope(tmp_path: Path) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        scope_overrides={"applies_to_exact_packaged_files": False},
    )

    message = _release_error(root)

    assert "field='redistribution_evidence.scope.applies_to_exact_packaged_files'" in (
        message
    )
    assert "must apply to exact packaged files" in message


def test_release_gate_rejects_incomplete_packaged_file_list(tmp_path: Path) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        scope_overrides={"packaged_files": ["substrate_map.csv"]},
    )

    message = _release_error(root)

    assert "field='redistribution_evidence.scope.packaged_files'" in message
    assert "missing=['ATTRIBUTION.md']" in message


def test_release_gate_rejects_extra_packaged_file_entry(tmp_path: Path) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        scope_overrides={
            "packaged_files": [
                "substrate_map.csv",
                "ATTRIBUTION.md",
                "extra.csv",
            ]
        },
    )

    message = _release_error(root)

    assert "field='redistribution_evidence.scope.packaged_files'" in message
    assert "extra=['extra.csv']" in message


def test_release_gate_rejects_duplicate_packaged_file_entry(tmp_path: Path) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        scope_overrides={
            "packaged_files": [
                "substrate_map.csv",
                "ATTRIBUTION.md",
                "ATTRIBUTION.md",
            ]
        },
    )

    message = _release_error(root)

    assert "field='redistribution_evidence.scope.packaged_files'" in message
    assert "must not contain duplicates" in message


@pytest.mark.parametrize(
    "packaged_file",
    ["C:/absolute.csv", "/absolute.csv", "../secret.csv", "nested\\file.csv"],
)
def test_release_gate_rejects_absolute_or_parent_traversal_packaged_file_path(
    tmp_path: Path,
    packaged_file: str,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        scope_overrides={
            "packaged_files": [
                "substrate_map.csv",
                "ATTRIBUTION.md",
                packaged_file,
            ]
        },
    )

    message = _release_error(root)

    assert "redistribution_evidence.scope.packaged_files" in message
    assert f"actual_value={packaged_file!r}" in message


def test_release_gate_rejects_missing_bundle_attribution_location(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        files_without_attribution=True,
        scope_overrides={"packaged_files": ["substrate_map.csv"]},
    )

    message = _release_error(root)

    assert "field='redistribution_evidence.attribution.bundle_attribution_path'" in (
        message
    )
    assert "must be included in manifest.files" in message


def test_release_gate_rejects_missing_repository_notice_file(tmp_path: Path) -> None:
    root = _write_manifest_bundle(tmp_path, omit_notice_file=True)

    message = _release_error(root)

    assert "field='redistribution_evidence.attribution.repository_notice_path'" in (
        message
    )
    assert "repository notice file does not exist" in message


def test_release_gate_rejects_missing_bundle_attribution_file(tmp_path: Path) -> None:
    root = _write_manifest_bundle(tmp_path, omit_bundle_attribution_file=True)

    message = _release_error(root)

    assert "field='redistribution_evidence.attribution.bundle_attribution_path'" in (
        message
    )
    assert "bundle attribution file does not exist" in message


def test_release_gate_rejects_future_bundle_scope(tmp_path: Path) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        scope_overrides={"applies_to_future_bundles": True},
    )

    message = _release_error(root)

    assert "field='redistribution_evidence.scope.applies_to_future_bundles'" in (
        message
    )
    assert "must not cover future bundles" in message


@pytest.mark.parametrize(
    ("manifest_overrides", "evidence_overrides", "expected_field", "phrase"),
    [
        (
            {"redistribution_notes": "This is not legal approval."},
            None,
            "redistribution_notes",
            "not legal approval",
        ),
        (
            None,
            {"notes": "The exact packaged bundle has not been approved."},
            "redistribution_evidence.notes",
            "the exact packaged bundle has not been approved",
        ),
        (
            {"limitations": ["redistribution remains unresolved."]},
            None,
            "limitations[0]",
            "redistribution remains unresolved",
        ),
        (
            None,
            {"notes": "Approval has not been independently verified."},
            "redistribution_evidence.notes",
            "approval has not been independently verified",
        ),
    ],
)
def test_release_gate_rejects_contradictory_approval_text_in_any_string_field(
    tmp_path: Path,
    manifest_overrides: dict[str, object] | None,
    evidence_overrides: dict[str, object] | None,
    expected_field: str,
    phrase: str,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides=manifest_overrides,
        evidence_overrides=evidence_overrides,
    )

    message = _release_error(root)

    assert f"field='{expected_field}'" in message
    assert f"contradictory approval text: {phrase!r}" in message


def test_manifest_rejects_unknown_top_level_field(tmp_path: Path) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"extra_release_note": "reviewed by unit test"},
    )

    message = _release_error(root)

    assert "unrecognized field" in message
    assert "extra_release_note" in message


def test_manifest_rejects_unknown_file_entry_field(tmp_path: Path) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        file_overrides={"extra_release_note": "reviewed by unit test"},
    )

    message = _release_error(root)

    assert "unrecognized field" in message
    assert "files[0].extra_release_note" in message


def test_unknown_top_level_field_cannot_hide_approval_contradiction(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "extra_release_note": "The exact packaged bundle has not been approved."
        },
    )

    message = _release_error(root)

    assert "unrecognized field" in message
    assert "extra_release_note" in message
    assert "contradictory approval text" not in message


def test_unknown_file_field_cannot_hide_approval_contradiction(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        file_overrides={
            "extra_release_note": "The exact packaged bundle has not been approved."
        },
    )

    message = _release_error(root)

    assert "unrecognized field" in message
    assert "files[0].extra_release_note" in message
    assert "contradictory approval text" not in message


def test_known_limitations_field_still_rejects_approval_contradiction(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"limitations": ["Redistribution remains unresolved."]},
    )

    message = _release_error(root)

    assert "field='limitations[0]'" in message
    assert "contradictory approval text: 'redistribution remains unresolved'" in message


def test_scientific_unknown_language_remains_allowed(tmp_path: Path) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "limitations": ["Coverage of unknown kinase substrates remains incomplete."]
        },
    )

    manifests = validate_bundled_reference_manifests(root)

    assert len(manifests) == 1


def test_exact_snapshot_scope_language_remains_allowed(tmp_path: Path) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        evidence_overrides={
            "notes": (
                "No independent direct permission from PhosphoSitePlus is claimed. "
                "Approval applies only to this exact PhosR-derived snapshot. "
                "This does not apply to future bundles."
            )
        },
    )

    manifests = validate_bundled_reference_manifests(root)

    assert len(manifests) == 1


def test_release_gate_accepts_limited_scope_and_no_independent_permission_language(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        evidence_overrides={
            "notes": (
                "No independent direct permission from PhosphoSitePlus is claimed. "
                "Approval applies only to this exact PhosR-derived snapshot. "
                "This does not apply to future bundles."
            )
        },
    )

    manifests = validate_bundled_reference_manifests(root)

    assert len(manifests) == 1


def test_release_gate_accepts_scientific_unknown_wording(tmp_path: Path) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "limitations": ["Coverage of unknown kinase substrates remains incomplete."]
        },
    )

    manifests = validate_bundled_reference_manifests(root)

    assert len(manifests) == 1


def test_release_gate_rejects_affirmative_database_permission_claim(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "redistribution_notes": (
                "Independent direct permission from PhosphoSitePlus is claimed."
            )
        },
    )

    message = _release_error(root)

    assert "field='redistribution_notes'" in message
    assert "claims independent direct database permission" in message


def test_release_gate_rejects_approved_raw_redistribution_allowed_false(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"redistribution_allowed": False},
    )

    message = _release_error(root)

    assert "field='redistribution_allowed'" in message
    assert "actual_value=False" in message
    assert "must be true when redistribution_status is 'approved'" in message


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
def test_release_gate_rejects_raw_redistribution_allowed_non_boolean_values(
    tmp_path: Path,
    raw_value: object,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={"redistribution_allowed": raw_value},
    )

    message = _release_error(root)

    assert "redistribution_allowed" in message
    assert "JSON Boolean" in message
    assert type(raw_value).__name__ in message


@pytest.mark.parametrize(
    "redistribution_status",
    [
        pytest.param("unresolved", id="unresolved"),
        pytest.param("external_only", id="external-only"),
    ],
)
def test_release_gate_rejects_non_releasable_status_with_true_raw_flag(
    tmp_path: Path,
    redistribution_status: str,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "license_name": None,
            "license_url": None,
            "redistribution_status": redistribution_status,
            "redistribution_allowed": True,
            "redistribution_notes": "redistribution review has not completed",
        },
        remove_redistribution_evidence=True,
    )

    message = _release_error(root)

    assert "field='redistribution_allowed'" in message
    assert "actual_value=True" in message
    assert "must be false for non-releasable redistribution_status values" in message


def test_release_gate_rejects_external_only_bundled_manifest(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        manifest_overrides={
            "license_name": None,
            "license_url": None,
            "redistribution_status": "external_only",
            "redistribution_allowed": False,
            "redistribution_notes": "source must be supplied externally",
        },
        remove_redistribution_evidence=True,
    )

    message = _release_error(root)

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
            "redistribution_allowed": False,
            "redistribution_notes": "redistribution review has not completed",
        },
        remove_redistribution_evidence=True,
    )

    message = _release_error(root)

    assert "field='redistribution_status'" in message
    assert "redistribution_status='unresolved'" in message
    assert "requires redistribution_status 'approved'" in message


def test_release_gate_error_identifies_reference_context_field_status_and_actual_value(
    tmp_path: Path,
) -> None:
    root = _write_manifest_bundle(
        tmp_path,
        rat_policy=True,
        upstream_overrides={"package_version": "1.22.0"},
    )

    message = _release_error(root)

    assert "Reference release validation failed:" in message
    assert "reference_id='l6_native'" in message
    assert "display_name='Rat L6 fixture'" in message
    assert "organism='Rattus norvegicus'" in message
    assert "namespace='display_site_id'" in message
    assert "field='redistribution_evidence.upstream_package.package_version'" in (
        message
    )
    assert "redistribution_status='approved'" in message
    assert "actual_value='1.22.0'" in message


def _release_error(root: Path) -> str:
    with pytest.raises(ReferenceManifestError) as exc_info:
        validate_bundled_reference_manifests(root)
    return str(exc_info.value)


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
    file_overrides: dict[str, object] | None = None,
    evidence_overrides: dict[str, object] | None = None,
    upstream_overrides: dict[str, object] | None = None,
    scope_overrides: dict[str, object] | None = None,
    attribution_overrides: dict[str, object] | None = None,
    remove_evidence_field: str | None = None,
    remove_redistribution_evidence: bool = False,
    files_without_attribution: bool = False,
    omit_notice_file: bool = False,
    omit_bundle_attribution_file: bool = False,
    rat_policy: bool = False,
) -> Path:
    if not omit_notice_file:
        (tmp_path / "NOTICE.md").write_text("Unit notice\n", encoding="utf-8")
    reference_id = "l6_native" if rat_policy else "unit_reference"
    bundle_root = tmp_path / "reference_bundles" / "rat" / reference_id
    bundle_root.mkdir(parents=True)
    data = "kinase,site_id\nAKT1,MAPK1;S123;\n"
    data_path = bundle_root / "substrate_map.csv"
    data_path.write_text(data, encoding="utf-8")
    attribution_path = bundle_root / "ATTRIBUTION.md"
    if not omit_bundle_attribution_file:
        attribution_path.write_text("Unit attribution\n", encoding="utf-8")
    file_payloads = [
        {
            "relative_path": "substrate_map.csv",
            "role": "kinase_substrate",
            "format": "csv",
            "sha256": sha256(data_path.read_bytes()).hexdigest(),
            "row_count": 1,
            "column_names": ["kinase", "site_id"],
        }
    ]
    if file_overrides is not None:
        file_payloads[0].update(file_overrides)
    if not files_without_attribution:
        attribution_hash = (
            sha256(attribution_path.read_bytes()).hexdigest()
            if attribution_path.exists()
            else "0" * 64
        )
        file_payloads.append(
            {
                "relative_path": "ATTRIBUTION.md",
                "role": "attribution",
                "format": "markdown",
                "sha256": attribution_hash,
                "row_count": None,
                "column_names": None,
            }
        )
    payload = _valid_manifest_payload(
        file_hash=file_payloads[0]["sha256"],
        files=file_payloads,
        rat_policy=rat_policy,
    )
    evidence = deepcopy(payload["redistribution_evidence"])
    assert isinstance(evidence, dict)
    if evidence_overrides is not None:
        evidence.update(evidence_overrides)
    if upstream_overrides is not None:
        upstream = evidence["upstream_package"]
        assert isinstance(upstream, dict)
        upstream.update(upstream_overrides)
    if scope_overrides is not None:
        scope = evidence["scope"]
        assert isinstance(scope, dict)
        scope.update(scope_overrides)
    if attribution_overrides is not None:
        attribution = evidence["attribution"]
        assert isinstance(attribution, dict)
        attribution.update(attribution_overrides)
    if remove_evidence_field is not None:
        evidence.pop(remove_evidence_field, None)
    payload["redistribution_evidence"] = evidence
    if manifest_overrides is not None:
        payload.update(manifest_overrides)
    if remove_redistribution_evidence:
        payload.pop("redistribution_evidence", None)
    (bundle_root / "manifest.json").write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )
    return tmp_path / "reference_bundles"


def _valid_manifest_payload(
    *,
    file_hash: object,
    files: list[dict[str, object]],
    rat_policy: bool,
) -> dict[str, object]:
    reference_id = "l6_native" if rat_policy else "unit_reference"
    reference_version = "bundled-snapshot-2026-04-16" if rat_policy else "v1"
    source_version = "PhosR 1.20.0" if rat_policy else "unit-source-v1"
    license_name = "GPL-3 + file LICENSE" if rat_policy else "MIT"
    license_url = (
        "https://github.com/PYangLab/PhosR/blob/master/LICENSE"
        if rat_policy
        else "https://example.test/license"
    )
    return {
        "reference_id": reference_id,
        "display_name": "Rat L6 fixture" if rat_policy else "Unit reference",
        "organism": "Rattus norvegicus",
        "taxonomy_id": 10116,
        "protein_namespace": "display_site_id",
        "reference_version": reference_version,
        "source_name": "PhosR" if rat_policy else "unit source",
        "source_version": source_version,
        "source_url": "https://example.test/reference",
        "retrieved_at": "2026-06-29",
        "table_sha256": file_hash,
        "source_publication": None,
        "license_name": license_name,
        "license_url": license_url,
        "redistribution_status": "approved",
        "redistribution_allowed": True,
        "redistribution_notes": "redistribution approved for exact fixture files",
        "redistribution_evidence": _valid_evidence_payload(
            reference_id=reference_id,
            reference_version=reference_version,
            package_name="PhosR" if rat_policy else "UnitPackage",
            package_version="1.20.0" if rat_policy else "1.0.0",
            license_name=license_name,
            license_url=license_url,
        ),
        "derived_from": ["unit test"],
        "generated_by": "unit test",
        "generated_at_utc": "2026-06-29T00:00:00Z",
        "manifest_schema_version": "1.1",
        "files": files,
        "sequence_context_policy": "centered phosphosite sequence window",
        "sequence_window_length": 3,
        "sequence_center_index": 1,
        "allowed_sequence_alphabet": "ACDEFGHIKLMNPQRSTVWY",
    }


def _valid_evidence_payload(
    *,
    reference_id: str,
    reference_version: str,
    package_name: str = "UnitPackage",
    package_version: str = "1.0.0",
    license_name: str = "MIT",
    license_url: str = "https://example.test/license",
) -> dict[str, object]:
    return {
        "evidence_type": "upstream_package_license",
        "upstream_package": {
            "package_name": package_name,
            "package_version": package_version,
            "license_name": license_name,
            "license_url": license_url,
        },
        "scope": {
            "reference_id": reference_id,
            "reference_version": reference_version,
            "applies_to_exact_packaged_files": True,
            "packaged_files": ["substrate_map.csv", "ATTRIBUTION.md"],
            "applies_to_future_bundles": False,
        },
        "attribution": {
            "repository_notice_path": "NOTICE.md",
            "bundle_attribution_path": "ATTRIBUTION.md",
        },
        "independent_database_permission_claimed": False,
        "evidence_url": "https://example.test/approval-record",
        "verified_at": "2026-06-29",
        "notes": "Approval applies only to this exact packaged snapshot.",
    }
