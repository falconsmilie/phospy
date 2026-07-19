from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[2]
POLICY_PATH = ROOT / "release" / "attestation-policy.json"
BUILDER_PATH = ROOT / "scripts" / "release" / "build_release_attestation.py"
VERIFIER_PATH = ROOT / "scripts" / "release" / "verify_release_attestation.py"
PREPARE_PATH = ROOT / "scripts" / "release" / "prepare_publish_directory.py"
SOURCE_IDENTITY_PATH = ROOT / "scripts" / "release" / "create_source_identity.py"

pytestmark = pytest.mark.release_gate


@dataclass(frozen=True, slots=True)
class EvidenceSet:
    root: Path
    policy: Path
    source_identity: Path
    build_manifest: Path
    wheel: Path
    sdist: Path
    source_reports: tuple[Path, ...]
    artifact_reports: tuple[Path, ...]
    attestation: Path


def _load_module(path: Path, name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _write_json(path: Path, payload: dict[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _prefixed_file_sha256(path: Path) -> str:
    return "sha256:" + _file_sha256(path)


def _complete_evidence(
    tmp_path: Path,
    *,
    source_identity_payload: dict[str, Any] | None = None,
) -> EvidenceSet:
    root = tmp_path
    reports = root / "build" / "reports"
    dist = root / "dist"
    dist.mkdir(parents=True)
    policy = _write_json(
        root / "release" / "attestation-policy.json", _read_json(POLICY_PATH)
    )

    source_identity = root / "build" / "reports" / "source-identity.json"
    _write_json(
        source_identity,
        source_identity_payload
        or {
            "schema": "phospy.source-identity/v1",
            "method": "source-tree",
            "source_tree": {
                "algorithm": "sha256",
                "digest": "sha256:" + ("a" * 64),
                "file_count": 3,
            },
            "package": {"name": "phospy", "version": "1.6.0"},
        },
    )
    source_identity_digest = _prefixed_file_sha256(source_identity)

    wheel = dist / "phospy-1.6.0-py3-none-any.whl"
    sdist = dist / "phospy-1.6.0.tar.gz"
    wheel.write_bytes(b"wheel bytes")
    sdist.write_bytes(b"sdist bytes")
    build_manifest = _write_json(
        reports / "build-manifest.json",
        {
            "schema": "phospy.build-manifest/v1",
            "source_identity_digest": source_identity_digest,
            "package_name": "phospy",
            "package_version": "1.6.0",
            "artifacts": [
                {
                    "kind": "wheel",
                    "filename": wheel.name,
                    "sha256": _file_sha256(wheel),
                },
                {
                    "kind": "sdist",
                    "filename": sdist.name,
                    "sha256": _file_sha256(sdist),
                },
            ],
        },
    )

    policy_payload = _read_json(policy)
    source_reports = tuple(
        _write_source_report(
            reports,
            suite_id=suite_id,
            report_classes=policy_payload["source_suite_report_classes"][suite_id],
            source_identity_digest=source_identity_digest,
        )
        for suite_id in policy_payload["required_source_suite_report_ids"]
    )
    artifact_reports = tuple(
        _write_artifact_report(
            reports,
            artifact_kind=cell["artifact_kind"],
            python_version=cell["python_version"],
            artifact_path=wheel if cell["artifact_kind"] == "wheel" else sdist,
            source_identity_digest=source_identity_digest,
        )
        for cell in policy_payload["required_artifact_verification_matrix"]
    )
    return EvidenceSet(
        root=root,
        policy=policy,
        source_identity=source_identity,
        build_manifest=build_manifest,
        wheel=wheel,
        sdist=sdist,
        source_reports=source_reports,
        artifact_reports=artifact_reports,
        attestation=root / "build" / "release" / "release-attestation.json",
    )


def _write_source_report(
    reports: Path,
    *,
    suite_id: str,
    report_classes: list[str],
    source_identity_digest: str,
    status: str = "success",
) -> Path:
    python_version = suite_id.rsplit("py", maxsplit=1)[1]
    return _write_json(
        reports / f"source-{suite_id}.json",
        {
            "schema": "phospy.source-check/v1",
            "status": status,
            "source_identity_digest": source_identity_digest,
            "package": {"name": "phospy", "version": "1.6.0"},
            "suite": {
                "id": suite_id,
                "report_classes": report_classes,
            },
            "python": python_version,
            "environment": {"python": python_version},
            "evidence": {
                "junit_xml": f"{suite_id}.xml",
                "junit_xml_sha256": "sha256:" + ("1" * 64),
                "summary": {"tests": 1, "failures": 0, "errors": 0, "skipped": 0},
            },
        },
    )


def _write_artifact_report(
    reports: Path,
    *,
    artifact_kind: str,
    python_version: str,
    artifact_path: Path,
    source_identity_digest: str,
    status: str = "success",
    package_version: str = "1.6.0",
    artifact_sha256: str | None = None,
) -> Path:
    check_names = _read_json(POLICY_PATH)["required_verifier_check_names"]
    digest = artifact_sha256 or _file_sha256(artifact_path)
    return _write_json(
        reports / f"{artifact_kind}-artifact-py{python_version}.json",
        {
            "schema": "phospy.artifact-verification/v1",
            "status": status,
            "source_identity_digest": source_identity_digest,
            "artifact": {
                "kind": artifact_kind,
                "filename": artifact_path.name,
                "sha256": digest,
            },
            "package": {"name": "phospy", "version": package_version},
            "environment": {
                "python": f"{python_version}.9",
                "platform": "test",
                "implementation": "CPython",
                "dependency_snapshot_sha256": "sha256:" + ("2" * 64),
            },
            "checks": {name: "pass" for name in check_names},
            "check_details": [
                {
                    "name": name,
                    "status": "pass",
                    "duration_seconds": 0.0,
                    "details": (
                        {
                            "bundle": "rat/l6_native",
                            "manifest_sha256": "3" * 64,
                            "manifest_file_count": 5,
                        }
                        if name == "packaged-scientific-resources"
                        else {}
                    ),
                }
                for name in check_names
            ],
            "distribution": {"name": "phospy", "version": package_version},
        },
    )


def _attest(evidence: EvidenceSet) -> None:
    builder = _load_module(BUILDER_PATH, "build_release_attestation_under_test")
    builder.build_release_attestation(
        policy_path=evidence.policy,
        source_identity_path=evidence.source_identity,
        build_manifest_path=evidence.build_manifest,
        artifact_paths=(evidence.wheel, evidence.sdist),
        source_check_report_paths=evidence.source_reports,
        artifact_verification_report_paths=evidence.artifact_reports,
        output_path=evidence.attestation,
        evidence_root=evidence.root,
        generated_at_utc="2026-07-19T00:00:00Z",
    )


def test_successful_aggregation_using_complete_evidence_set(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)

    _attest(evidence)

    payload = _read_json(evidence.attestation)
    assert payload["schema"] == "phospy.release-attestation/v1"
    assert payload["status"] == "success"
    assert payload["package"] == {"name": "phospy", "version": "1.6.0"}
    assert payload["artifacts"]["wheel"]["sha256"] == _file_sha256(evidence.wheel)
    assert payload["artifacts"]["sdist"]["sha256"] == _file_sha256(evidence.sdist)
    assert len(payload["artifact_verification_matrix"]) == 6
    assert payload["scientific_resources"]["manifest_digests"] == ["3" * 64]


def test_successful_archive_based_aggregation_without_git(tmp_path: Path) -> None:
    source_identity = _load_module(SOURCE_IDENTITY_PATH, "source_identity_under_test")
    archive = tmp_path / "phospy-1.6.0.tar.gz"
    archive.write_bytes(b"source archive")
    source_record = tmp_path / "build" / "reports" / "source-identity.json"
    source_identity.write_source_identity(
        repository_root=tmp_path,
        output_path=source_record,
        source_archive=archive,
        package_name="phospy",
        package_version="1.6.0",
    )
    evidence = _complete_evidence(
        tmp_path, source_identity_payload=_read_json(source_record)
    )

    _attest(evidence)

    payload = _read_json(evidence.attestation)
    assert payload["source_identity"]["method"] == "source-archive"
    assert (
        payload["source_identity"]["source_archive_digest"]
        == _read_json(source_record)["source_archive"]["digest"]
    )


def test_git_source_identity_records_revision_and_dirty_state(tmp_path: Path) -> None:
    source_identity = _load_module(
        SOURCE_IDENTITY_PATH, "source_identity_git_under_test"
    )
    repo = tmp_path / "repo"
    repo.mkdir()
    (repo / "pyproject.toml").write_text(
        '[project]\nname = "phospy"\nversion = "1.6.0"\n',
        encoding="utf-8",
    )
    subprocess.run(("git", "init"), cwd=repo, check=True, capture_output=True)
    subprocess.run(
        ("git", "config", "user.email", "test@example.com"), cwd=repo, check=True
    )
    subprocess.run(("git", "config", "user.name", "Test User"), cwd=repo, check=True)
    subprocess.run(("git", "add", "pyproject.toml"), cwd=repo, check=True)
    subprocess.run(
        ("git", "commit", "-m", "initial"), cwd=repo, check=True, capture_output=True
    )
    revision = subprocess.run(
        ("git", "rev-parse", "HEAD"),
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    (repo / "pyproject.toml").write_text(
        '[project]\nname = "phospy"\nversion = "1.6.1"\n',
        encoding="utf-8",
    )

    output = repo / "build" / "reports" / "source-identity.json"
    source_identity.write_source_identity(repository_root=repo, output_path=output)

    payload = _read_json(output)
    assert payload["method"] == "git"
    assert payload["git"]["revision"] == revision
    assert payload["git"]["dirty"] is True
    assert any("pyproject.toml" in entry for entry in payload["git"]["dirty_paths"])


def test_dirty_git_source_identity_is_rejected_by_public_policy(tmp_path: Path) -> None:
    evidence = _complete_evidence(
        tmp_path,
        source_identity_payload={
            "schema": "phospy.source-identity/v1",
            "method": "git",
            "git": {
                "revision": "a" * 40,
                "dirty": True,
                "dirty_paths": ["M pyproject.toml"],
                "tracked_file_count": 3,
                "untracked_file_count": 0,
            },
            "source_tree": {
                "algorithm": "sha256",
                "digest": "sha256:" + ("a" * 64),
                "file_count": 3,
            },
            "package": {"name": "phospy", "version": "1.6.0"},
        },
    )
    builder = _load_module(BUILDER_PATH, "attestation_dirty_git")

    with pytest.raises(builder.ReleaseAttestationError, match="dirty Git"):
        _attest_with_builder(builder, evidence)


def test_missing_source_suite_report_is_rejected(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    builder = _load_module(BUILDER_PATH, "attestation_missing_source")

    with pytest.raises(
        builder.ReleaseAttestationError, match="missing required source-check"
    ):
        builder.build_release_attestation(
            policy_path=evidence.policy,
            source_identity_path=evidence.source_identity,
            build_manifest_path=evidence.build_manifest,
            artifact_paths=(evidence.wheel, evidence.sdist),
            source_check_report_paths=evidence.source_reports[:-1],
            artifact_verification_report_paths=evidence.artifact_reports,
            output_path=evidence.attestation,
            evidence_root=evidence.root,
        )
    assert not evidence.attestation.exists()


def test_malformed_evidence_report_is_rejected(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    evidence.source_reports[0].write_text("{not-json\n", encoding="utf-8")
    builder = _load_module(BUILDER_PATH, "attestation_malformed_source")

    with pytest.raises(builder.ReleaseAttestationError, match="malformed JSON"):
        _attest_with_builder(builder, evidence)
    assert not evidence.attestation.exists()


def test_missing_wheel_python_matrix_entry_is_rejected(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    builder = _load_module(BUILDER_PATH, "attestation_missing_wheel")
    reports = tuple(
        path
        for path in evidence.artifact_reports
        if path.name != "wheel-artifact-py3.10.json"
    )

    with pytest.raises(builder.ReleaseAttestationError, match="wheel/py3.10"):
        builder.build_release_attestation(
            policy_path=evidence.policy,
            source_identity_path=evidence.source_identity,
            build_manifest_path=evidence.build_manifest,
            artifact_paths=(evidence.wheel, evidence.sdist),
            source_check_report_paths=evidence.source_reports,
            artifact_verification_report_paths=reports,
            output_path=evidence.attestation,
            evidence_root=evidence.root,
        )


def test_missing_sdist_python_matrix_entry_is_rejected(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    builder = _load_module(BUILDER_PATH, "attestation_missing_sdist")
    reports = tuple(
        path
        for path in evidence.artifact_reports
        if path.name != "sdist-artifact-py3.12.json"
    )

    with pytest.raises(builder.ReleaseAttestationError, match="sdist/py3.12"):
        builder.build_release_attestation(
            policy_path=evidence.policy,
            source_identity_path=evidence.source_identity,
            build_manifest_path=evidence.build_manifest,
            artifact_paths=(evidence.wheel, evidence.sdist),
            source_check_report_paths=evidence.source_reports,
            artifact_verification_report_paths=reports,
            output_path=evidence.attestation,
            evidence_root=evidence.root,
        )


def test_failed_artifact_verification_report_is_rejected(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    payload = _read_json(evidence.artifact_reports[0])
    payload["status"] = "failure"
    _write_json(evidence.artifact_reports[0], payload)
    builder = _load_module(BUILDER_PATH, "attestation_failed_artifact")

    with pytest.raises(builder.ReleaseAttestationError, match="did not report success"):
        _attest_with_builder(builder, evidence)


def test_unknown_or_duplicate_report_identity_is_rejected(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    builder = _load_module(BUILDER_PATH, "attestation_duplicate_source")

    with pytest.raises(builder.ReleaseAttestationError, match="duplicate source-check"):
        builder.build_release_attestation(
            policy_path=evidence.policy,
            source_identity_path=evidence.source_identity,
            build_manifest_path=evidence.build_manifest,
            artifact_paths=(evidence.wheel, evidence.sdist),
            source_check_report_paths=(
                *evidence.source_reports,
                evidence.source_reports[0],
            ),
            artifact_verification_report_paths=evidence.artifact_reports,
            output_path=evidence.attestation,
            evidence_root=evidence.root,
        )

    payload = _read_json(evidence.source_reports[0])
    payload["suite"]["id"] = "unknown-suite"
    _write_json(evidence.source_reports[0], payload)
    with pytest.raises(builder.ReleaseAttestationError, match="unknown source-check"):
        _attest_with_builder(builder, evidence)


def test_source_identity_mismatch_is_rejected(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    payload = _read_json(evidence.source_reports[0])
    payload["source_identity_digest"] = "sha256:" + ("4" * 64)
    _write_json(evidence.source_reports[0], payload)
    builder = _load_module(BUILDER_PATH, "attestation_source_mismatch")

    with pytest.raises(
        builder.ReleaseAttestationError, match="source identity mismatch"
    ):
        _attest_with_builder(builder, evidence)


def test_package_version_mismatch_is_rejected(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    payload = _read_json(evidence.artifact_reports[0])
    payload["package"]["version"] = "9.9.9"
    payload["distribution"]["version"] = "9.9.9"
    _write_json(evidence.artifact_reports[0], payload)
    builder = _load_module(BUILDER_PATH, "attestation_package_mismatch")

    with pytest.raises(
        builder.ReleaseAttestationError, match="package identity mismatch"
    ):
        _attest_with_builder(builder, evidence)


def test_build_manifest_artifact_digest_mismatch_is_rejected(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    payload = _read_json(evidence.build_manifest)
    payload["artifacts"][0]["sha256"] = "5" * 64
    _write_json(evidence.build_manifest, payload)
    builder = _load_module(BUILDER_PATH, "attestation_manifest_digest")

    with pytest.raises(builder.ReleaseAttestationError, match="digest mismatch"):
        _attest_with_builder(builder, evidence)


def test_artifact_verification_digest_mismatch_is_rejected(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    payload = _read_json(evidence.artifact_reports[0])
    payload["artifact"]["sha256"] = "6" * 64
    _write_json(evidence.artifact_reports[0], payload)
    builder = _load_module(BUILDER_PATH, "attestation_artifact_digest")

    with pytest.raises(
        builder.ReleaseAttestationError, match="digest does not match actual"
    ):
        _attest_with_builder(builder, evidence)


def test_wheel_tampering_is_detected(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    evidence.wheel.write_bytes(b"tampered wheel")
    builder = _load_module(BUILDER_PATH, "attestation_wheel_tamper")

    with pytest.raises(builder.ReleaseAttestationError, match="wheel digest mismatch"):
        _attest_with_builder(builder, evidence)


def test_sdist_tampering_is_detected(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    evidence.sdist.write_bytes(b"tampered sdist")
    builder = _load_module(BUILDER_PATH, "attestation_sdist_tamper")

    with pytest.raises(builder.ReleaseAttestationError, match="sdist digest mismatch"):
        _attest_with_builder(builder, evidence)


def test_source_check_report_tampering_invalidates_verification(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    _attest(evidence)
    payload = _read_json(evidence.source_reports[0])
    payload["evidence"]["summary"]["tests"] = 2
    _write_json(evidence.source_reports[0], payload)
    verifier = _load_module(VERIFIER_PATH, "verify_attestation_source_tamper")

    with pytest.raises(
        verifier.ReleaseAttestationVerificationError, match="digest changed"
    ):
        verifier.verify_release_attestation(
            attestation_path=evidence.attestation,
            evidence_root=evidence.root,
            artifact_dir=evidence.wheel.parent,
        )


def test_scientific_report_tampering_invalidates_verification(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    _attest(evidence)
    scientific_report = next(
        path
        for path in evidence.source_reports
        if "reference-manifest-py3.10" in path.name
    )
    payload = _read_json(scientific_report)
    payload["evidence"]["summary"]["tests"] = 99
    _write_json(scientific_report, payload)
    verifier = _load_module(VERIFIER_PATH, "verify_attestation_science_tamper")

    with pytest.raises(
        verifier.ReleaseAttestationVerificationError, match="digest changed"
    ):
        verifier.verify_release_attestation(
            attestation_path=evidence.attestation,
            evidence_root=evidence.root,
            artifact_dir=evidence.wheel.parent,
        )


def test_stale_final_attestation_is_removed_before_validation(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    evidence.attestation.parent.mkdir(parents=True)
    evidence.attestation.write_text('{"status":"success"}\n', encoding="utf-8")
    builder = _load_module(BUILDER_PATH, "attestation_stale_removal")

    with pytest.raises(builder.ReleaseAttestationError):
        builder.build_release_attestation(
            policy_path=evidence.policy,
            source_identity_path=evidence.source_identity,
            build_manifest_path=evidence.build_manifest,
            artifact_paths=(evidence.wheel, evidence.sdist),
            source_check_report_paths=evidence.source_reports[:-1],
            artifact_verification_report_paths=evidence.artifact_reports,
            output_path=evidence.attestation,
            evidence_root=evidence.root,
        )
    assert not evidence.attestation.exists()


def test_final_file_is_written_only_after_all_validation_succeeds(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    evidence = _complete_evidence(tmp_path)
    builder = _load_module(BUILDER_PATH, "attestation_atomic_after_validation")
    called = False

    def fail_if_called(path: Path, payload: dict[str, Any]) -> None:
        nonlocal called
        called = True
        raise AssertionError("atomic writer should not run")

    monkeypatch.setattr(builder, "_atomic_write_json", fail_if_called)

    with pytest.raises(builder.ReleaseAttestationError):
        builder.build_release_attestation(
            policy_path=evidence.policy,
            source_identity_path=evidence.source_identity,
            build_manifest_path=evidence.build_manifest,
            artifact_paths=(evidence.wheel, evidence.sdist),
            source_check_report_paths=evidence.source_reports[:-1],
            artifact_verification_report_paths=evidence.artifact_reports,
            output_path=evidence.attestation,
            evidence_root=evidence.root,
        )
    assert called is False
    assert not evidence.attestation.exists()


def test_publish_preparation_rejects_unlisted_artifact(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    _attest(evidence)
    (evidence.wheel.parent / "phospy-1.6.0-extra-py3-none-any.whl").write_bytes(
        b"extra"
    )
    prepare = _load_module(PREPARE_PATH, "prepare_publish_extra")

    with pytest.raises(prepare.PublishPreparationError, match="unattested"):
        prepare.prepare_publish_directory(
            attestation_path=evidence.attestation,
            dist_dir=evidence.wheel.parent,
            output_dir=evidence.root / "build" / "publish",
        )


def test_publish_preparation_rejects_changed_listed_artifact(tmp_path: Path) -> None:
    evidence = _complete_evidence(tmp_path)
    _attest(evidence)
    evidence.wheel.write_bytes(b"changed")
    prepare = _load_module(PREPARE_PATH, "prepare_publish_changed")

    with pytest.raises(prepare.PublishPreparationError, match="digest changed"):
        prepare.prepare_publish_directory(
            attestation_path=evidence.attestation,
            dist_dir=evidence.wheel.parent,
            output_dir=evidence.root / "build" / "publish",
        )


def test_publish_preparation_copies_exactly_attested_wheel_and_sdist(
    tmp_path: Path,
) -> None:
    evidence = _complete_evidence(tmp_path)
    _attest(evidence)
    prepare = _load_module(PREPARE_PATH, "prepare_publish_success")
    output = evidence.root / "build" / "publish"

    prepare.prepare_publish_directory(
        attestation_path=evidence.attestation,
        dist_dir=evidence.wheel.parent,
        output_dir=output,
    )

    assert sorted(path.name for path in output.iterdir()) == sorted(
        [evidence.sdist.name, evidence.wheel.name]
    )


def test_attestation_builder_does_not_run_pytest_or_import_runtime_package() -> None:
    source = BUILDER_PATH.read_text(encoding="utf-8")

    assert "pytest" not in source.lower()
    assert "import phospy" not in source
    assert "from phospy" not in source


def _attest_with_builder(builder: ModuleType, evidence: EvidenceSet) -> None:
    builder.build_release_attestation(
        policy_path=evidence.policy,
        source_identity_path=evidence.source_identity,
        build_manifest_path=evidence.build_manifest,
        artifact_paths=(evidence.wheel, evidence.sdist),
        source_check_report_paths=evidence.source_reports,
        artifact_verification_report_paths=evidence.artifact_reports,
        output_path=evidence.attestation,
        evidence_root=evidence.root,
        generated_at_utc="2026-07-19T00:00:00Z",
    )
