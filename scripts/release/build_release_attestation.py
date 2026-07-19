from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ATTESTATION_SCHEMA = "phospy.release-attestation/v1"
POLICY_SCHEMA = "phospy.attestation-policy/v1"
SOURCE_IDENTITY_SCHEMA = "phospy.source-identity/v1"
SOURCE_CHECK_REPORT_SCHEMA = "phospy.source-check/v1"
ARTIFACT_VERIFICATION_SCHEMA = "phospy.artifact-verification/v1"
BUILD_MANIFEST_SCHEMA = "phospy.build-manifest/v1"

_ARTIFACT_SUFFIXES = {
    "wheel": ".whl",
    "sdist": ".tar.gz",
}


class ReleaseAttestationError(RuntimeError):
    """Raised when release evidence cannot support a success attestation."""


@dataclass(frozen=True, slots=True)
class _JsonEvidence:
    path: Path
    filename: str
    digest: str
    payload: dict[str, Any]


@dataclass(frozen=True, slots=True)
class _ArtifactEvidence:
    kind: str
    path: Path
    filename: str
    sha256: str


@dataclass(frozen=True, slots=True)
class _Policy:
    path: Path
    digest: str
    payload: dict[str, Any]
    supported_python_versions: tuple[str, ...]
    required_artifact_kinds: tuple[str, ...]
    required_matrix: tuple[tuple[str, str], ...]
    required_verifier_check_names: tuple[str, ...]
    required_source_suite_report_ids: tuple[str, ...]
    source_suite_report_classes: dict[str, tuple[str, ...]]
    required_report_classes: tuple[str, ...]
    allow_dirty_git: bool


def build_release_attestation(
    *,
    policy_path: Path,
    source_identity_path: Path,
    build_manifest_path: Path,
    artifact_paths: Sequence[Path],
    source_check_report_paths: Sequence[Path],
    artifact_verification_report_paths: Sequence[Path],
    output_path: Path,
    evidence_root: Path | None = None,
    workflow_run_id: str | None = None,
    generated_at_utc: str | None = None,
) -> Path:
    _remove_stale_output(output_path)
    attestation = _build_attestation_payload(
        policy_path=policy_path,
        source_identity_path=source_identity_path,
        build_manifest_path=build_manifest_path,
        artifact_paths=artifact_paths,
        source_check_report_paths=source_check_report_paths,
        artifact_verification_report_paths=artifact_verification_report_paths,
        evidence_root=(evidence_root or Path.cwd()).resolve(),
        workflow_run_id=workflow_run_id,
        generated_at_utc=generated_at_utc,
    )
    _atomic_write_json(output_path, attestation)
    return output_path


def _build_attestation_payload(
    *,
    policy_path: Path,
    source_identity_path: Path,
    build_manifest_path: Path,
    artifact_paths: Sequence[Path],
    source_check_report_paths: Sequence[Path],
    artifact_verification_report_paths: Sequence[Path],
    evidence_root: Path,
    workflow_run_id: str | None,
    generated_at_utc: str | None,
) -> dict[str, Any]:
    policy = _load_policy(policy_path, evidence_root=evidence_root)
    source_identity = _load_json_evidence(
        source_identity_path, evidence_root=evidence_root
    )
    _require_schema(
        source_identity.payload,
        SOURCE_IDENTITY_SCHEMA,
        label="source identity record",
    )
    _validate_source_identity(source_identity, policy)

    build_manifest = _load_json_evidence(
        build_manifest_path, evidence_root=evidence_root
    )
    package = _validate_build_manifest(
        build_manifest=build_manifest,
        source_identity_digest=source_identity.digest,
    )
    _validate_source_package(source_identity, package)
    artifacts = _validate_artifacts(
        artifact_paths=artifact_paths,
        build_manifest=build_manifest.payload,
        policy=policy,
    )
    source_checks = _validate_source_checks(
        report_paths=source_check_report_paths,
        policy=policy,
        source_identity_digest=source_identity.digest,
        package=package,
        evidence_root=evidence_root,
    )
    artifact_matrix = _validate_artifact_reports(
        report_paths=artifact_verification_report_paths,
        policy=policy,
        source_identity_digest=source_identity.digest,
        package=package,
        artifacts=artifacts,
        evidence_root=evidence_root,
    )
    report_class_coverage = _report_class_coverage(source_checks, policy)

    return {
        "schema": ATTESTATION_SCHEMA,
        "status": "success",
        "package": package,
        "generated_at_utc": generated_at_utc or _utc_timestamp(),
        "source_identity": _source_identity_summary(source_identity),
        "release_evidence_policy": {
            "filename": policy.payload["_filename"],
            "digest": policy.digest,
        },
        "build_manifest": {
            "filename": build_manifest.filename,
            "digest": build_manifest.digest,
        },
        "source_checks": source_checks,
        "report_class_coverage": report_class_coverage,
        "artifacts": {
            artifact.kind: {
                "filename": artifact.filename,
                "sha256": artifact.sha256,
            }
            for artifact in sorted(artifacts.values(), key=lambda item: item.kind)
        },
        "artifact_verification_matrix": artifact_matrix,
        "scientific_resources": _scientific_resource_summary(artifact_matrix),
        "environment_reports": _environment_summary(source_checks, artifact_matrix),
        **(
            {}
            if workflow_run_id is None
            else {"workflow_run_id": workflow_run_id.strip()}
        ),
    }


def _load_policy(path: Path, *, evidence_root: Path) -> _Policy:
    evidence = _load_json_evidence(path, evidence_root=evidence_root)
    payload = evidence.payload
    _require_schema(payload, POLICY_SCHEMA, label="attestation policy")
    attestation_schema = _required_text(
        payload.get("attestation_schema"),
        field_name="attestation policy attestation_schema",
    )
    _require(
        attestation_schema == ATTESTATION_SCHEMA,
        "attestation policy schema target mismatch: "
        f"expected {ATTESTATION_SCHEMA!r}, got {attestation_schema!r}",
    )
    supported_python_versions = tuple(
        _unique_text_list(
            payload.get("supported_python_versions"),
            field_name="supported_python_versions",
        )
    )
    required_artifact_kinds = tuple(
        _unique_text_list(
            payload.get("required_artifact_kinds"),
            field_name="required_artifact_kinds",
        )
    )
    for kind in required_artifact_kinds:
        _require(
            kind in _ARTIFACT_SUFFIXES, f"unsupported policy artifact kind: {kind!r}"
        )

    required_matrix = _required_matrix(payload)
    expected_matrix = {
        (kind, python_version)
        for kind in required_artifact_kinds
        for python_version in supported_python_versions
    }
    _require(
        set(required_matrix) == expected_matrix,
        "attestation policy artifact-verification matrix must be the full "
        "artifact kind and supported Python cross-product",
    )

    source_suite_report_ids = tuple(
        _unique_text_list(
            payload.get("required_source_suite_report_ids"),
            field_name="required_source_suite_report_ids",
        )
    )
    source_suite_report_classes = _source_suite_report_classes(payload)
    _require(
        set(source_suite_report_ids) == set(source_suite_report_classes),
        "source suite report class mapping must cover exactly the required source "
        "suite report IDs",
    )
    required_report_classes = tuple(
        _unique_text_list(
            payload.get("required_report_classes"),
            field_name="required_report_classes",
        )
    )
    for suite_id, classes in source_suite_report_classes.items():
        unknown_classes = sorted(set(classes) - set(required_report_classes))
        _require(
            unknown_classes == [],
            f"source suite {suite_id!r} declares unknown report classes: "
            + ", ".join(unknown_classes),
        )

    payload["_filename"] = evidence.filename
    return _Policy(
        path=evidence.path,
        digest=evidence.digest,
        payload=payload,
        supported_python_versions=supported_python_versions,
        required_artifact_kinds=required_artifact_kinds,
        required_matrix=required_matrix,
        required_verifier_check_names=tuple(
            _unique_text_list(
                payload.get("required_verifier_check_names"),
                field_name="required_verifier_check_names",
            )
        ),
        required_source_suite_report_ids=source_suite_report_ids,
        source_suite_report_classes=source_suite_report_classes,
        required_report_classes=required_report_classes,
        allow_dirty_git=bool(payload.get("allow_dirty_git", False)),
    )


def _required_matrix(payload: Mapping[str, object]) -> tuple[tuple[str, str], ...]:
    matrix = payload.get("required_artifact_verification_matrix")
    _require(
        isinstance(matrix, list),
        "required_artifact_verification_matrix must be an array",
    )
    parsed: list[tuple[str, str]] = []
    for position, entry in enumerate(matrix):
        _require(
            isinstance(entry, Mapping),
            f"required_artifact_verification_matrix[{position}] must be an object",
        )
        kind = _required_text(
            entry.get("artifact_kind"),
            field_name=f"required_artifact_verification_matrix[{position}].artifact_kind",
        )
        python_version = _required_text(
            entry.get("python_version"),
            field_name=f"required_artifact_verification_matrix[{position}].python_version",
        )
        parsed.append((kind, python_version))
    _require(
        len(parsed) == len(set(parsed)),
        "required_artifact_verification_matrix contains duplicate cells",
    )
    return tuple(sorted(parsed))


def _source_suite_report_classes(
    payload: Mapping[str, object],
) -> dict[str, tuple[str, ...]]:
    value = payload.get("source_suite_report_classes")
    _require(
        isinstance(value, Mapping), "source_suite_report_classes must be an object"
    )
    result: dict[str, tuple[str, ...]] = {}
    for raw_suite_id, raw_classes in value.items():
        suite_id = _required_text(raw_suite_id, field_name="source suite report id")
        result[suite_id] = tuple(
            _unique_text_list(
                raw_classes,
                field_name=f"source_suite_report_classes[{suite_id}]",
            )
        )
    return result


def _validate_source_identity(source_identity: _JsonEvidence, policy: _Policy) -> None:
    payload = source_identity.payload
    method = _required_text(payload.get("method"), field_name="source identity method")
    if method == "git":
        git = _required_mapping(payload.get("git"), field_name="source identity git")
        dirty = git.get("dirty")
        _require(isinstance(dirty, bool), "source identity git.dirty must be boolean")
        _required_text(git.get("revision"), field_name="source identity git.revision")
        if dirty and not policy.allow_dirty_git:
            raise ReleaseAttestationError(
                "dirty Git source identity cannot be attested by this public "
                "release policy"
            )
        tree = _required_mapping(
            payload.get("source_tree"),
            field_name="source identity source_tree",
        )
        _require_sha256_digest(
            _required_text(tree.get("digest"), field_name="source tree digest"),
            field_name="source tree digest",
        )
    elif method == "source-tree":
        tree = _required_mapping(
            payload.get("source_tree"),
            field_name="source identity source_tree",
        )
        _require_sha256_digest(
            _required_text(tree.get("digest"), field_name="source tree digest"),
            field_name="source tree digest",
        )
    elif method == "source-archive":
        archive = _required_mapping(
            payload.get("source_archive"),
            field_name="source identity source_archive",
        )
        _require_sha256_digest(
            _required_text(archive.get("digest"), field_name="source archive digest"),
            field_name="source archive digest",
        )
    else:
        raise ReleaseAttestationError(f"unsupported source identity method: {method!r}")


def _validate_build_manifest(
    *,
    build_manifest: _JsonEvidence,
    source_identity_digest: str,
) -> dict[str, str]:
    payload = build_manifest.payload
    _require_schema(payload, BUILD_MANIFEST_SCHEMA, label="build manifest")
    observed_source_digest = _required_text(
        payload.get("source_identity_digest"),
        field_name="build manifest source_identity_digest",
    )
    _require(
        observed_source_digest == source_identity_digest,
        "build manifest source identity mismatch",
    )
    name = _required_text(
        payload.get("package_name", "phospy"),
        field_name="build manifest package_name",
    )
    version = _required_text(
        payload.get("package_version"),
        field_name="build manifest package_version",
    )
    return {"name": name, "version": version}


def _validate_source_package(
    source_identity: _JsonEvidence,
    package: Mapping[str, str],
) -> None:
    source_package = source_identity.payload.get("package")
    _require(source_package is not None, "source identity package metadata is required")
    _require_package(
        source_package,
        expected=package,
        label=source_identity.filename,
    )


def _validate_artifacts(
    *,
    artifact_paths: Sequence[Path],
    build_manifest: Mapping[str, object],
    policy: _Policy,
) -> dict[str, _ArtifactEvidence]:
    artifacts = build_manifest.get("artifacts")
    _require(isinstance(artifacts, list), "build manifest artifacts must be an array")
    manifest_artifacts: dict[str, Mapping[str, object]] = {}
    for position, entry in enumerate(artifacts):
        _require(
            isinstance(entry, Mapping),
            f"build manifest artifacts[{position}] must be an object",
        )
        kind = _required_text(
            entry.get("kind"),
            field_name=f"build manifest artifacts[{position}].kind",
        )
        _require(
            kind in policy.required_artifact_kinds,
            f"unknown artifact kind in build manifest: {kind!r}",
        )
        _require(
            kind not in manifest_artifacts,
            f"duplicate artifact kind in build manifest: {kind!r}",
        )
        manifest_artifacts[kind] = entry
    _require(
        set(manifest_artifacts) == set(policy.required_artifact_kinds),
        "build manifest must contain exactly the policy-required artifact kinds",
    )

    supplied_artifacts: dict[str, _ArtifactEvidence] = {}
    for raw_path in artifact_paths:
        path = raw_path.resolve()
        _require(path.is_file(), f"artifact file is missing: {path.as_posix()}")
        kind = _artifact_kind_from_filename(path.name)
        _require(
            kind in policy.required_artifact_kinds,
            f"unknown supplied artifact kind: {kind!r}",
        )
        _require(
            kind not in supplied_artifacts,
            f"duplicate supplied artifact kind: {kind!r}",
        )
        supplied_artifacts[kind] = _ArtifactEvidence(
            kind=kind,
            path=path,
            filename=path.name,
            sha256=_file_sha256(path),
        )
    _require(
        set(supplied_artifacts) == set(policy.required_artifact_kinds),
        "actual artifact files must contain exactly one required wheel and one "
        "required sdist",
    )

    for kind, artifact in supplied_artifacts.items():
        manifest_entry = manifest_artifacts[kind]
        manifest_filename = _required_text(
            manifest_entry.get("filename"),
            field_name=f"build manifest {kind}.filename",
        )
        manifest_sha256 = _required_text(
            manifest_entry.get("sha256"),
            field_name=f"build manifest {kind}.sha256",
        )
        _require_sha256_hex(manifest_sha256, field_name=f"build manifest {kind}.sha256")
        _require(
            manifest_filename == artifact.filename,
            f"{kind} filename mismatch between build manifest and artifact file",
        )
        _require(
            manifest_sha256 == artifact.sha256,
            f"{kind} digest mismatch between build manifest and artifact file",
        )
    return supplied_artifacts


def _validate_source_checks(
    *,
    report_paths: Sequence[Path],
    policy: _Policy,
    source_identity_digest: str,
    package: Mapping[str, str],
    evidence_root: Path,
) -> list[dict[str, Any]]:
    seen: dict[str, _JsonEvidence] = {}
    records: list[dict[str, Any]] = []
    for report_path in report_paths:
        evidence = _load_json_evidence(report_path, evidence_root=evidence_root)
        payload = evidence.payload
        _require_schema(
            payload, SOURCE_CHECK_REPORT_SCHEMA, label="source-check report"
        )
        _require_status_success(
            payload, label=f"source-check report {evidence.filename}"
        )
        _require(
            _required_text(
                payload.get("source_identity_digest"),
                field_name="source-check source_identity_digest",
            )
            == source_identity_digest,
            f"source-check report {evidence.filename} source identity mismatch",
        )
        _require_package(
            payload.get("package"), expected=package, label=evidence.filename
        )
        suite = _required_mapping(payload.get("suite"), field_name="source-check suite")
        suite_id = _required_text(suite.get("id"), field_name="source-check suite.id")
        _require(
            suite_id in policy.required_source_suite_report_ids,
            f"unknown source-check report identity: {suite_id!r}",
        )
        _require(
            suite_id not in seen,
            f"duplicate source-check report identity: {suite_id!r}",
        )
        report_classes = tuple(
            _unique_text_list(
                suite.get("report_classes"),
                field_name=f"source-check {suite_id} report_classes",
            )
        )
        expected_classes = policy.source_suite_report_classes[suite_id]
        _require(
            set(report_classes) == set(expected_classes),
            f"source-check report {suite_id!r} report classes do not match policy",
        )
        seen[suite_id] = evidence
        records.append(
            {
                "id": suite_id,
                "status": "success",
                "report_classes": list(report_classes),
                "python": _optional_text(payload.get("python")),
                "filename": evidence.filename,
                "digest": evidence.digest,
                "environment": _optional_mapping(payload.get("environment")),
                "evidence": _optional_mapping(payload.get("evidence")),
            }
        )

    missing = sorted(set(policy.required_source_suite_report_ids) - set(seen))
    _require(
        missing == [], "missing required source-check reports: " + ", ".join(missing)
    )
    return sorted(records, key=lambda item: str(item["id"]))


def _validate_artifact_reports(
    *,
    report_paths: Sequence[Path],
    policy: _Policy,
    source_identity_digest: str,
    package: Mapping[str, str],
    artifacts: Mapping[str, _ArtifactEvidence],
    evidence_root: Path,
) -> list[dict[str, Any]]:
    seen: dict[tuple[str, str], _JsonEvidence] = {}
    records: list[dict[str, Any]] = []
    for report_path in report_paths:
        evidence = _load_json_evidence(report_path, evidence_root=evidence_root)
        payload = evidence.payload
        _require_schema(
            payload,
            ARTIFACT_VERIFICATION_SCHEMA,
            label="artifact-verification report",
        )
        _require_status_success(
            payload,
            label=f"artifact-verification report {evidence.filename}",
        )
        _require(
            _required_text(
                payload.get("source_identity_digest"),
                field_name="artifact-verification source_identity_digest",
            )
            == source_identity_digest,
            f"artifact-verification report {evidence.filename} source identity mismatch",
        )
        _require_package(
            payload.get("package"), expected=package, label=evidence.filename
        )
        _require_distribution_version(
            payload, expected=package, label=evidence.filename
        )

        artifact_payload = _required_mapping(
            payload.get("artifact"),
            field_name="artifact-verification artifact",
        )
        kind = _required_text(
            artifact_payload.get("kind"),
            field_name="artifact-verification artifact.kind",
        )
        filename = _required_text(
            artifact_payload.get("filename"),
            field_name="artifact-verification artifact.filename",
        )
        artifact_sha256 = _required_text(
            artifact_payload.get("sha256"),
            field_name="artifact-verification artifact.sha256",
        )
        _require_sha256_hex(
            artifact_sha256,
            field_name="artifact-verification artifact.sha256",
        )
        _require(
            kind in artifacts,
            f"artifact report references unknown artifact kind: {kind!r}",
        )
        artifact = artifacts[kind]
        _require(
            filename == artifact.filename,
            f"artifact report {evidence.filename} filename does not match attested {kind}",
        )
        _require(
            artifact_sha256 == artifact.sha256,
            f"artifact report {evidence.filename} digest does not match actual {kind}",
        )

        environment = _required_mapping(
            payload.get("environment"),
            field_name="artifact-verification environment",
        )
        python_version = _python_minor_version(
            _required_text(
                environment.get("python"),
                field_name="artifact-verification environment.python",
            )
        )
        key = (kind, python_version)
        _require(
            key in policy.required_matrix,
            f"unknown artifact-verification matrix cell: {key!r}",
        )
        _require(
            key not in seen, f"duplicate artifact-verification matrix cell: {key!r}"
        )
        seen[key] = evidence
        check_statuses = _validate_required_verifier_checks(payload, policy)
        scientific_resources = _scientific_resources_from_report(payload)
        _require(
            scientific_resources is not None,
            f"artifact report {evidence.filename} lacks scientific resource summary",
        )
        records.append(
            {
                "artifact_kind": kind,
                "python_version": python_version,
                "status": "success",
                "filename": evidence.filename,
                "digest": evidence.digest,
                "artifact_filename": artifact.filename,
                "artifact_sha256": artifact.sha256,
                "checks": check_statuses,
                "environment": {
                    "python": environment.get("python"),
                    "platform": environment.get("platform"),
                    "implementation": environment.get("implementation"),
                    "dependency_snapshot_sha256": environment.get(
                        "dependency_snapshot_sha256"
                    ),
                },
                "scientific_resources": scientific_resources,
            }
        )

    missing = sorted(set(policy.required_matrix) - set(seen))
    _require(
        missing == [],
        "missing required artifact-verification reports: "
        + ", ".join(f"{kind}/py{python}" for kind, python in missing),
    )
    return sorted(
        records,
        key=lambda item: (str(item["artifact_kind"]), str(item["python_version"])),
    )


def _validate_required_verifier_checks(
    payload: Mapping[str, object],
    policy: _Policy,
) -> dict[str, str]:
    details = payload.get("check_details")
    _require(
        isinstance(details, list),
        "artifact-verification check_details must be an array",
    )
    observed: dict[str, str] = {}
    for position, entry in enumerate(details):
        _require(
            isinstance(entry, Mapping),
            f"artifact-verification check_details[{position}] must be an object",
        )
        name = _required_text(
            entry.get("name"),
            field_name=f"artifact-verification check_details[{position}].name",
        )
        status = _required_text(
            entry.get("status"),
            field_name=f"artifact-verification check_details[{position}].status",
        )
        _require(
            name in policy.required_verifier_check_names,
            f"unknown artifact verifier check name: {name!r}",
        )
        _require(
            name not in observed, f"duplicate artifact verifier check name: {name!r}"
        )
        _require(status == "pass", f"artifact verifier check did not pass: {name}")
        observed[name] = status
    missing = sorted(set(policy.required_verifier_check_names) - set(observed))
    _require(missing == [], "missing artifact verifier checks: " + ", ".join(missing))
    return {name: "pass" for name in policy.required_verifier_check_names}


def _scientific_resources_from_report(
    payload: Mapping[str, object],
) -> dict[str, Any] | None:
    details = payload.get("check_details")
    if not isinstance(details, list):
        return None
    for entry in details:
        if not isinstance(entry, Mapping):
            continue
        if entry.get("name") != "packaged-scientific-resources":
            continue
        detail_payload = entry.get("details")
        if not isinstance(detail_payload, Mapping):
            return None
        summary = dict(detail_payload)
        manifest_digest = summary.get("manifest_sha256")
        if isinstance(manifest_digest, str):
            _require_sha256_hex(
                manifest_digest,
                field_name="scientific resource manifest_sha256",
            )
        else:
            raise ReleaseAttestationError(
                "packaged-scientific-resources check must report manifest_sha256"
            )
        return summary
    return None


def _report_class_coverage(
    source_checks: Sequence[Mapping[str, object]],
    policy: _Policy,
) -> dict[str, list[str]]:
    coverage: dict[str, list[str]] = {key: [] for key in policy.required_report_classes}
    for record in source_checks:
        report_id = str(record["id"])
        classes = record.get("report_classes")
        if not isinstance(classes, list):
            continue
        for class_name in classes:
            if isinstance(class_name, str) and class_name in coverage:
                coverage[class_name].append(report_id)
    missing = sorted(key for key, report_ids in coverage.items() if not report_ids)
    _require(missing == [], "missing required report classes: " + ", ".join(missing))
    return {key: sorted(value) for key, value in coverage.items()}


def _source_identity_summary(source_identity: _JsonEvidence) -> dict[str, Any]:
    payload = source_identity.payload
    method = _required_text(payload.get("method"), field_name="source identity method")
    summary: dict[str, Any] = {
        "filename": source_identity.filename,
        "digest": source_identity.digest,
        "method": method,
    }
    if method in {"git", "source-tree"}:
        tree = _required_mapping(
            payload.get("source_tree"),
            field_name="source identity source_tree",
        )
        summary["source_tree_digest"] = tree["digest"]
    if method == "source-archive":
        archive = _required_mapping(
            payload.get("source_archive"),
            field_name="source identity source_archive",
        )
        summary["source_archive_digest"] = archive["digest"]
    if method == "git":
        git = _required_mapping(payload.get("git"), field_name="source identity git")
        summary["git_revision"] = git["revision"]
        summary["git_dirty"] = git["dirty"]
    return summary


def _scientific_resource_summary(
    artifact_matrix: Sequence[Mapping[str, object]],
) -> dict[str, Any]:
    manifest_digests = sorted(
        {
            str(resources["manifest_sha256"])
            for item in artifact_matrix
            if isinstance(item.get("scientific_resources"), Mapping)
            for resources in [item["scientific_resources"]]
            if isinstance(resources.get("manifest_sha256"), str)
        }
    )
    summaries = [
        {
            "artifact_kind": item["artifact_kind"],
            "python_version": item["python_version"],
            "summary": item.get("scientific_resources"),
        }
        for item in artifact_matrix
        if item.get("scientific_resources") is not None
    ]
    return {
        "manifest_digests": manifest_digests,
        "verification_summaries": summaries,
    }


def _environment_summary(
    source_checks: Sequence[Mapping[str, object]],
    artifact_matrix: Sequence[Mapping[str, object]],
) -> dict[str, Any]:
    source = [
        {
            "id": item["id"],
            "python": item.get("python"),
            "environment": item.get("environment"),
        }
        for item in source_checks
        if item.get("environment") is not None or item.get("python") is not None
    ]
    artifacts = [
        {
            "artifact_kind": item["artifact_kind"],
            "python_version": item["python_version"],
            "environment": item.get("environment"),
        }
        for item in artifact_matrix
    ]
    return {"source_checks": source, "artifact_verification": artifacts}


def _remove_stale_output(output_path: Path) -> None:
    if output_path.exists():
        _require(
            output_path.is_file() or output_path.is_symlink(),
            f"final attestation output exists and is not a file: {output_path.as_posix()}",
        )
        output_path.unlink()


def _atomic_write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        "w",
        delete=False,
        dir=path.parent,
        encoding="utf-8",
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    temp_path = Path(handle.name)
    try:
        with handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        temp_path.replace(path)
    except BaseException:
        try:
            temp_path.unlink()
        except FileNotFoundError:
            pass
        raise


def _load_json_evidence(path: Path, *, evidence_root: Path) -> _JsonEvidence:
    resolved = path.resolve()
    _require(resolved.is_file(), f"evidence file is missing: {resolved.as_posix()}")
    try:
        payload = json.loads(resolved.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ReleaseAttestationError(
            f"evidence file is malformed JSON: {resolved.as_posix()}: {exc}"
        ) from exc
    _require(
        isinstance(payload, dict),
        f"evidence file must contain a JSON object: {resolved.as_posix()}",
    )
    return _JsonEvidence(
        path=resolved,
        filename=_display_path(resolved, evidence_root),
        digest=_prefixed_file_sha256(resolved),
        payload=payload,
    )


def _display_path(path: Path, evidence_root: Path) -> str:
    try:
        return path.resolve().relative_to(evidence_root.resolve()).as_posix()
    except ValueError:
        return path.name


def _artifact_kind_from_filename(filename: str) -> str:
    if filename.endswith(_ARTIFACT_SUFFIXES["wheel"]):
        return "wheel"
    if filename.endswith(_ARTIFACT_SUFFIXES["sdist"]):
        return "sdist"
    raise ReleaseAttestationError(f"unsupported artifact filename: {filename!r}")


def _require_schema(
    payload: Mapping[str, object], expected: str, *, label: str
) -> None:
    observed = _required_text(payload.get("schema"), field_name=f"{label} schema")
    _require(
        observed == expected,
        f"{label} schema mismatch: expected {expected!r}, got {observed!r}",
    )


def _require_status_success(payload: Mapping[str, object], *, label: str) -> None:
    status = _required_text(payload.get("status"), field_name=f"{label} status")
    _require(status == "success", f"{label} did not report success: {status!r}")


def _require_package(
    value: object,
    *,
    expected: Mapping[str, str],
    label: str,
) -> None:
    package = _required_mapping(value, field_name=f"{label} package")
    observed = {
        "name": _required_text(package.get("name"), field_name=f"{label} package.name"),
        "version": _required_text(
            package.get("version"),
            field_name=f"{label} package.version",
        ),
    }
    _require(observed == dict(expected), f"{label} package identity mismatch")


def _require_distribution_version(
    payload: Mapping[str, object],
    *,
    expected: Mapping[str, str],
    label: str,
) -> None:
    distribution = payload.get("distribution")
    if distribution is None:
        return
    distribution_mapping = _required_mapping(
        distribution,
        field_name=f"{label} distribution",
    )
    name = distribution_mapping.get("name")
    if name is not None:
        _require(
            _required_text(name, field_name=f"{label} distribution.name")
            == expected["name"],
            f"{label} installed distribution name mismatch",
        )
    version = distribution_mapping.get("version")
    if version is not None:
        _require(
            _required_text(version, field_name=f"{label} distribution.version")
            == expected["version"],
            f"{label} installed distribution version mismatch",
        )


def _python_minor_version(version: str) -> str:
    parts = version.split(".")
    _require(len(parts) >= 2, f"Python version must include major.minor: {version!r}")
    return f"{parts[0]}.{parts[1]}"


def _required_mapping(value: object, *, field_name: str) -> dict[str, Any]:
    _require(isinstance(value, Mapping), f"{field_name} must be an object")
    return dict(value)


def _optional_mapping(value: object) -> dict[str, Any] | None:
    if value is None:
        return None
    _require(isinstance(value, Mapping), "optional evidence value must be an object")
    return dict(value)


def _optional_text(value: object) -> str | None:
    if value is None:
        return None
    return _required_text(value, field_name="optional text")


def _unique_text_list(value: object, *, field_name: str) -> list[str]:
    _require(isinstance(value, list), f"{field_name} must be an array")
    normalized: list[str] = []
    seen: set[str] = set()
    for position, item in enumerate(value):
        text = _required_text(item, field_name=f"{field_name}[{position}]")
        _require(text not in seen, f"{field_name} contains duplicate value: {text!r}")
        normalized.append(text)
        seen.add(text)
    _require(normalized != [], f"{field_name} must not be empty")
    return normalized


def _required_text(value: object, *, field_name: str) -> str:
    _require(
        isinstance(value, str) and bool(value.strip()),
        f"{field_name} must be non-empty text",
    )
    return value.strip()


def _require_sha256_digest(value: str, *, field_name: str) -> None:
    _require(value.startswith("sha256:"), f"{field_name} must use sha256: prefix")
    _require_sha256_hex(value.removeprefix("sha256:"), field_name=field_name)


def _require_sha256_hex(value: str, *, field_name: str) -> None:
    _require(
        len(value) == 64
        and all(character in "0123456789abcdef" for character in value),
        f"{field_name} must be 64 lowercase hexadecimal characters",
    )


def _prefixed_file_sha256(path: Path) -> str:
    return "sha256:" + _file_sha256(path)


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _utc_timestamp() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ReleaseAttestationError(message)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    path = build_release_attestation(
        policy_path=Path(args.policy),
        source_identity_path=Path(args.source_identity),
        build_manifest_path=Path(args.build_manifest),
        artifact_paths=tuple(Path(item) for item in _flatten(args.artifact)),
        source_check_report_paths=tuple(
            Path(item) for item in _flatten(args.source_check_report)
        ),
        artifact_verification_report_paths=tuple(
            Path(item) for item in _flatten(args.artifact_verification_report)
        ),
        output_path=Path(args.output),
        evidence_root=Path(args.evidence_root),
        workflow_run_id=args.workflow_run_id,
    )
    print(path.as_posix())
    return 0


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate verified release evidence into a success attestation."
    )
    parser.add_argument(
        "--policy",
        default="release/attestation-policy.json",
        help="Checked-in release evidence policy JSON.",
    )
    parser.add_argument(
        "--source-identity",
        required=True,
        help="Source identity record JSON.",
    )
    parser.add_argument("--build-manifest", required=True, help="Build manifest JSON.")
    parser.add_argument(
        "--artifact",
        action="append",
        nargs="+",
        required=True,
        help="Actual wheel or sdist file. Pass once for each required artifact.",
    )
    parser.add_argument(
        "--source-check-report",
        action="append",
        nargs="+",
        required=True,
        help="Successful source-check report JSON. Pass once per required report.",
    )
    parser.add_argument(
        "--artifact-verification-report",
        action="append",
        nargs="+",
        required=True,
        help=(
            "Successful installed-artifact verification report JSON. Pass once "
            "per required matrix cell."
        ),
    )
    parser.add_argument(
        "--output",
        default="build/release/release-attestation.json",
        help="Final success attestation path.",
    )
    parser.add_argument(
        "--evidence-root",
        default=".",
        help="Root used to store relative evidence filenames in the attestation.",
    )
    parser.add_argument(
        "--workflow-run-id",
        default=os.environ.get("GITHUB_RUN_ID"),
        help="Optional workflow or release-run identifier.",
    )
    return parser.parse_args(argv)


def _flatten(values: Sequence[Sequence[str]]) -> list[str]:
    return [item for group in values for item in group]


if __name__ == "__main__":
    raise SystemExit(main())
