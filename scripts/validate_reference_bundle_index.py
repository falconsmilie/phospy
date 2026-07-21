"""Validate staged reference-bundle files against staged manifests."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

REFERENCE_BUNDLES_ROOT = PurePosixPath("src/phospy/data/reference_bundles")
MANIFEST_FILENAME = "manifest.json"

_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_WINDOWS_DRIVE_PATTERN = re.compile(r"^[A-Za-z]:")


class ReferenceBundleIndexError(RuntimeError):
    """A staged reference-bundle file does not match its staged manifest."""


@dataclass(frozen=True, slots=True)
class StagedReferenceBundleFile:
    manifest_path: str
    reference_id: str | None
    relative_path: str
    file_path: str
    expected_sha256: str
    actual_sha256: str


@dataclass(frozen=True, slots=True)
class _DeclaredIndexFile:
    relative_path: str
    expected_sha256: str


@dataclass(frozen=True, slots=True)
class _ValidationIssue:
    manifest_path: str
    reference_id: str | None
    affected_file: str
    expected_digest: str
    actual_digest: str
    reason: str

    def __str__(self) -> str:
        reference_id = self.reference_id or "unknown"
        return (
            f"{self.reason}: "
            f"manifest path={self.manifest_path}; "
            f"reference ID={reference_id}; "
            f"affected file={self.affected_file}; "
            f"expected digest={self.expected_digest}; "
            f"actual digest={self.actual_digest}"
        )


def validate_reference_bundle_index(
    *,
    repo_root: str | Path | None = None,
    reference_root: str | PurePosixPath = REFERENCE_BUNDLES_ROOT,
) -> tuple[StagedReferenceBundleFile, ...]:
    """Validate staged manifest-declared reference files against staged blobs."""

    resolved_repo_root = _resolve_repo_root(repo_root)
    root = _normalize_reference_root(reference_root)
    issues, validated_files = _collect_index_issues(
        repo_root=resolved_repo_root,
        reference_root=root,
    )
    if issues:
        raise ReferenceBundleIndexError(
            "Reference bundle index validation failed:\n"
            + "\n".join(str(issue) for issue in issues)
        )
    return tuple(validated_files)


def _collect_index_issues(
    *,
    repo_root: Path,
    reference_root: PurePosixPath,
) -> tuple[list[_ValidationIssue], list[StagedReferenceBundleFile]]:
    index_files = _list_index_files(repo_root=repo_root, reference_root=reference_root)
    manifests = sorted(
        path
        for path in index_files
        if _is_reference_bundle_manifest(PurePosixPath(path), reference_root)
    )
    if not manifests:
        return [
            _ValidationIssue(
                manifest_path=(reference_root / MANIFEST_FILENAME).as_posix(),
                reference_id=None,
                affected_file=MANIFEST_FILENAME,
                expected_digest="at least one staged reference-bundle manifest",
                actual_digest="missing",
                reason="Git index contains no tracked reference-bundle manifests",
            )
        ], []

    issues: list[_ValidationIssue] = []
    validated_files: list[StagedReferenceBundleFile] = []
    required_paths: set[str] = set()
    for manifest_path in manifests:
        manifest_issues, manifest_files, manifest_required_paths = (
            _validate_manifest_from_index(
                repo_root=repo_root,
                manifest_path=manifest_path,
                index_files=index_files,
            )
        )
        issues.extend(manifest_issues)
        validated_files.extend(manifest_files)
        required_paths.update(manifest_required_paths)
    for extra_path in sorted(index_files - required_paths):
        issues.append(
            _ValidationIssue(
                manifest_path=str(PurePosixPath(extra_path).parent / MANIFEST_FILENAME),
                reference_id=None,
                affected_file=extra_path,
                expected_digest="manifest-declared reference-bundle file",
                actual_digest="extra",
                reason="Git index contains undeclared reference-bundle file",
            )
        )
    return issues, validated_files


def _validate_manifest_from_index(
    *,
    repo_root: Path,
    manifest_path: str,
    index_files: set[str],
) -> tuple[list[_ValidationIssue], list[StagedReferenceBundleFile], set[str]]:
    required_paths = {manifest_path}
    manifest_bytes, manifest_issue = _read_staged_blob(
        repo_root=repo_root,
        path=manifest_path,
        manifest_path=manifest_path,
        reference_id=None,
        affected_file=MANIFEST_FILENAME,
        expected_digest="staged manifest blob",
    )
    if manifest_issue is not None:
        return [manifest_issue], [], required_paths

    payload, payload_issue = _load_manifest_payload(
        raw_payload=manifest_bytes,
        manifest_path=manifest_path,
    )
    if payload_issue is not None:
        return [payload_issue], [], required_paths

    reference_id = _reference_id(payload)
    declared_files, declaration_issues = _declared_manifest_files(
        payload=payload,
        manifest_path=manifest_path,
        reference_id=reference_id,
    )
    issues = list(declaration_issues)
    validated_files: list[StagedReferenceBundleFile] = []
    bundle_path = PurePosixPath(manifest_path).parent

    for declared_file in declared_files:
        staged_path = bundle_path.joinpath(declared_file.relative_path).as_posix()
        required_paths.add(staged_path)
        if staged_path not in index_files:
            issues.append(
                _ValidationIssue(
                    manifest_path=manifest_path,
                    reference_id=reference_id,
                    affected_file=declared_file.relative_path,
                    expected_digest=declared_file.expected_sha256,
                    actual_digest="missing",
                    reason="manifest-listed file is missing from Git index",
                )
            )
            continue
        staged_bytes, blob_issue = _read_staged_blob(
            repo_root=repo_root,
            path=staged_path,
            manifest_path=manifest_path,
            reference_id=reference_id,
            affected_file=declared_file.relative_path,
            expected_digest=declared_file.expected_sha256,
        )
        if blob_issue is not None:
            issues.append(blob_issue)
            continue
        actual_sha256 = hashlib.sha256(staged_bytes).hexdigest()
        if actual_sha256 != declared_file.expected_sha256:
            issues.append(
                _ValidationIssue(
                    manifest_path=manifest_path,
                    reference_id=reference_id,
                    affected_file=declared_file.relative_path,
                    expected_digest=declared_file.expected_sha256,
                    actual_digest=actual_sha256,
                    reason="reference bundle staged blob digest mismatch",
                )
            )
            continue
        validated_files.append(
            StagedReferenceBundleFile(
                manifest_path=manifest_path,
                reference_id=reference_id,
                relative_path=declared_file.relative_path,
                file_path=staged_path,
                expected_sha256=declared_file.expected_sha256,
                actual_sha256=actual_sha256,
            )
        )

    return issues, validated_files, required_paths


def _load_manifest_payload(
    *,
    raw_payload: bytes,
    manifest_path: str,
) -> tuple[Mapping[str, Any], _ValidationIssue | None]:
    try:
        payload = json.loads(raw_payload.decode("utf-8"))
    except UnicodeDecodeError:
        return {}, _ValidationIssue(
            manifest_path=manifest_path,
            reference_id=None,
            affected_file=MANIFEST_FILENAME,
            expected_digest="valid UTF-8 JSON object",
            actual_digest="invalid",
            reason="staged reference-bundle manifest is not valid UTF-8",
        )
    except json.JSONDecodeError:
        return {}, _ValidationIssue(
            manifest_path=manifest_path,
            reference_id=None,
            affected_file=MANIFEST_FILENAME,
            expected_digest="valid JSON object",
            actual_digest="invalid",
            reason="staged reference-bundle manifest is not valid JSON",
        )
    if not isinstance(payload, Mapping):
        return {}, _ValidationIssue(
            manifest_path=manifest_path,
            reference_id=None,
            affected_file=MANIFEST_FILENAME,
            expected_digest="JSON object",
            actual_digest=type(payload).__name__,
            reason="staged reference-bundle manifest must decode to an object",
        )
    return payload, None


def _declared_manifest_files(
    *,
    payload: Mapping[str, Any],
    manifest_path: str,
    reference_id: str | None,
) -> tuple[list[_DeclaredIndexFile], list[_ValidationIssue]]:
    files = payload.get("files")
    if not isinstance(files, list) or not files:
        return [], [
            _ValidationIssue(
                manifest_path=manifest_path,
                reference_id=reference_id,
                affected_file="files",
                expected_digest="non-empty manifest files array",
                actual_digest="missing",
                reason="staged reference-bundle manifest has no declared files",
            )
        ]

    declared_files: list[_DeclaredIndexFile] = []
    issues: list[_ValidationIssue] = []
    seen_paths: set[str] = set()
    for index, item in enumerate(files):
        if not isinstance(item, Mapping):
            issues.append(
                _ValidationIssue(
                    manifest_path=manifest_path,
                    reference_id=reference_id,
                    affected_file=f"files[{index}]",
                    expected_digest="file manifest object",
                    actual_digest=type(item).__name__,
                    reason="manifest files entry is not an object",
                )
            )
            continue
        relative_path = _normalize_relative_manifest_path(
            item.get("relative_path"),
            manifest_path=manifest_path,
            reference_id=reference_id,
            affected_file=f"files[{index}].relative_path",
        )
        if isinstance(relative_path, _ValidationIssue):
            issues.append(relative_path)
            continue
        expected_sha256 = item.get("sha256")
        if not isinstance(expected_sha256, str) or not _SHA256_PATTERN.fullmatch(
            expected_sha256
        ):
            issues.append(
                _ValidationIssue(
                    manifest_path=manifest_path,
                    reference_id=reference_id,
                    affected_file=relative_path,
                    expected_digest="lowercase sha256 hex digest",
                    actual_digest=str(expected_sha256),
                    reason="manifest files entry has invalid sha256",
                )
            )
            continue
        if relative_path in seen_paths:
            issues.append(
                _ValidationIssue(
                    manifest_path=manifest_path,
                    reference_id=reference_id,
                    affected_file=relative_path,
                    expected_digest=expected_sha256,
                    actual_digest="duplicate",
                    reason="manifest files entry is duplicated",
                )
            )
            continue
        seen_paths.add(relative_path)
        declared_files.append(
            _DeclaredIndexFile(
                relative_path=relative_path,
                expected_sha256=expected_sha256,
            )
        )
    return declared_files, issues


def _normalize_relative_manifest_path(
    value: object,
    *,
    manifest_path: str,
    reference_id: str | None,
    affected_file: str,
) -> str | _ValidationIssue:
    if not isinstance(value, str) or not value.strip():
        return _ValidationIssue(
            manifest_path=manifest_path,
            reference_id=reference_id,
            affected_file=affected_file,
            expected_digest="relative POSIX path",
            actual_digest="missing",
            reason="manifest file path is missing or blank",
        )
    raw_path = value.strip()
    path = PurePosixPath(raw_path)
    if (
        "\\" in raw_path
        or _WINDOWS_DRIVE_PATTERN.match(raw_path) is not None
        or path.is_absolute()
        or ".." in path.parts
    ):
        return _ValidationIssue(
            manifest_path=manifest_path,
            reference_id=reference_id,
            affected_file=raw_path,
            expected_digest="relative POSIX path inside bundle",
            actual_digest="invalid",
            reason="manifest file path is not a safe bundle-relative POSIX path",
        )
    normalized = path.as_posix()
    if normalized in ("", "."):
        return _ValidationIssue(
            manifest_path=manifest_path,
            reference_id=reference_id,
            affected_file=raw_path,
            expected_digest="relative POSIX file path",
            actual_digest="invalid",
            reason="manifest file path must identify a file",
        )
    return normalized


def _read_staged_blob(
    *,
    repo_root: Path,
    path: str,
    manifest_path: str,
    reference_id: str | None,
    affected_file: str,
    expected_digest: str,
) -> tuple[bytes, _ValidationIssue | None]:
    result = _run_git(repo_root, ("show", f":{path}"))
    if result.returncode == 0:
        return result.stdout, None
    return b"", _ValidationIssue(
        manifest_path=manifest_path,
        reference_id=reference_id,
        affected_file=affected_file,
        expected_digest=expected_digest,
        actual_digest="missing",
        reason="Git index blob is missing",
    )


def _list_index_files(*, repo_root: Path, reference_root: PurePosixPath) -> set[str]:
    result = _run_git(repo_root, ("ls-files", "-z", "--", reference_root.as_posix()))
    if result.returncode != 0:
        stderr = result.stderr.decode("utf-8", errors="replace").strip()
        raise ReferenceBundleIndexError(
            "failed to list staged reference bundle files"
            + (f": {stderr}" if stderr else "")
        )
    return {
        item.decode("utf-8", errors="surrogateescape")
        for item in result.stdout.split(b"\0")
        if item
    }


def _run_git(
    repo_root: Path,
    args: Sequence[str],
) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ("git", *args),
        cwd=repo_root,
        capture_output=True,
        check=False,
    )


def _resolve_repo_root(repo_root: str | Path | None) -> Path:
    if repo_root is not None:
        return Path(repo_root)
    result = subprocess.run(
        ("git", "rev-parse", "--show-toplevel"),
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        stderr = result.stderr.decode("utf-8", errors="replace").strip()
        raise ReferenceBundleIndexError(
            "failed to resolve Git repository root" + (f": {stderr}" if stderr else "")
        )
    return Path(result.stdout.decode("utf-8", errors="replace").strip())


def _normalize_reference_root(value: str | PurePosixPath) -> PurePosixPath:
    path = PurePosixPath(str(value))
    if path.is_absolute() or ".." in path.parts:
        raise ReferenceBundleIndexError(
            f"reference root must be a repository-relative POSIX path: {value}"
        )
    return path


def _reference_id(payload: Mapping[str, Any]) -> str | None:
    value = payload.get("reference_id")
    if not isinstance(value, str) or not value.strip():
        return None
    return value.strip()


def _is_reference_bundle_manifest(
    path: PurePosixPath,
    reference_root: PurePosixPath,
) -> bool:
    root_parts = reference_root.parts
    if path.name != MANIFEST_FILENAME:
        return False
    if len(path.parts) <= len(root_parts):
        return False
    return path.parts[: len(root_parts)] == root_parts


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate staged reference-bundle files against staged manifests."
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=None,
        help="Git repository root. Defaults to git rev-parse --show-toplevel.",
    )
    parser.add_argument(
        "--reference-root",
        default=REFERENCE_BUNDLES_ROOT.as_posix(),
        help="Repository-relative reference bundle root.",
    )
    args = parser.parse_args(argv)

    try:
        validated_files = validate_reference_bundle_index(
            repo_root=args.repo_root,
            reference_root=args.reference_root,
        )
    except ReferenceBundleIndexError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    for item in validated_files:
        print(
            "validated staged reference bundle file: "
            f"file={item.file_path} "
            f"sha256={item.actual_sha256} "
            f"manifest={item.manifest_path}"
        )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
