"""Validate packaged reference-bundle files inside built distribution archives."""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import re
import subprocess
import sys
import tarfile
import zipfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

WHEEL_REFERENCE_BUNDLES_ROOT = PurePosixPath("phospy/data/reference_bundles")
SOURCE_REFERENCE_BUNDLES_ROOT = PurePosixPath("src/phospy/data/reference_bundles")
MANIFEST_FILENAME = "manifest.json"
DEFAULT_BUNDLE_ATTRIBUTION_PATH = "ATTRIBUTION.md"

_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


class ReferenceBundleDistributionError(RuntimeError):
    """A built distribution archive does not match bundled reference manifests."""


@dataclass(frozen=True, slots=True)
class _ArchiveReferenceFile:
    logical_path: str
    archive_entry: str
    data: bytes


@dataclass(frozen=True, slots=True)
class _ArchiveContents:
    files: dict[str, _ArchiveReferenceFile]
    issues: list[_ValidationIssue]


@dataclass(frozen=True, slots=True)
class _DeclaredArchiveFile:
    relative_path: str
    expected_sha256: str


@dataclass(frozen=True, slots=True)
class _ValidationIssue:
    archive_path: Path
    bundle_path: str
    reference_id: str | None
    affected_file: str
    expected_digest: str
    actual_digest: str
    reason: str

    def __str__(self) -> str:
        reference_id = self.reference_id or "unknown"
        return (
            f"{self.reason}: "
            f"archive path={self.archive_path}; "
            f"bundle path={self.bundle_path}; "
            f"reference ID={reference_id}; "
            f"affected file={self.affected_file}; "
            f"expected digest={self.expected_digest}; "
            f"actual digest={self.actual_digest}"
        )


def validate_reference_bundle_wheel(wheel_path: str | Path) -> None:
    """Validate all packaged reference-bundle manifests in one wheel."""

    validate_reference_bundle_archive(wheel_path)


def validate_reference_bundle_wheels(wheel_paths: Sequence[str | Path]) -> None:
    """Validate all supplied wheel paths and report every detected issue."""

    validate_reference_bundle_archives(wheel_paths)


def validate_reference_bundle_archive(
    archive_path: str | Path,
    *,
    compare_git_index: bool = False,
    repo_root: str | Path | None = None,
) -> None:
    """Validate all packaged reference-bundle manifests in one archive."""

    issues = _collect_archive_path_issues(
        Path(archive_path),
        compare_git_index=compare_git_index,
        repo_root=repo_root,
    )
    if issues:
        raise ReferenceBundleDistributionError(
            "Reference bundle distribution validation failed:\n"
            + "\n".join(str(issue) for issue in issues)
        )


def validate_reference_bundle_archives(
    archive_paths: Sequence[str | Path],
    *,
    compare_git_index: bool = False,
    repo_root: str | Path | None = None,
) -> None:
    """Validate all supplied distribution archives and report every issue."""

    issues: list[_ValidationIssue] = []
    for archive_path in archive_paths:
        issues.extend(
            _collect_archive_path_issues(
                Path(archive_path),
                compare_git_index=compare_git_index,
                repo_root=repo_root,
            )
        )
    if issues:
        raise ReferenceBundleDistributionError(
            "Reference bundle distribution validation failed:\n"
            + "\n".join(str(issue) for issue in issues)
        )


def _collect_archive_path_issues(
    archive_path: Path,
    *,
    compare_git_index: bool,
    repo_root: str | Path | None,
) -> list[_ValidationIssue]:
    if not archive_path.is_file():
        return [
            _ValidationIssue(
                archive_path=archive_path,
                bundle_path=SOURCE_REFERENCE_BUNDLES_ROOT.as_posix(),
                reference_id=None,
                affected_file=MANIFEST_FILENAME,
                expected_digest="present",
                actual_digest="missing",
                reason="distribution archive does not exist",
            )
        ]

    contents = _read_reference_files_from_archive(archive_path)
    issues = list(contents.issues)
    issues.extend(_collect_packaged_manifest_issues(archive_path, contents.files))
    if compare_git_index:
        issues.extend(
            _collect_git_index_comparison_issues(
                archive_path=archive_path,
                archive_files=contents.files,
                repo_root=repo_root,
            )
        )
    return issues


def _read_reference_files_from_archive(archive_path: Path) -> _ArchiveContents:
    try:
        if zipfile.is_zipfile(archive_path):
            with zipfile.ZipFile(archive_path) as archive:
                return _archive_contents_from_zip(archive_path, archive)
        if tarfile.is_tarfile(archive_path):
            with tarfile.open(archive_path) as archive:
                return _archive_contents_from_tar(archive_path, archive)
    except (tarfile.TarError, zipfile.BadZipFile, OSError) as exc:
        return _ArchiveContents(
            files={},
            issues=[
                _ValidationIssue(
                    archive_path=archive_path,
                    bundle_path=SOURCE_REFERENCE_BUNDLES_ROOT.as_posix(),
                    reference_id=None,
                    affected_file=MANIFEST_FILENAME,
                    expected_digest="valid wheel, sdist, or source zip archive",
                    actual_digest=type(exc).__name__,
                    reason="distribution archive could not be read",
                )
            ],
        )

    return _ArchiveContents(
        files={},
        issues=[
            _ValidationIssue(
                archive_path=archive_path,
                bundle_path=SOURCE_REFERENCE_BUNDLES_ROOT.as_posix(),
                reference_id=None,
                affected_file=MANIFEST_FILENAME,
                expected_digest="valid wheel, sdist, or source zip archive",
                actual_digest="invalid",
                reason="distribution archive is not a supported archive format",
            )
        ],
    )


def _archive_contents_from_zip(
    archive_path: Path,
    archive: zipfile.ZipFile,
) -> _ArchiveContents:
    files: dict[str, _ArchiveReferenceFile] = {}
    issues: list[_ValidationIssue] = []
    for info in archive.infolist():
        if info.is_dir():
            continue
        reference_file = _archive_reference_file(info.filename, archive.read(info))
        if reference_file is None:
            continue
        _add_reference_file(
            archive_path=archive_path,
            files=files,
            issues=issues,
            reference_file=reference_file,
        )
    return _ArchiveContents(files=files, issues=issues)


def _archive_contents_from_tar(
    archive_path: Path,
    archive: tarfile.TarFile,
) -> _ArchiveContents:
    files: dict[str, _ArchiveReferenceFile] = {}
    issues: list[_ValidationIssue] = []
    for member in archive.getmembers():
        if not member.isfile():
            continue
        extracted = archive.extractfile(member)
        if extracted is None:
            continue
        reference_file = _archive_reference_file(member.name, extracted.read())
        if reference_file is None:
            continue
        _add_reference_file(
            archive_path=archive_path,
            files=files,
            issues=issues,
            reference_file=reference_file,
        )
    return _ArchiveContents(files=files, issues=issues)


def _archive_reference_file(
    archive_entry: str,
    data: bytes,
) -> _ArchiveReferenceFile | None:
    normalized_entry = archive_entry.replace("\\", "/").lstrip("/")
    logical_path = _logical_source_reference_path(PurePosixPath(normalized_entry))
    if logical_path is None:
        return None
    return _ArchiveReferenceFile(
        logical_path=logical_path.as_posix(),
        archive_entry=normalized_entry,
        data=data,
    )


def _logical_source_reference_path(path: PurePosixPath) -> PurePosixPath | None:
    if _path_starts_with(path, WHEEL_REFERENCE_BUNDLES_ROOT):
        relative_parts = path.parts[len(WHEEL_REFERENCE_BUNDLES_ROOT.parts) :]
        return SOURCE_REFERENCE_BUNDLES_ROOT.joinpath(*relative_parts)

    source_parts = SOURCE_REFERENCE_BUNDLES_ROOT.parts
    for index in range(0, len(path.parts) - len(source_parts) + 1):
        if path.parts[index : index + len(source_parts)] == source_parts:
            relative_parts = path.parts[index + len(source_parts) :]
            return SOURCE_REFERENCE_BUNDLES_ROOT.joinpath(*relative_parts)
    return None


def _path_starts_with(path: PurePosixPath, root: PurePosixPath) -> bool:
    return (
        len(path.parts) >= len(root.parts)
        and path.parts[: len(root.parts)] == root.parts
    )


def _add_reference_file(
    *,
    archive_path: Path,
    files: dict[str, _ArchiveReferenceFile],
    issues: list[_ValidationIssue],
    reference_file: _ArchiveReferenceFile,
) -> None:
    existing = files.get(reference_file.logical_path)
    if existing is not None:
        issues.append(
            _ValidationIssue(
                archive_path=archive_path,
                bundle_path=str(PurePosixPath(reference_file.logical_path).parent),
                reference_id=None,
                affected_file=reference_file.logical_path,
                expected_digest=existing.archive_entry,
                actual_digest=reference_file.archive_entry,
                reason="duplicate packaged reference-bundle file",
            )
        )
        return
    files[reference_file.logical_path] = reference_file


def _collect_packaged_manifest_issues(
    archive_path: Path,
    archive_files: Mapping[str, _ArchiveReferenceFile],
) -> list[_ValidationIssue]:
    manifests = sorted(
        path
        for path in archive_files
        if _is_source_reference_bundle_manifest(PurePosixPath(path))
    )
    if not manifests:
        return [
            _ValidationIssue(
                archive_path=archive_path,
                bundle_path=SOURCE_REFERENCE_BUNDLES_ROOT.as_posix(),
                reference_id=None,
                affected_file=MANIFEST_FILENAME,
                expected_digest="at least one packaged manifest",
                actual_digest="missing",
                reason="distribution archive contains no reference-bundle manifests",
            )
        ]

    issues: list[_ValidationIssue] = []
    for manifest_path in manifests:
        issues.extend(
            _collect_manifest_issues(
                archive_path=archive_path,
                archive_files=archive_files,
                manifest_path=manifest_path,
            )
        )
    return issues


def _collect_manifest_issues(
    *,
    archive_path: Path,
    archive_files: Mapping[str, _ArchiveReferenceFile],
    manifest_path: str,
) -> list[_ValidationIssue]:
    manifest_file = archive_files[manifest_path]
    bundle_path = PurePosixPath(manifest_path).parent
    bundle_path_text = bundle_path.as_posix()
    payload, manifest_issue = _load_manifest_payload(
        raw_payload=manifest_file.data,
        archive_path=archive_path,
        bundle_path=bundle_path_text,
    )
    if manifest_issue is not None:
        return [manifest_issue]

    reference_id = _reference_id(payload)
    declared_files, issues = _declared_manifest_files(
        payload=payload,
        archive_path=archive_path,
        bundle_path=bundle_path_text,
        reference_id=reference_id,
    )
    attribution_path = _required_attribution_path(payload)
    normalized_attribution = _normalize_relative_manifest_path(
        attribution_path,
        archive_path=archive_path,
        bundle_path=bundle_path_text,
        reference_id=reference_id,
        affected_file="redistribution_evidence.attribution.bundle_attribution_path",
    )
    if isinstance(normalized_attribution, _ValidationIssue):
        issues.append(normalized_attribution)
        attribution_relative_path = DEFAULT_BUNDLE_ATTRIBUTION_PATH
    else:
        attribution_relative_path = normalized_attribution

    declared_by_path = {item.relative_path: item for item in declared_files}
    if attribution_relative_path not in declared_by_path:
        issues.append(
            _ValidationIssue(
                archive_path=archive_path,
                bundle_path=bundle_path_text,
                reference_id=reference_id,
                affected_file=attribution_relative_path,
                expected_digest="declared sha256",
                actual_digest="missing",
                reason=(
                    "required bundle-local attribution file is not declared in "
                    "manifest files"
                ),
            )
        )

    for declared_file in declared_files:
        logical_path = bundle_path.joinpath(declared_file.relative_path).as_posix()
        archive_file = archive_files.get(logical_path)
        if archive_file is None:
            reason = (
                "required bundle-local attribution file is missing from archive"
                if declared_file.relative_path == attribution_relative_path
                else "manifest-listed file is missing from archive"
            )
            issues.append(
                _ValidationIssue(
                    archive_path=archive_path,
                    bundle_path=bundle_path_text,
                    reference_id=reference_id,
                    affected_file=declared_file.relative_path,
                    expected_digest=declared_file.expected_sha256,
                    actual_digest="missing",
                    reason=reason,
                )
            )
            continue
        actual_sha256 = _sha256_bytes(archive_file.data)
        if actual_sha256 != declared_file.expected_sha256:
            issues.append(
                _ValidationIssue(
                    archive_path=archive_path,
                    bundle_path=bundle_path_text,
                    reference_id=reference_id,
                    affected_file=declared_file.relative_path,
                    expected_digest=declared_file.expected_sha256,
                    actual_digest=actual_sha256,
                    reason="reference bundle file digest mismatch",
                )
            )

    return issues


def _collect_git_index_comparison_issues(
    *,
    archive_path: Path,
    archive_files: Mapping[str, _ArchiveReferenceFile],
    repo_root: str | Path | None,
) -> list[_ValidationIssue]:
    root = _resolve_repo_root(repo_root)
    index_files, index_issue = _list_git_reference_files(
        archive_path=archive_path,
        repo_root=root,
    )
    if index_issue is not None:
        return [index_issue]

    required_paths, issues = _required_git_reference_paths(
        archive_path=archive_path,
        repo_root=root,
        index_files=index_files,
    )

    for logical_path in sorted(required_paths):
        archive_file = archive_files.get(logical_path)
        if archive_file is None:
            issues.append(
                _ValidationIssue(
                    archive_path=archive_path,
                    bundle_path=str(PurePosixPath(logical_path).parent),
                    reference_id=None,
                    affected_file=logical_path,
                    expected_digest="Git index blob",
                    actual_digest="missing",
                    reason=(
                        "committed reference-bundle file is missing from "
                        "distribution archive"
                    ),
                )
            )
            continue
        git_bytes, git_issue = _read_git_index_blob(
            archive_path=archive_path,
            repo_root=root,
            path=logical_path,
        )
        if git_issue is not None:
            issues.append(git_issue)
            continue
        expected_sha256 = _sha256_bytes(git_bytes)
        actual_sha256 = _sha256_bytes(archive_file.data)
        if actual_sha256 != expected_sha256:
            issues.append(
                _ValidationIssue(
                    archive_path=archive_path,
                    bundle_path=str(PurePosixPath(logical_path).parent),
                    reference_id=None,
                    affected_file=logical_path,
                    expected_digest=expected_sha256,
                    actual_digest=actual_sha256,
                    reason=(
                        "distribution reference-bundle file does not reproduce "
                        "committed Git index blob"
                    ),
                )
            )

    return issues


def _list_git_reference_files(
    *,
    archive_path: Path,
    repo_root: Path,
) -> tuple[set[str], _ValidationIssue | None]:
    result = _run_git(
        repo_root,
        ("ls-files", "-z", "--", SOURCE_REFERENCE_BUNDLES_ROOT.as_posix()),
    )
    if result.returncode == 0:
        return {
            item.decode("utf-8", errors="surrogateescape")
            for item in result.stdout.split(b"\0")
            if item
        }, None

    stderr = result.stderr.decode("utf-8", errors="replace").strip()
    return set(), _ValidationIssue(
        archive_path=archive_path,
        bundle_path=SOURCE_REFERENCE_BUNDLES_ROOT.as_posix(),
        reference_id=None,
        affected_file=MANIFEST_FILENAME,
        expected_digest="Git index file list",
        actual_digest=stderr or "unavailable",
        reason="failed to list committed reference-bundle files from Git index",
    )


def _required_git_reference_paths(
    *,
    archive_path: Path,
    repo_root: Path,
    index_files: set[str],
) -> tuple[set[str], list[_ValidationIssue]]:
    manifests = sorted(
        path
        for path in index_files
        if _is_source_reference_bundle_manifest(PurePosixPath(path))
    )
    if not manifests:
        return set(), [
            _ValidationIssue(
                archive_path=archive_path,
                bundle_path=SOURCE_REFERENCE_BUNDLES_ROOT.as_posix(),
                reference_id=None,
                affected_file=MANIFEST_FILENAME,
                expected_digest="at least one Git-index reference manifest",
                actual_digest="missing",
                reason="Git index contains no reference-bundle manifests",
            )
        ]

    required_paths: set[str] = set(manifests)
    issues: list[_ValidationIssue] = []
    for manifest_path in manifests:
        manifest_bytes, manifest_issue = _read_git_index_blob(
            archive_path=archive_path,
            repo_root=repo_root,
            path=manifest_path,
        )
        if manifest_issue is not None:
            issues.append(manifest_issue)
            continue
        payload, payload_issue = _load_manifest_payload(
            raw_payload=manifest_bytes,
            archive_path=archive_path,
            bundle_path=str(PurePosixPath(manifest_path).parent),
        )
        if payload_issue is not None:
            issues.append(payload_issue)
            continue
        reference_id = _reference_id(payload)
        declared_files, declaration_issues = _declared_manifest_files(
            payload=payload,
            archive_path=archive_path,
            bundle_path=str(PurePosixPath(manifest_path).parent),
            reference_id=reference_id,
        )
        issues.extend(declaration_issues)
        bundle_path = PurePosixPath(manifest_path).parent
        for declared_file in declared_files:
            required_paths.add(
                bundle_path.joinpath(declared_file.relative_path).as_posix()
            )
    return required_paths, issues


def _read_git_index_blob(
    *,
    archive_path: Path,
    repo_root: Path,
    path: str,
) -> tuple[bytes, _ValidationIssue | None]:
    result = _run_git(repo_root, ("show", f":{path}"))
    if result.returncode == 0:
        return result.stdout, None
    stderr = result.stderr.decode("utf-8", errors="replace").strip()
    return b"", _ValidationIssue(
        archive_path=archive_path,
        bundle_path=str(PurePosixPath(path).parent),
        reference_id=None,
        affected_file=path,
        expected_digest="Git index blob",
        actual_digest=stderr or "missing",
        reason="failed to read committed reference-bundle blob from Git index",
    )


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
    if result.returncode == 0:
        return Path(result.stdout.decode("utf-8", errors="replace").strip())
    stderr = result.stderr.decode("utf-8", errors="replace").strip()
    raise ReferenceBundleDistributionError(
        "failed to resolve Git repository root" + (f": {stderr}" if stderr else "")
    )


def _load_manifest_payload(
    *,
    raw_payload: bytes,
    archive_path: Path,
    bundle_path: str,
) -> tuple[Mapping[str, Any], _ValidationIssue | None]:
    try:
        payload = json.loads(raw_payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        return {}, _ValidationIssue(
            archive_path=archive_path,
            bundle_path=bundle_path,
            reference_id=None,
            affected_file=MANIFEST_FILENAME,
            expected_digest="valid JSON object",
            actual_digest="invalid",
            reason="reference-bundle manifest is not valid UTF-8 JSON",
        )
    if not isinstance(payload, Mapping):
        return {}, _ValidationIssue(
            archive_path=archive_path,
            bundle_path=bundle_path,
            reference_id=None,
            affected_file=MANIFEST_FILENAME,
            expected_digest="JSON object",
            actual_digest=type(payload).__name__,
            reason="reference-bundle manifest must decode to an object",
        )
    return payload, None


def _declared_manifest_files(
    *,
    payload: Mapping[str, Any],
    archive_path: Path,
    bundle_path: str,
    reference_id: str | None,
) -> tuple[list[_DeclaredArchiveFile], list[_ValidationIssue]]:
    files = payload.get("files")
    if not isinstance(files, list) or not files:
        return [], [
            _ValidationIssue(
                archive_path=archive_path,
                bundle_path=bundle_path,
                reference_id=reference_id,
                affected_file="files",
                expected_digest="non-empty manifest files array",
                actual_digest="missing",
                reason="reference-bundle manifest has no declared files",
            )
        ]

    declared_files: list[_DeclaredArchiveFile] = []
    issues: list[_ValidationIssue] = []
    seen_paths: set[str] = set()
    for index, item in enumerate(files):
        if not isinstance(item, Mapping):
            issues.append(
                _ValidationIssue(
                    archive_path=archive_path,
                    bundle_path=bundle_path,
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
            archive_path=archive_path,
            bundle_path=bundle_path,
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
                    archive_path=archive_path,
                    bundle_path=bundle_path,
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
                    archive_path=archive_path,
                    bundle_path=bundle_path,
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
            _DeclaredArchiveFile(
                relative_path=relative_path,
                expected_sha256=expected_sha256,
            )
        )
    return declared_files, issues


def _normalize_relative_manifest_path(
    value: object,
    *,
    archive_path: Path,
    bundle_path: str,
    reference_id: str | None,
    affected_file: str,
) -> str | _ValidationIssue:
    if not isinstance(value, str) or not value.strip():
        return _ValidationIssue(
            archive_path=archive_path,
            bundle_path=bundle_path,
            reference_id=reference_id,
            affected_file=affected_file,
            expected_digest="relative POSIX path",
            actual_digest="missing",
            reason="manifest file path is missing or blank",
        )
    raw_path = value.strip()
    path = PurePosixPath(raw_path)
    if "\\" in raw_path or path.is_absolute() or ".." in path.parts:
        return _ValidationIssue(
            archive_path=archive_path,
            bundle_path=bundle_path,
            reference_id=reference_id,
            affected_file=raw_path,
            expected_digest="relative POSIX path inside bundle",
            actual_digest="invalid",
            reason="manifest file path is not a safe bundle-relative POSIX path",
        )
    normalized = path.as_posix()
    if normalized in ("", "."):
        return _ValidationIssue(
            archive_path=archive_path,
            bundle_path=bundle_path,
            reference_id=reference_id,
            affected_file=raw_path,
            expected_digest="relative POSIX file path",
            actual_digest="invalid",
            reason="manifest file path must identify a file",
        )
    return normalized


def _required_attribution_path(payload: Mapping[str, Any]) -> str:
    evidence = payload.get("redistribution_evidence")
    if isinstance(evidence, Mapping):
        attribution = evidence.get("attribution")
        if isinstance(attribution, Mapping):
            bundle_path = attribution.get("bundle_attribution_path")
            if isinstance(bundle_path, str) and bundle_path.strip():
                return bundle_path
    return DEFAULT_BUNDLE_ATTRIBUTION_PATH


def _reference_id(payload: Mapping[str, Any]) -> str | None:
    value = payload.get("reference_id")
    if not isinstance(value, str) or not value.strip():
        return None
    return value.strip()


def _is_source_reference_bundle_manifest(path: PurePosixPath) -> bool:
    root_parts = SOURCE_REFERENCE_BUNDLES_ROOT.parts
    if path.name != MANIFEST_FILENAME:
        return False
    if len(path.parts) <= len(root_parts):
        return False
    return path.parts[: len(root_parts)] == root_parts


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _expand_archive_args(arguments: Sequence[str]) -> tuple[Path, ...]:
    paths: list[Path] = []
    for argument in arguments:
        matches = sorted(glob.glob(argument))
        if matches:
            file_matches = [Path(match) for match in matches if Path(match).is_file()]
            paths.extend(file_matches or (Path(match) for match in matches))
        else:
            paths.append(Path(argument))
    return tuple(paths)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Validate packaged reference-bundle manifests inside built wheel "
            "and source distribution archives."
        )
    )
    parser.add_argument(
        "archives",
        nargs="+",
        help="Distribution archive paths or glob patterns to validate.",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=None,
        help="Git repository root for Git-index comparison.",
    )
    parser.add_argument(
        "--no-git-index-compare",
        action="store_true",
        help="Skip comparison against committed Git-index reference files.",
    )
    args = parser.parse_args(argv)

    archive_paths = _expand_archive_args(args.archives)
    try:
        validate_reference_bundle_archives(
            archive_paths,
            compare_git_index=not args.no_git_index_compare,
            repo_root=args.repo_root,
        )
    except ReferenceBundleDistributionError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    for archive_path in archive_paths:
        print(f"validated reference bundles in {archive_path}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
