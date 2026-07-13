"""Validate packaged reference-bundle files inside built wheel archives."""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import re
import sys
import zipfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

REFERENCE_BUNDLES_ROOT = PurePosixPath("phospy/data/reference_bundles")
MANIFEST_FILENAME = "manifest.json"
DEFAULT_BUNDLE_ATTRIBUTION_PATH = "ATTRIBUTION.md"

_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


class ReferenceBundleDistributionError(RuntimeError):
    """A built wheel does not match bundled reference manifests."""


@dataclass(frozen=True, slots=True)
class _DeclaredWheelFile:
    relative_path: str
    expected_sha256: str


@dataclass(frozen=True, slots=True)
class _ValidationIssue:
    wheel_path: Path
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
            f"wheel path={self.wheel_path}; "
            f"bundle path={self.bundle_path}; "
            f"reference ID={reference_id}; "
            f"affected file={self.affected_file}; "
            f"expected digest={self.expected_digest}; "
            f"actual digest={self.actual_digest}"
        )


def validate_reference_bundle_wheel(wheel_path: str | Path) -> None:
    """Validate all packaged reference-bundle manifests in one wheel."""

    path = Path(wheel_path)
    issues = _collect_wheel_issues(path)
    if issues:
        raise ReferenceBundleDistributionError(
            "Reference bundle distribution validation failed:\n"
            + "\n".join(str(issue) for issue in issues)
        )


def validate_reference_bundle_wheels(wheel_paths: Sequence[str | Path]) -> None:
    """Validate all supplied wheel paths and report every detected issue."""

    issues: list[_ValidationIssue] = []
    for wheel_path in wheel_paths:
        issues.extend(_collect_wheel_issues(Path(wheel_path)))
    if issues:
        raise ReferenceBundleDistributionError(
            "Reference bundle distribution validation failed:\n"
            + "\n".join(str(issue) for issue in issues)
        )


def _collect_wheel_issues(wheel_path: Path) -> list[_ValidationIssue]:
    if not wheel_path.is_file():
        return [
            _ValidationIssue(
                wheel_path=wheel_path,
                bundle_path=REFERENCE_BUNDLES_ROOT.as_posix(),
                reference_id=None,
                affected_file=MANIFEST_FILENAME,
                expected_digest="present",
                actual_digest="missing",
                reason="wheel file does not exist",
            )
        ]
    try:
        with zipfile.ZipFile(wheel_path) as archive:
            return _collect_archive_issues(wheel_path, archive)
    except zipfile.BadZipFile:
        return [
            _ValidationIssue(
                wheel_path=wheel_path,
                bundle_path=REFERENCE_BUNDLES_ROOT.as_posix(),
                reference_id=None,
                affected_file=MANIFEST_FILENAME,
                expected_digest="valid wheel zip archive",
                actual_digest="invalid",
                reason="wheel file is not a valid ZIP archive",
            )
        ]


def _collect_archive_issues(
    wheel_path: Path,
    archive: zipfile.ZipFile,
) -> list[_ValidationIssue]:
    issues: list[_ValidationIssue] = []
    archive_files = {info.filename for info in archive.infolist() if not info.is_dir()}
    manifests = sorted(
        name
        for name in archive_files
        if _is_reference_bundle_manifest(PurePosixPath(name))
    )
    if not manifests:
        return [
            _ValidationIssue(
                wheel_path=wheel_path,
                bundle_path=REFERENCE_BUNDLES_ROOT.as_posix(),
                reference_id=None,
                affected_file=MANIFEST_FILENAME,
                expected_digest="at least one packaged manifest",
                actual_digest="missing",
                reason="wheel contains no packaged reference-bundle manifests",
            )
        ]

    for manifest_entry in manifests:
        issues.extend(
            _collect_manifest_issues(
                wheel_path=wheel_path,
                archive=archive,
                archive_files=archive_files,
                manifest_entry=manifest_entry,
            )
        )
    return issues


def _collect_manifest_issues(
    *,
    wheel_path: Path,
    archive: zipfile.ZipFile,
    archive_files: set[str],
    manifest_entry: str,
) -> list[_ValidationIssue]:
    bundle_path = PurePosixPath(manifest_entry).parent
    bundle_path_text = bundle_path.as_posix()
    payload, manifest_issue = _load_manifest_payload(
        wheel_path=wheel_path,
        archive=archive,
        manifest_entry=manifest_entry,
        bundle_path=bundle_path_text,
    )
    if manifest_issue is not None:
        return [manifest_issue]

    reference_id = _reference_id(payload)
    declared_files, issues = _declared_manifest_files(
        payload=payload,
        wheel_path=wheel_path,
        bundle_path=bundle_path_text,
        reference_id=reference_id,
    )
    attribution_path = _required_attribution_path(payload)
    normalized_attribution = _normalize_relative_manifest_path(
        attribution_path,
        wheel_path=wheel_path,
        bundle_path=bundle_path_text,
        reference_id=reference_id,
        affected_file="redistribution_evidence.attribution.bundle_attribution_path",
    )
    if isinstance(normalized_attribution, _ValidationIssue):
        issues.append(normalized_attribution)
        attribution_relative_path: str | None = DEFAULT_BUNDLE_ATTRIBUTION_PATH
    else:
        attribution_relative_path = normalized_attribution

    declared_by_path = {item.relative_path: item for item in declared_files}
    if attribution_relative_path not in declared_by_path:
        issues.append(
            _ValidationIssue(
                wheel_path=wheel_path,
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
        archive_entry = bundle_path.joinpath(declared_file.relative_path).as_posix()
        if archive_entry not in archive_files:
            reason = (
                "required bundle-local attribution file is missing from wheel"
                if declared_file.relative_path == attribution_relative_path
                else "manifest-listed file is missing from wheel"
            )
            issues.append(
                _ValidationIssue(
                    wheel_path=wheel_path,
                    bundle_path=bundle_path_text,
                    reference_id=reference_id,
                    affected_file=declared_file.relative_path,
                    expected_digest=declared_file.expected_sha256,
                    actual_digest="missing",
                    reason=reason,
                )
            )
            continue
        actual_sha256 = _sha256_bytes(archive.read(archive_entry))
        if actual_sha256 != declared_file.expected_sha256:
            issues.append(
                _ValidationIssue(
                    wheel_path=wheel_path,
                    bundle_path=bundle_path_text,
                    reference_id=reference_id,
                    affected_file=declared_file.relative_path,
                    expected_digest=declared_file.expected_sha256,
                    actual_digest=actual_sha256,
                    reason="reference bundle file digest mismatch",
                )
            )

    if attribution_relative_path is not None:
        attribution_entry = bundle_path.joinpath(attribution_relative_path).as_posix()
        if (
            attribution_relative_path not in declared_by_path
            and attribution_entry not in archive_files
        ):
            issues.append(
                _ValidationIssue(
                    wheel_path=wheel_path,
                    bundle_path=bundle_path_text,
                    reference_id=reference_id,
                    affected_file=attribution_relative_path,
                    expected_digest="present",
                    actual_digest="missing",
                    reason="required bundle-local attribution file is missing from wheel",
                )
            )

    return issues


def _load_manifest_payload(
    *,
    wheel_path: Path,
    archive: zipfile.ZipFile,
    manifest_entry: str,
    bundle_path: str,
) -> tuple[Mapping[str, Any], _ValidationIssue | None]:
    try:
        raw_payload = archive.read(manifest_entry)
    except KeyError:
        return {}, _ValidationIssue(
            wheel_path=wheel_path,
            bundle_path=bundle_path,
            reference_id=None,
            affected_file=MANIFEST_FILENAME,
            expected_digest="present",
            actual_digest="missing",
            reason="packaged reference-bundle manifest is missing from wheel",
        )
    try:
        payload = json.loads(raw_payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        return {}, _ValidationIssue(
            wheel_path=wheel_path,
            bundle_path=bundle_path,
            reference_id=None,
            affected_file=MANIFEST_FILENAME,
            expected_digest="valid JSON object",
            actual_digest="invalid",
            reason="packaged reference-bundle manifest is not valid UTF-8 JSON",
        )
    if not isinstance(payload, Mapping):
        return {}, _ValidationIssue(
            wheel_path=wheel_path,
            bundle_path=bundle_path,
            reference_id=None,
            affected_file=MANIFEST_FILENAME,
            expected_digest="JSON object",
            actual_digest=type(payload).__name__,
            reason="packaged reference-bundle manifest must decode to an object",
        )
    return payload, None


def _declared_manifest_files(
    *,
    payload: Mapping[str, Any],
    wheel_path: Path,
    bundle_path: str,
    reference_id: str | None,
) -> tuple[list[_DeclaredWheelFile], list[_ValidationIssue]]:
    files = payload.get("files")
    if not isinstance(files, list) or not files:
        return [], [
            _ValidationIssue(
                wheel_path=wheel_path,
                bundle_path=bundle_path,
                reference_id=reference_id,
                affected_file="files",
                expected_digest="non-empty manifest files array",
                actual_digest="missing",
                reason="packaged reference-bundle manifest has no declared files",
            )
        ]

    declared_files: list[_DeclaredWheelFile] = []
    issues: list[_ValidationIssue] = []
    seen_paths: set[str] = set()
    for index, item in enumerate(files):
        if not isinstance(item, Mapping):
            issues.append(
                _ValidationIssue(
                    wheel_path=wheel_path,
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
            wheel_path=wheel_path,
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
                    wheel_path=wheel_path,
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
                    wheel_path=wheel_path,
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
            _DeclaredWheelFile(
                relative_path=relative_path,
                expected_sha256=expected_sha256,
            )
        )
    return declared_files, issues


def _normalize_relative_manifest_path(
    value: object,
    *,
    wheel_path: Path,
    bundle_path: str,
    reference_id: str | None,
    affected_file: str,
) -> str | _ValidationIssue:
    if not isinstance(value, str) or not value.strip():
        return _ValidationIssue(
            wheel_path=wheel_path,
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
            wheel_path=wheel_path,
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
            wheel_path=wheel_path,
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


def _is_reference_bundle_manifest(path: PurePosixPath) -> bool:
    root_parts = REFERENCE_BUNDLES_ROOT.parts
    if path.name != MANIFEST_FILENAME:
        return False
    if len(path.parts) <= len(root_parts):
        return False
    return path.parts[: len(root_parts)] == root_parts


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _expand_wheel_args(arguments: Sequence[str]) -> tuple[Path, ...]:
    paths: list[Path] = []
    for argument in arguments:
        matches = sorted(glob.glob(argument))
        if matches:
            paths.extend(Path(match) for match in matches)
        else:
            paths.append(Path(argument))
    return tuple(paths)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate packaged reference-bundle manifests inside wheels."
    )
    parser.add_argument(
        "wheels",
        nargs="+",
        help="Wheel paths or glob patterns to validate.",
    )
    args = parser.parse_args(argv)

    wheel_paths = _expand_wheel_args(args.wheels)
    try:
        validate_reference_bundle_wheels(wheel_paths)
    except ReferenceBundleDistributionError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    for wheel_path in wheel_paths:
        print(f"validated reference bundles in {wheel_path}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
