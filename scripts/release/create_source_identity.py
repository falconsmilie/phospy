from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 fallback
    import tomli as tomllib

SOURCE_IDENTITY_SCHEMA = "phospy.source-identity/v1"
SOURCE_DIGEST_ALGORITHM = "sha256"

_DEFAULT_TREE_EXCLUDED_DIRECTORIES = {
    ".git",
    ".hg",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    ".tox",
    ".venv",
    "__pycache__",
    "build",
    "dist",
    "site",
    "venv",
}


class SourceIdentityError(RuntimeError):
    """Raised when a source identity cannot be created."""


def write_source_identity(
    *,
    repository_root: Path,
    output_path: Path,
    source_archive: Path | None = None,
    source_archive_digest: str | None = None,
    package_name: str | None = None,
    package_version: str | None = None,
) -> Path:
    payload = build_source_identity(
        repository_root=repository_root,
        output_path=output_path,
        source_archive=source_archive,
        source_archive_digest=source_archive_digest,
        package_name=package_name,
        package_version=package_version,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return output_path


def build_source_identity(
    *,
    repository_root: Path,
    output_path: Path | None = None,
    source_archive: Path | None = None,
    source_archive_digest: str | None = None,
    package_name: str | None = None,
    package_version: str | None = None,
) -> dict[str, Any]:
    root = repository_root.resolve()
    _require(root.is_dir(), f"repository root does not exist: {root.as_posix()}")

    package = _package_metadata(
        root,
        package_name=package_name,
        package_version=package_version,
    )
    if _is_git_worktree(root):
        payload = _git_source_identity(root)
    elif source_archive is not None or source_archive_digest is not None:
        payload = _source_archive_identity(
            root=root,
            source_archive=source_archive,
            source_archive_digest=source_archive_digest,
        )
    else:
        excluded_paths: set[Path] = set()
        if output_path is not None:
            excluded_paths.add(output_path.resolve())
        payload = _normalized_tree_source_identity(
            root,
            excluded_paths=excluded_paths,
        )
    if package is not None:
        payload["package"] = package
    return payload


def source_identity_digest(path: Path) -> str:
    return f"{SOURCE_DIGEST_ALGORITHM}:{_file_sha256(path)}"


def _source_archive_identity(
    *,
    root: Path,
    source_archive: Path | None,
    source_archive_digest: str | None,
) -> dict[str, Any]:
    digest = source_archive_digest
    archive_payload: dict[str, Any] = {}
    if source_archive is not None:
        archive_path = source_archive.resolve()
        _require(
            archive_path.is_file(),
            f"source archive does not exist: {archive_path.as_posix()}",
        )
        digest = source_identity_digest(archive_path)
        archive_payload["filename"] = archive_path.name
    _require(digest is not None, "source archive digest is required")
    _require_sha256_digest(digest, field_name="source_archive_digest")
    return {
        "schema": SOURCE_IDENTITY_SCHEMA,
        "method": "source-archive",
        "source_archive": {
            **archive_payload,
            "digest": digest,
        },
        "repository": {
            "root_name": root.name,
        },
    }


def _git_source_identity(root: Path) -> dict[str, Any]:
    revision = _git_stdout(root, "rev-parse", "HEAD")
    dirty_entries = _git_status_entries(root)
    tracked_paths = _git_paths(root, "ls-files", "-z")
    untracked_paths = _git_paths(
        root, "ls-files", "-z", "--others", "--exclude-standard"
    )
    digest_paths = sorted(
        set(tracked_paths).union(untracked_paths),
        key=lambda path: path.as_posix(),
    )
    return {
        "schema": SOURCE_IDENTITY_SCHEMA,
        "method": "git",
        "git": {
            "revision": revision,
            "dirty": bool(dirty_entries),
            "dirty_paths": dirty_entries,
            "tracked_file_count": len(tracked_paths),
            "untracked_file_count": len(untracked_paths),
        },
        "source_tree": {
            "digest": _source_tree_digest(root, digest_paths),
            "file_count": len(digest_paths),
            "algorithm": SOURCE_DIGEST_ALGORITHM,
            "file_selection": "git tracked plus untracked non-ignored files",
        },
    }


def _normalized_tree_source_identity(
    root: Path,
    *,
    excluded_paths: set[Path],
) -> dict[str, Any]:
    paths = _normalized_tree_paths(root, excluded_paths=excluded_paths)
    return {
        "schema": SOURCE_IDENTITY_SCHEMA,
        "method": "source-tree",
        "source_tree": {
            "digest": _source_tree_digest(root, paths),
            "file_count": len(paths),
            "algorithm": SOURCE_DIGEST_ALGORITHM,
            "file_selection": "normalized filesystem tree excluding generated directories",
        },
    }


def _source_tree_digest(root: Path, relative_paths: Iterable[Path]) -> str:
    digest = hashlib.sha256()
    for relative_path in sorted(relative_paths, key=lambda path: path.as_posix()):
        normalized = relative_path.as_posix()
        path = root / relative_path
        digest.update(normalized.encode("utf-8"))
        digest.update(b"\0")
        if path.is_file():
            file_digest = _file_sha256(path)
            size = path.stat().st_size
            digest.update(str(size).encode("ascii"))
            digest.update(b"\0")
            digest.update(file_digest.encode("ascii"))
        else:
            digest.update(b"missing")
        digest.update(b"\n")
    return f"{SOURCE_DIGEST_ALGORITHM}:{digest.hexdigest()}"


def _normalized_tree_paths(root: Path, *, excluded_paths: set[Path]) -> list[Path]:
    paths: list[Path] = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        resolved = path.resolve()
        if resolved in excluded_paths:
            continue
        relative = path.relative_to(root)
        if _is_generated_or_vcs_path(relative):
            continue
        paths.append(relative)
    return sorted(paths, key=lambda path: path.as_posix())


def _is_generated_or_vcs_path(relative_path: Path) -> bool:
    return any(
        part in _DEFAULT_TREE_EXCLUDED_DIRECTORIES for part in relative_path.parts
    )


def _is_git_worktree(root: Path) -> bool:
    try:
        result = subprocess.run(
            ("git", "rev-parse", "--is-inside-work-tree"),
            cwd=root,
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError:
        return False
    return result.returncode == 0 and result.stdout.strip() == "true"


def _git_stdout(root: Path, *args: str) -> str:
    result = subprocess.run(
        ("git", *args),
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def _git_paths(root: Path, *args: str) -> list[Path]:
    result = subprocess.run(
        ("git", *args),
        cwd=root,
        check=True,
        capture_output=True,
    )
    return sorted(
        (Path(raw.decode("utf-8")) for raw in result.stdout.split(b"\0") if raw),
        key=lambda path: path.as_posix(),
    )


def _git_status_entries(root: Path) -> list[str]:
    result = subprocess.run(
        ("git", "status", "--porcelain=v1", "--untracked-files=normal"),
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    return [line.rstrip() for line in result.stdout.splitlines() if line.strip()]


def _package_metadata(
    root: Path,
    *,
    package_name: str | None,
    package_version: str | None,
) -> dict[str, str] | None:
    pyproject_package = _pyproject_package(root)
    name = package_name or (
        None if pyproject_package is None else pyproject_package.get("name")
    )
    version = package_version or (
        None if pyproject_package is None else pyproject_package.get("version")
    )
    if name is None and version is None:
        return None
    _require(
        isinstance(name, str) and bool(name.strip()),
        "package name is required when package metadata is recorded",
    )
    _require(
        isinstance(version, str) and bool(version.strip()),
        "package version is required when package metadata is recorded",
    )
    return {"name": name.strip(), "version": version.strip()}


def _pyproject_package(root: Path) -> dict[str, str] | None:
    pyproject = root / "pyproject.toml"
    if not pyproject.is_file():
        return None
    with pyproject.open("rb") as handle:
        payload = tomllib.load(handle)
    project = payload.get("project")
    if not isinstance(project, Mapping):
        return None
    name = project.get("name")
    version = project.get("version")
    if not isinstance(name, str) or not isinstance(version, str):
        return None
    return {"name": name, "version": version}


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require_sha256_digest(value: str, *, field_name: str) -> None:
    _require(value.startswith("sha256:"), f"{field_name} must use sha256: prefix")
    _require_sha256_hex(value.removeprefix("sha256:"), field_name=field_name)


def _require_sha256_hex(value: str, *, field_name: str) -> None:
    _require(
        len(value) == 64
        and all(character in "0123456789abcdef" for character in value),
        f"{field_name} must be 64 lowercase hexadecimal characters",
    )


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise SourceIdentityError(message)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    path = write_source_identity(
        repository_root=Path(args.repository_root),
        output_path=Path(args.output),
        source_archive=None
        if args.source_archive is None
        else Path(args.source_archive),
        source_archive_digest=args.source_archive_digest,
        package_name=args.package_name,
        package_version=args.package_version,
    )
    print(path.as_posix())
    return 0


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a deterministic PhosPy release source identity record."
    )
    parser.add_argument(
        "--repository-root",
        default=".",
        help="Source checkout or unpacked source tree root.",
    )
    parser.add_argument(
        "--output",
        default="build/reports/source-identity.json",
        help="Path where the source identity JSON record is written.",
    )
    parser.add_argument(
        "--source-archive",
        default=None,
        help="Optional source archive to bind instead of inspecting Git metadata.",
    )
    parser.add_argument(
        "--source-archive-digest",
        default=None,
        help="Optional precomputed sha256:<hex> source archive digest.",
    )
    parser.add_argument(
        "--package-name",
        default=None,
        help="Package name override when pyproject.toml is unavailable.",
    )
    parser.add_argument(
        "--package-version",
        default=None,
        help="Package version override when pyproject.toml is unavailable.",
    )
    return parser.parse_args(argv)


if __name__ == "__main__":
    raise SystemExit(main())
