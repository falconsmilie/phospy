from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from collections.abc import Sequence
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 fallback
    import tomli as tomllib

BUILD_MANIFEST_SCHEMA = "phospy.build-manifest/v1"


def write_build_manifest(
    *,
    dist_dir: Path,
    output_path: Path,
    repository_root: Path,
    source_identity_path: Path | None = None,
    package_name: str | None = None,
    package_version: str | None = None,
) -> Path:
    root = repository_root.resolve()
    package = _package_metadata(
        root,
        package_name=package_name,
        package_version=package_version,
    )
    artifacts = _distribution_artifacts(dist_dir.resolve())
    payload = {
        "schema": BUILD_MANIFEST_SCHEMA,
        "source_identity_digest": _source_identity_digest(root, source_identity_path),
        "package_name": package["name"],
        "package_version": package["version"],
        "artifacts": artifacts,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return output_path


def _distribution_artifacts(dist_dir: Path) -> list[dict[str, str]]:
    wheels = sorted(dist_dir.glob("*.whl"))
    sdists = sorted(dist_dir.glob("*.tar.gz"))
    _require(len(wheels) == 1, f"expected exactly one wheel in {dist_dir.as_posix()}")
    _require(len(sdists) == 1, f"expected exactly one sdist in {dist_dir.as_posix()}")
    return [
        _artifact_payload("wheel", wheels[0]),
        _artifact_payload("sdist", sdists[0]),
    ]


def _artifact_payload(kind: str, path: Path) -> dict[str, str]:
    return {
        "kind": kind,
        "filename": path.name,
        "sha256": _file_sha256(path),
    }


def _source_identity_digest(
    repository_root: Path,
    source_identity_path: Path | None = None,
) -> str:
    if source_identity_path is not None:
        path = source_identity_path.resolve()
        _require(
            path.is_file(), f"source identity record is missing: {path.as_posix()}"
        )
        return "sha256:" + _file_sha256(path)

    digest = hashlib.sha256()
    for relative_path in _git_tracked_files(repository_root):
        path = repository_root / relative_path
        digest.update(relative_path.as_posix().encode("utf-8"))
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
    return "sha256:" + digest.hexdigest()


def _git_tracked_files(repository_root: Path) -> list[Path]:
    result = subprocess.run(
        ("git", "ls-files", "-z"),
        cwd=repository_root,
        check=True,
        capture_output=True,
    )
    paths = [Path(raw.decode("utf-8")) for raw in result.stdout.split(b"\0") if raw]
    _require(paths, "git tracked source file list is empty")
    return sorted(paths, key=lambda path: path.as_posix())


def _package_metadata(
    repository_root: Path,
    *,
    package_name: str | None = None,
    package_version: str | None = None,
) -> dict[str, str]:
    project: dict[str, object] = {}
    if package_name is None or package_version is None:
        with (repository_root / "pyproject.toml").open("rb") as handle:
            payload = tomllib.load(handle)
        project = payload["project"]
    name = package_name or project["name"]
    version = package_version or project["version"]
    _require(isinstance(name, str) and bool(name.strip()), "package name missing")
    _require(
        isinstance(version, str) and bool(version.strip()), "package version missing"
    )
    return {"name": name.strip(), "version": version.strip()}


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    path = write_build_manifest(
        dist_dir=Path(args.dist_dir),
        output_path=Path(args.output),
        repository_root=Path(args.repository_root),
        source_identity_path=(
            None if args.source_identity is None else Path(args.source_identity)
        ),
        package_name=args.package_name,
        package_version=args.package_version,
    )
    print(path.as_posix())
    return 0


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write a PhosPy build manifest for distribution artifacts."
    )
    parser.add_argument(
        "--dist-dir",
        default="dist",
        help="Directory containing exactly one wheel and one sdist.",
    )
    parser.add_argument(
        "--output",
        default="build/reports/build-manifest.json",
        help="Path where the build manifest JSON is written.",
    )
    parser.add_argument(
        "--repository-root",
        default=".",
        help="Repository root used to compute the source identity digest.",
    )
    parser.add_argument(
        "--package-version",
        default=None,
        help="Package version override; defaults to project.version from pyproject.toml.",
    )
    parser.add_argument(
        "--package-name",
        default=None,
        help="Package name override; defaults to project.name from pyproject.toml.",
    )
    parser.add_argument(
        "--source-identity",
        default=None,
        help=(
            "Source identity record to bind into the build manifest. When omitted, "
            "the legacy Git-tracked tree digest is used."
        ),
    )
    return parser.parse_args(argv)


if __name__ == "__main__":
    raise SystemExit(main())
