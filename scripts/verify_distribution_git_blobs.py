"""Compare packaged production files with the committed Git candidate bytes.

The verifier checks logical production members only. Wheel and sdist container
metadata, timestamps, compression, ordering, and other outer-archive
reproducibility properties are deliberately outside this gate.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tarfile
import zipfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

PRODUCTION_PREFIX = "src/phospy/"
WHEEL_PREFIX = "phospy/"


class DistributionGitBlobVerificationError(RuntimeError):
    """Raised when packaged production bytes differ from the Git candidate."""


@dataclass(frozen=True, slots=True)
class ArtifactProductionFiles:
    artifact: Path
    kind: str
    files: Mapping[str, bytes]


def compare_production_files(
    expected: Mapping[str, bytes],
    observed: Mapping[str, bytes],
) -> tuple[str, ...]:
    """Return exact path and byte differences for one distribution."""

    expected_paths = set(expected)
    observed_paths = set(observed)
    missing = expected_paths - observed_paths
    additional = observed_paths - expected_paths
    problems: list[str] = []

    unmatched_additional = set(additional)
    for missing_path in sorted(missing):
        folded = missing_path.casefold()
        renamed = sorted(
            path for path in unmatched_additional if path.casefold() == folded
        )
        if len(renamed) == 1:
            actual = renamed[0]
            unmatched_additional.remove(actual)
            problems.append(f"renamed production member: {missing_path} -> {actual}")
        else:
            problems.append(f"missing production member: {missing_path}")
    problems.extend(
        f"additional production member: {path}" for path in sorted(unmatched_additional)
    )
    for path in sorted(expected_paths & observed_paths):
        if expected[path] != observed[path]:
            problems.append(f"byte-different production member: {path}")
    return tuple(problems)


def git_production_blobs(repo_root: Path, revision: str) -> dict[str, bytes]:
    """Read every tracked ``src/phospy`` blob from one committed revision."""

    root = repo_root.resolve()
    listed = _git(
        root,
        "ls-tree",
        "-r",
        "-z",
        "--name-only",
        revision,
        "--",
        "src/phospy",
    )
    paths = [path for path in listed.decode("utf-8").split("\0") if path]
    expected: dict[str, bytes] = {}
    for path in paths:
        if not path.startswith(PRODUCTION_PREFIX):
            raise DistributionGitBlobVerificationError(
                f"unexpected Git production path: {path}"
            )
        _require_canonical_member(path)
        expected[path] = _git(root, "show", f"{revision}:{path}")
    if not expected:
        raise DistributionGitBlobVerificationError(
            f"revision {revision!r} has no tracked src/phospy files"
        )
    return expected


def read_artifact_production_files(path: Path) -> ArtifactProductionFiles:
    """Read and canonicalise every production member in a wheel or sdist."""

    artifact = path.resolve()
    if artifact.name.endswith(".whl"):
        return ArtifactProductionFiles(
            artifact=artifact,
            kind="wheel",
            files=_read_wheel(artifact),
        )
    if artifact.name.endswith(".tar.gz"):
        return ArtifactProductionFiles(
            artifact=artifact,
            kind="sdist",
            files=_read_sdist(artifact),
        )
    raise DistributionGitBlobVerificationError(
        f"unsupported distribution artifact: {artifact.name}"
    )


def verify_distribution_git_blobs(
    *,
    repo_root: Path,
    revision: str,
    artifacts: Sequence[Path],
) -> tuple[ArtifactProductionFiles, ...]:
    """Verify one wheel and one sdist against committed production blobs."""

    resolved_revision = _git(repo_root.resolve(), "rev-parse", f"{revision}^{{commit}}")
    commit = resolved_revision.decode("ascii").strip()
    expected = git_production_blobs(repo_root, commit)
    observed = tuple(read_artifact_production_files(path) for path in artifacts)
    kinds = [item.kind for item in observed]
    if kinds.count("wheel") != 1 or kinds.count("sdist") != 1 or len(kinds) != 2:
        raise DistributionGitBlobVerificationError(
            "expected exactly one wheel and one sdist artifact"
        )
    failures: list[str] = []
    for item in observed:
        failures.extend(
            f"{item.artifact.name}: {problem}"
            for problem in compare_production_files(expected, item.files)
        )
    if failures:
        raise DistributionGitBlobVerificationError("\n".join(failures))
    return observed


def _read_wheel(path: Path) -> dict[str, bytes]:
    files: dict[str, bytes] = {}
    with zipfile.ZipFile(path) as archive:
        for info in archive.infolist():
            name = info.filename
            if info.is_dir() or not name.startswith(WHEEL_PREFIX):
                continue
            canonical = f"src/{name}"
            _add_member(files, canonical, archive.read(info), artifact=path)
    return files


def _read_sdist(path: Path) -> dict[str, bytes]:
    files: dict[str, bytes] = {}
    with tarfile.open(path, mode="r:gz") as archive:
        roots = {
            PurePosixPath(member.name).parts[0]
            for member in archive.getmembers()
            if member.name and PurePosixPath(member.name).parts
        }
        if len(roots) != 1:
            raise DistributionGitBlobVerificationError(
                f"{path.name}: expected one sdist root directory, found {sorted(roots)}"
            )
        root = next(iter(roots))
        prefix = f"{root}/{PRODUCTION_PREFIX}"
        for member in archive.getmembers():
            if not member.isfile() or not member.name.startswith(prefix):
                continue
            canonical = member.name.removeprefix(f"{root}/")
            handle = archive.extractfile(member)
            if handle is None:
                raise DistributionGitBlobVerificationError(
                    f"{path.name}: cannot read production member {member.name}"
                )
            _add_member(files, canonical, handle.read(), artifact=path)
    return files


def _add_member(
    files: dict[str, bytes],
    path: str,
    content: bytes,
    *,
    artifact: Path,
) -> None:
    _require_canonical_member(path)
    if path in files:
        raise DistributionGitBlobVerificationError(
            f"{artifact.name}: duplicate production member {path}"
        )
    files[path] = content


def _require_canonical_member(path: str) -> None:
    pure = PurePosixPath(path)
    if (
        "\\" in path
        or path.startswith("/")
        or ".." in pure.parts
        or pure.as_posix() != path
        or not path.startswith(PRODUCTION_PREFIX)
    ):
        raise DistributionGitBlobVerificationError(
            f"non-canonical production member path: {path}"
        )


def _git(repo_root: Path, *arguments: str) -> bytes:
    try:
        return subprocess.run(
            ["git", *arguments],
            cwd=repo_root,
            check=True,
            capture_output=True,
        ).stdout
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.decode("utf-8", errors="replace").strip()
        raise DistributionGitBlobVerificationError(
            f"git {' '.join(arguments)} failed: {stderr}"
        ) from exc


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts", nargs="+", type=Path)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument("--revision", default="HEAD")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    try:
        artifacts = verify_distribution_git_blobs(
            repo_root=args.repo_root,
            revision=args.revision,
            artifacts=args.artifacts,
        )
    except DistributionGitBlobVerificationError as exc:
        print(f"distribution-to-Git verification failed: {exc}", file=sys.stderr)
        return 1
    counts = ", ".join(f"{item.kind}={len(item.files)}" for item in artifacts)
    print(
        "distribution-to-Git production bytes verified "
        f"for {args.revision} ({counts}); outer archive reproducibility is out of scope"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
