from __future__ import annotations

import subprocess
from collections.abc import Iterable
from pathlib import Path


def changed_python_paths(
    repo: Path, *, base_commit: str, candidate_commit: str
) -> list[str]:
    completed = subprocess.run(
        (
            "git",
            "-C",
            str(repo),
            "diff",
            "--name-only",
            "--no-renames",
            "-z",
            base_commit,
            candidate_commit,
            "--",
            "src/phospy",
        ),
        check=True,
        capture_output=True,
    )
    return sorted(
        value.decode("utf-8")
        for value in completed.stdout.split(b"\0")
        if value and value.endswith(b".py")
    )


def source_maps(
    repo: Path,
    paths: Iterable[str],
    *,
    base_commit: str,
    candidate_commit: str,
) -> tuple[dict[str, bytes], dict[str, bytes]]:
    base: dict[str, bytes] = {}
    candidate: dict[str, bytes] = {}
    for path in paths:
        for revision, destination in (
            (base_commit, base),
            (candidate_commit, candidate),
        ):
            completed = subprocess.run(
                ("git", "-C", str(repo), "show", f"{revision}:{path}"),
                check=False,
                capture_output=True,
            )
            if completed.returncode == 0:
                destination[path] = completed.stdout
    return base, candidate
