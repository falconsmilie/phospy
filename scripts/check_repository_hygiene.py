"""Enforce permanent repository hygiene rules."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

FORBIDDEN_TRACKED_PREFIXES = ("artifacts/", "docs/refactor/")
LARGE_FILE_REVIEW_BYTES = 50 * 1024 * 1024
GITHUB_FILE_LIMIT_BYTES = 100 * 1024 * 1024


def _tracked_existing_paths(repo_root: Path) -> tuple[str, ...]:
    completed = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=repo_root,
        check=True,
        capture_output=True,
    )
    tracked = (item.decode("utf-8") for item in completed.stdout.split(b"\0") if item)
    return tuple(path for path in tracked if (repo_root / path).is_file())


def hygiene_problems(repo_root: Path) -> list[str]:
    tracked = _tracked_existing_paths(repo_root)
    problems = [
        f"forbidden tracked path: {path}"
        for path in tracked
        if path.startswith(FORBIDDEN_TRACKED_PREFIXES)
    ]
    if "/artifacts/" not in (repo_root / ".gitignore").read_text(encoding="utf-8"):
        problems.append(".gitignore must contain /artifacts/")

    for path in tracked:
        size = (repo_root / path).stat().st_size
        if size > GITHUB_FILE_LIMIT_BYTES:
            problems.append(f"tracked file exceeds GitHub limit: {path} ({size} bytes)")
        elif size > LARGE_FILE_REVIEW_BYTES:
            problems.append(f"tracked file exceeds 50 MiB: {path} ({size} bytes)")
    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    args = parser.parse_args(argv)
    problems = hygiene_problems(args.repo_root.resolve())
    if problems:
        print("\n".join(problems))
        return 1
    print("repository hygiene: clean")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
