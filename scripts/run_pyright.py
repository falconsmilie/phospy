#!/usr/bin/env python3
"""Run pyright against the active interpreter."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def _candidate_interpreters(repo_root: Path) -> list[str]:
    candidates: list[str] = []
    override = os.environ.get("PHOSPY_PYRIGHT_PYTHON")
    if override:
        candidates.append(override)

    venv_dirs = [".venv-ci", ".venv", "venv"]
    windows_rel = Path("Scripts") / "python.exe"
    unix_rel = Path("bin") / "python"
    for name in venv_dirs:
        candidates.append(str(repo_root / name / windows_rel))
        candidates.append(str(repo_root / name / unix_rel))

    candidates.append(sys.executable)
    seen: set[str] = set()
    unique: list[str] = []
    for candidate in candidates:
        normalized = str(Path(candidate))
        if normalized in seen:
            continue
        seen.add(normalized)
        unique.append(normalized)
    return unique


def _supports_pyright(python_bin: str) -> bool:
    check = subprocess.run(
        [
            python_bin,
            "-c",
            "import pyright, numpy, pandas, scipy, sklearn",
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return check.returncode == 0


def _resolve_python_bin(repo_root: Path) -> tuple[str | None, list[str]]:
    checked: list[str] = []
    for candidate in _candidate_interpreters(repo_root):
        checked.append(candidate)
        if not Path(candidate).exists():
            continue
        if _supports_pyright(candidate):
            return candidate, checked
    return None, checked


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    python_bin, checked = _resolve_python_bin(repo_root)
    if python_bin is None:
        sys.stderr.write(
            "Unable to find a Python interpreter with required packages for Pyright.\n"
            "Checked:\n"
        )
        for candidate in checked:
            sys.stderr.write(f"  - {candidate}\n")
        sys.stderr.write(
            "Install dev dependencies in a project venv (for example: "
            'python -m pip install -e ".[dev]") and retry.\n'
        )
        return 1

    cmd = [
        python_bin,
        "-m",
        "pyright",
        "--warnings",
        "--pythonpath",
        python_bin,
        *sys.argv[1:],
    ]
    return subprocess.call(cmd)


if __name__ == "__main__":
    raise SystemExit(main())
