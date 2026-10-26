from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
PACKAGE_INITIALIZER = ROOT / "src" / "phospy" / "__init__.py"
METADATA_ABSENT = "__metadata_absent__"

_VERSION_PROBE = """
import importlib.metadata as metadata
import json
import sys

if sys.argv[2] == "__metadata_absent__":
    def missing_version(name):
        raise metadata.PackageNotFoundError(name)

    metadata.version = missing_version
else:
    metadata.version = lambda name: sys.argv[2]

sys.path.insert(0, sys.argv[1])
import phospy

print(json.dumps({"version": phospy.__version__}))
"""


def _copied_source_package(checkout: Path) -> Path:
    package_directory = checkout / "src" / "phospy"
    package_directory.mkdir(parents=True)
    shutil.copy2(PACKAGE_INITIALIZER, package_directory / "__init__.py")
    return checkout / "src"


def _probe_version(source_directory: Path, metadata_version: str) -> str:
    result = subprocess.run(
        [
            sys.executable,
            "-I",
            "-c",
            _VERSION_PROBE,
            str(source_directory),
            metadata_version,
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return str(json.loads(result.stdout)["version"])


def test_installed_distribution_metadata_precedes_checkout_version(
    tmp_path: Path,
) -> None:
    checkout = tmp_path / "checkout"
    source_directory = _copied_source_package(checkout)
    (checkout / "pyproject.toml").write_text(
        '[project]\nname = "phospy"\nversion = "9.8.7"\n',
        encoding="utf-8",
    )

    assert _probe_version(source_directory, "4.5.6") == "4.5.6"


def test_metadata_absent_source_checkout_uses_verified_project_file(
    tmp_path: Path,
) -> None:
    checkout = tmp_path / "checkout"
    source_directory = _copied_source_package(checkout)
    (checkout / "pyproject.toml").write_text(
        '[project]\nname = "phospy"\nversion = "9.8.7"\n',
        encoding="utf-8",
    )

    assert _probe_version(source_directory, METADATA_ABSENT) == "9.8.7"


@pytest.mark.parametrize(
    "pyproject_content",
    [
        None,
        '[project\nname = "phospy"\nversion = "9.8.7"\n',
        '[project]\nname = "another-project"\nversion = "9.8.7"\n',
        '[project]\nname = "phospy"\nversion = 987\n',
    ],
    ids=["missing", "malformed", "wrong-project", "non-string-version"],
)
def test_invalid_source_metadata_has_deterministic_unknown_version(
    tmp_path: Path,
    pyproject_content: str | None,
) -> None:
    checkout = tmp_path / "checkout"
    source_directory = _copied_source_package(checkout)
    if pyproject_content is not None:
        (checkout / "pyproject.toml").write_text(
            pyproject_content,
            encoding="utf-8",
        )

    assert _probe_version(source_directory, METADATA_ABSENT) == "0+unknown"


def test_unrelated_ancestor_project_file_is_never_consulted(tmp_path: Path) -> None:
    ancestor = tmp_path / "ancestor"
    ancestor.mkdir()
    (ancestor / "pyproject.toml").write_text(
        '[project]\nname = "phospy"\nversion = "7.7.7"\n',
        encoding="utf-8",
    )
    source_directory = _copied_source_package(ancestor / "nested" / "checkout")

    assert _probe_version(source_directory, METADATA_ABSENT) == "0+unknown"
