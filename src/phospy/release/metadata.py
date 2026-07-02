"""Machine-readable metadata for release-gate runs."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import os
import platform as runtime_platform
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType
from typing import Any, cast

from phospy.provenance.environment import (
    DEFAULT_ENVIRONMENT_DEPENDENCIES,
    collect_environment_provenance,
)

DEFAULT_RELEASE_GATE_METADATA_PATH = Path(
    "build/release-gate/release_gate_metadata.json"
)
DEFAULT_TEST_COMMAND = "make test-release-gate"
DEFAULT_TEST_MARKERS = (
    "not parity and not performance and not release_gate",
    "release_gate and (reproducibility or golden)",
    "release_gate",
    "parity and not parity_diagnostic",
    "performance or release_gate",
)
DEFAULT_TEST_STEPS = (
    'python -m pytest -m "not parity and not performance and not release_gate"',
    (
        "python -m pytest tests/golden tests/unit/test_provenance_regressions.py "
        "tests/integration/test_kinase_workflow_integration.py::"
        "test_kinase_public_predmat_provenance_matches_golden_contract "
        "tests/integration/test_signalome_workflow_integration.py::"
        "test_signalome_l6_provenance_matches_golden_contract "
        '-m "release_gate and (reproducibility or golden)"'
    ),
    'python -m pytest tests/release -m "release_gate"',
    'python -m pytest tests/parity -m "parity and not parity_diagnostic" -s',
    'python -m pytest tests/performance -m "performance or release_gate" -q',
)
RELEASE_GATE_DEPENDENCIES = (
    *DEFAULT_ENVIRONMENT_DEPENDENCIES,
    "pytest",
    "hypothesis",
)
METADATA_SCHEMA_VERSION = 1


def build_release_gate_metadata(
    *,
    test_command: str = DEFAULT_TEST_COMMAND,
    test_markers: Sequence[str] = DEFAULT_TEST_MARKERS,
    test_steps: Sequence[str] = DEFAULT_TEST_STEPS,
    project_root: Path | None = None,
    generated_at_utc: str | None = None,
) -> dict[str, Any]:
    """Build the release-gate metadata payload without writing it."""

    root = project_root or _project_root_from_module()
    environment = collect_environment_provenance(
        dependency_names=_unique_dependencies(RELEASE_GATE_DEPENDENCIES),
        use_cache=False,
    )
    return {
        "metadata_schema_version": METADATA_SCHEMA_VERSION,
        "phospy_version": _project_metadata_version(root)
        or environment.package_version,
        "python_version": environment.python_version,
        "platform": runtime_platform.platform(),
        "dependency_snapshot": dict(environment.dependency_versions),
        "test_command": test_command,
        "test_markers": list(test_markers),
        "test_steps": list(test_steps),
        "parity_fixture_versions": _collect_parity_fixture_versions(root),
        "generated_at_utc": generated_at_utc or _utc_timestamp(),
    }


def write_release_gate_metadata(
    output_path: str | Path = DEFAULT_RELEASE_GATE_METADATA_PATH,
    *,
    test_command: str = DEFAULT_TEST_COMMAND,
    test_markers: Sequence[str] = DEFAULT_TEST_MARKERS,
    test_steps: Sequence[str] = DEFAULT_TEST_STEPS,
    project_root: Path | None = None,
    generated_at_utc: str | None = None,
) -> Path:
    """Write the release-gate metadata JSON artifact and return its path."""

    path = Path(output_path)
    payload = build_release_gate_metadata(
        test_command=test_command,
        test_markers=test_markers,
        test_steps=test_steps,
        project_root=project_root,
        generated_at_utc=generated_at_utc,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Write release-gate metadata as a JSON artifact."
    )
    parser.add_argument(
        "--output",
        default=os.environ.get(
            "RELEASE_GATE_METADATA_PATH",
            DEFAULT_RELEASE_GATE_METADATA_PATH.as_posix(),
        ),
        help="Output JSON path.",
    )
    parser.add_argument(
        "--test-command",
        default=DEFAULT_TEST_COMMAND,
        help="Top-level release-gate command that was run.",
    )
    parser.add_argument(
        "--test-marker",
        action="append",
        dest="test_markers",
        help="Pytest marker expression included in the release gate.",
    )
    parser.add_argument(
        "--test-step",
        action="append",
        dest="test_steps",
        help="Concrete test command included in the release gate.",
    )
    args = parser.parse_args(argv)

    path = write_release_gate_metadata(
        args.output,
        test_command=args.test_command,
        test_markers=tuple(args.test_markers or DEFAULT_TEST_MARKERS),
        test_steps=tuple(args.test_steps or DEFAULT_TEST_STEPS),
    )
    print(path.as_posix())
    return 0


def _collect_parity_fixture_versions(project_root: Path) -> dict[str, Any]:
    fixture_versions: dict[str, Any] = {}
    rewrite_root = project_root / "tests" / "fixtures" / "rewrite_parity"
    public_workflow_root = (
        project_root / "tests" / "fixtures" / "public_workflow_reference"
    )
    reference_bundle_root = (
        project_root / "src" / "phospy" / "data" / "reference_bundles"
    )

    if rewrite_root.is_dir():
        for provenance_file in sorted(rewrite_root.rglob("PROVENANCE.md")):
            directory = provenance_file.parent
            fixture_versions[_relative_path(directory, project_root)] = {
                "content_sha256": _directory_sha256(directory),
                "file_count": _directory_file_count(directory),
                "provenance_sha256": _file_sha256(provenance_file),
                "declared_versions": _extract_declared_versions(
                    provenance_file.read_text(encoding="utf-8")
                ),
            }

    if public_workflow_root.is_dir():
        fixture_versions[_relative_path(public_workflow_root, project_root)] = {
            "content_sha256": _directory_sha256(public_workflow_root),
            "file_count": _directory_file_count(public_workflow_root),
            "contracts": _public_workflow_contract_versions(public_workflow_root),
            "provenance_sha256": _optional_file_sha256(
                public_workflow_root / "PROVENANCE.md"
            ),
        }

    if reference_bundle_root.is_dir():
        for manifest_file in sorted(reference_bundle_root.rglob("manifest.json")):
            directory = manifest_file.parent
            fixture_versions[_relative_path(directory, project_root)] = {
                "content_sha256": _directory_sha256(directory),
                "file_count": _directory_file_count(directory),
                "manifest": _reference_manifest_version(manifest_file),
            }

    return fixture_versions


def _public_workflow_contract_versions(directory: Path) -> dict[str, Any]:
    contracts: dict[str, Any] = {}
    for contract_file in sorted(directory.glob("*_contract.json")):
        with contract_file.open("rb") as handle:
            payload = json.load(handle)
        contracts[contract_file.name] = {
            "sha256": _file_sha256(contract_file),
            "fixture_set_id": _string_or_none(payload.get("fixture_set_id")),
            "generation_date": _string_or_none(payload.get("generation_date")),
        }
    return contracts


def _reference_manifest_version(manifest_file: Path) -> dict[str, Any]:
    with manifest_file.open("rb") as handle:
        payload = json.load(handle)
    return {
        "sha256": _file_sha256(manifest_file),
        "reference_id": _string_or_none(payload.get("reference_id")),
        "reference_version": _string_or_none(payload.get("reference_version")),
        "source_name": _string_or_none(payload.get("source_name")),
        "source_version": _string_or_none(payload.get("source_version")),
        "generated_at_utc": _string_or_none(payload.get("generated_at_utc")),
        "manifest_schema_version": _string_or_none(
            payload.get("manifest_schema_version")
        ),
    }


def _extract_declared_versions(provenance_text: str) -> dict[str, str]:
    declared_versions: dict[str, str] = {}
    for line in provenance_text.splitlines():
        stripped = line.strip()
        lowered = stripped.lower()
        if lowered.startswith("generated with r version"):
            declared_versions["r"] = stripped.removeprefix(
                "Generated with R version"
            ).strip()
        elif lowered.startswith("limma version:"):
            declared_versions["limma"] = stripped.split(":", maxsplit=1)[1].strip()
        elif lowered.startswith("seed:"):
            declared_versions["seed"] = stripped.split(":", maxsplit=1)[1].strip()
    return declared_versions


def _directory_sha256(directory: Path) -> str:
    hasher = hashlib.sha256()
    for file_path in _directory_files(directory):
        hasher.update(file_path.relative_to(directory).as_posix().encode("utf-8"))
        hasher.update(b"\0")
        hasher.update(file_path.read_bytes())
        hasher.update(b"\0")
    return hasher.hexdigest()


def _directory_file_count(directory: Path) -> int:
    return sum(1 for _ in _directory_files(directory))


def _directory_files(directory: Path) -> tuple[Path, ...]:
    return tuple(sorted(path for path in directory.rglob("*") if path.is_file()))


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _optional_file_sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    return _file_sha256(path)


def _relative_path(path: Path, project_root: Path) -> str:
    return path.relative_to(project_root).as_posix()


def _string_or_none(value: object) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip()
    if not normalized:
        return None
    return normalized


def _project_metadata_version(project_root: Path) -> str | None:
    toml_parser = _toml_parser()
    if toml_parser is None:
        return None
    try:
        payload = cast(
            Mapping[str, object],
            toml_parser.loads(
                (project_root / "pyproject.toml").read_text(encoding="utf-8")
            ),
        )
    except Exception:
        return None
    project = payload.get("project")
    if not isinstance(project, Mapping):
        return None
    name = _string_or_none(project.get("name"))
    if name != "phospy":
        return None
    return _string_or_none(project.get("version"))


def _toml_parser() -> ModuleType | None:
    for module_name in ("tomllib", "tomli"):
        try:
            return importlib.import_module(module_name)
        except ModuleNotFoundError:
            continue
    return None


def _unique_dependencies(dependencies: Sequence[str]) -> tuple[str, ...]:
    unique: list[str] = []
    seen: set[str] = set()
    for dependency in dependencies:
        normalized = dependency.strip()
        if not normalized or normalized in seen:
            continue
        unique.append(normalized)
        seen.add(normalized)
    return tuple(unique)


def _utc_timestamp() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _project_root_from_module() -> Path:
    return Path(__file__).resolve().parents[3]


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


__all__ = [
    "DEFAULT_RELEASE_GATE_METADATA_PATH",
    "DEFAULT_TEST_COMMAND",
    "DEFAULT_TEST_MARKERS",
    "DEFAULT_TEST_STEPS",
    "METADATA_SCHEMA_VERSION",
    "build_release_gate_metadata",
    "main",
    "write_release_gate_metadata",
]
