from __future__ import annotations

import json
import platform
import re
from datetime import datetime
from pathlib import Path
from typing import Any

import pytest

from phospy.release.metadata import (
    DEFAULT_RELEASE_GATE_METADATA_PATH,
    DEFAULT_TEST_COMMAND,
    DEFAULT_TEST_MARKERS,
    DEFAULT_TEST_STEPS,
    write_release_gate_metadata,
)

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 fallback
    import tomli as tomllib


ROOT = Path(__file__).resolve().parents[2]

pytestmark = pytest.mark.release_gate

RELEASE_PARITY_FILES = (
    Path("tests/parity/test_differential_analysis_parity.py"),
    Path("tests/parity/test_differential_limma_parity.py"),
    Path("tests/parity/test_kinase_workflow_parity.py"),
    Path("tests/parity/test_prediction_science_parity.py"),
    Path("tests/parity/test_l6_prediction_parity.py"),
    Path("tests/parity/test_public_predmat_parity.py"),
    Path("tests/parity/test_activity_stage_parity.py"),
    Path("tests/parity/test_signalome_workflow_parity.py"),
    Path("tests/parity/test_signalome_clustering_backend_parity.py"),
)


def _read(relative_path: str | Path) -> str:
    return (ROOT / relative_path).read_text(encoding="utf-8").replace("\r\n", "\n")


def _load_pyproject() -> dict[str, Any]:
    with (ROOT / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)


def _pytest_config() -> dict[str, Any]:
    pyproject = _load_pyproject()
    return pyproject["tool"]["pytest"]["ini_options"]


def _pytest_markers() -> set[str]:
    markers = _pytest_config()["markers"]
    return {marker.split(":", maxsplit=1)[0] for marker in markers}


def _make_target_body(target_name: str) -> str:
    lines = _read("Makefile").splitlines()
    start_index = next(
        index for index, line in enumerate(lines) if line.startswith(f"{target_name}:")
    )
    body: list[str] = []
    for line in lines[start_index + 1 :]:
        if line.startswith("\t"):
            body.append(line[1:])
            continue
        if body:
            break
    return "\n".join(body)


def _workflow_job_block(workflow_text: str, job_name: str) -> str:
    match = re.search(
        rf"(?ms)^  {re.escape(job_name)}:\n(?P<body>.*?)(?=^  [A-Za-z0-9_-]+:|\Z)",
        workflow_text,
    )
    assert match is not None, f"Workflow job is missing: {job_name}"
    return match.group("body")


def _has_needs(job_block: str, dependency: str) -> bool:
    return (
        re.search(
            rf"(?m)^\s+needs:\s*\n\s+- {re.escape(dependency)}\s*$",
            job_block,
        )
        is not None
    )


def test_default_pytest_keeps_parity_out_of_fast_local_loop() -> None:
    assert _pytest_config()["addopts"] == '-m "not parity"'


def test_publish_workflow_cannot_publish_without_scientific_release_gate() -> None:
    workflow = _read(".github/workflows/publish.yml")
    release_gate = _workflow_job_block(workflow, "release-gate")
    build = _workflow_job_block(workflow, "build")
    testpypi = _workflow_job_block(workflow, "publish-to-testpypi")
    pypi = _workflow_job_block(workflow, "publish-to-pypi")

    assert "run: make test-release-gate" in release_gate
    assert _has_needs(build, "release-gate")
    assert _has_needs(testpypi, "build")
    assert _has_needs(pypi, "build")


def test_make_release_gate_covers_declared_blocking_marker_policy() -> None:
    markers = _pytest_markers()
    assert markers == {
        "unit",
        "integration",
        "parity",
        "performance",
        "slow",
        "release_gate",
        "reproducibility",
        "golden",
        "activity_parity",
        "parity_diagnostic",
    }

    body = _make_target_body("test-release-gate")
    required_fragments = (
        '$(PYTEST) -m "not parity and not performance and not release_gate"',
        "tests/golden",
        "tests/unit/test_provenance_regressions.py",
        (
            "tests/integration/test_kinase_workflow_integration.py::"
            "test_kinase_public_predmat_provenance_matches_golden_contract"
        ),
        (
            "tests/integration/test_signalome_workflow_integration.py::"
            "test_signalome_l6_provenance_matches_golden_contract"
        ),
        '-m "release_gate and (reproducibility or golden)"',
        '$(PYTEST) tests/release -m "release_gate"',
        '$(PYTEST) tests/parity -m "parity and not parity_diagnostic" -s',
        '$(PYTEST) tests/performance -m "performance or release_gate" -q',
    )
    for fragment in required_fragments:
        assert fragment in body


def test_make_release_gate_writes_metadata_before_pytest_commands() -> None:
    body = _make_target_body("test-release-gate")
    body_lines = body.splitlines()

    metadata_fragment = "$(PYTHON) -m phospy.release.metadata"

    assert (
        f"RELEASE_GATE_METADATA_PATH ?= {DEFAULT_RELEASE_GATE_METADATA_PATH.as_posix()}"
    ) in _read("Makefile")
    assert metadata_fragment in body
    assert '--output "$(RELEASE_GATE_METADATA_PATH)"' in body
    metadata_index = next(
        index for index, line in enumerate(body_lines) if metadata_fragment in line
    )
    first_pytest_index = next(
        index
        for index, line in enumerate(body_lines)
        if line.startswith(
            '$(PYTEST) -m "not parity and not performance and not release_gate"'
        )
    )
    assert metadata_index < first_pytest_index


def test_release_gate_writes_metadata_artifact(tmp_path: Path) -> None:
    metadata_path = tmp_path / "release_gate_metadata.json"

    written_path = write_release_gate_metadata(
        metadata_path,
        project_root=ROOT,
        generated_at_utc="2026-07-02T00:00:00Z",
    )

    assert written_path == metadata_path
    assert metadata_path.is_file()

    payload = json.loads(metadata_path.read_text(encoding="utf-8"))
    serialized_payload = json.dumps(payload)
    assert {
        "phospy_version",
        "python_version",
        "platform",
        "dependency_snapshot",
        "test_command",
        "test_markers",
        "parity_fixture_versions",
        "generated_at_utc",
    }.issubset(payload)
    assert str(ROOT) not in serialized_payload
    assert ROOT.as_posix() not in serialized_payload


def test_release_gate_metadata_contains_runtime_and_package_versions(
    tmp_path: Path,
) -> None:
    metadata_path = tmp_path / "release_gate_metadata.json"

    write_release_gate_metadata(
        metadata_path,
        project_root=ROOT,
        generated_at_utc="2026-07-02T00:00:00Z",
    )

    payload = json.loads(metadata_path.read_text(encoding="utf-8"))
    generated_at = datetime.fromisoformat(
        payload["generated_at_utc"].replace("Z", "+00:00")
    )

    assert payload["phospy_version"] == _load_pyproject()["project"]["version"]
    assert payload["python_version"] == platform.python_version()
    assert payload["platform"] == platform.platform()
    assert generated_at.tzinfo is not None
    assert payload["test_command"] == DEFAULT_TEST_COMMAND
    assert payload["test_markers"] == list(DEFAULT_TEST_MARKERS)
    assert payload["test_steps"] == list(DEFAULT_TEST_STEPS)

    dependency_snapshot = payload["dependency_snapshot"]
    assert {
        "numpy",
        "pandas",
        "scipy",
        "scikit-learn",
        "pytest",
        "hypothesis",
    }.issubset(dependency_snapshot)
    assert all(
        version is None or isinstance(version, str)
        for version in dependency_snapshot.values()
    )

    parity_fixture_versions = payload["parity_fixture_versions"]
    assert isinstance(parity_fixture_versions, dict)
    assert "tests/fixtures/rewrite_parity/differential_r_reference" in (
        parity_fixture_versions
    )
    assert "tests/fixtures/public_workflow_reference" in parity_fixture_versions
    assert "src/phospy/data/reference_bundles/rat/l6_native" in (
        parity_fixture_versions
    )
    assert (
        parity_fixture_versions[
            "tests/fixtures/rewrite_parity/differential_r_reference"
        ]["declared_versions"]["limma"]
        == "3.66.0"
    )


def test_all_current_release_gate_marked_files_are_collected_by_make_gate() -> None:
    body = _make_target_body("test-release-gate")
    covered_roots = ("tests/golden/", "tests/performance/", "tests/release/")
    explicitly_collected_files = {
        "tests/unit/test_provenance_regressions.py",
        "tests/integration/test_kinase_workflow_integration.py",
        "tests/integration/test_signalome_workflow_integration.py",
    }

    release_gate_files = {
        path.relative_to(ROOT).as_posix()
        for path in (ROOT / "tests").rglob("test_*.py")
        if "release_gate" in path.read_text(encoding="utf-8")
    }
    assert release_gate_files

    for relative_path in sorted(release_gate_files):
        covered_by_root = relative_path.startswith(covered_roots)
        covered_explicitly = (
            relative_path in explicitly_collected_files and relative_path in body
        )
        assert covered_by_root or covered_explicitly, relative_path


def test_ci_keeps_diagnostic_parity_non_blocking() -> None:
    workflow = _read(".github/workflows/ci.yml")
    hard_parity = _workflow_job_block(workflow, "parity-tests")
    diagnostics = _workflow_job_block(workflow, "parity-diagnostics")

    assert '-m "parity and not parity_diagnostic"' in hard_parity
    assert "continue-on-error: true" not in hard_parity
    assert "continue-on-error: true" in diagnostics
    assert '-m "parity_diagnostic"' in diagnostics


@pytest.mark.parametrize(
    "relative_path",
    RELEASE_PARITY_FILES,
    ids=lambda path: path.as_posix(),
)
def test_release_gate_includes_required_threshold_bearing_parity_lanes(
    relative_path: Path,
) -> None:
    body = _make_target_body("test-release-gate")
    source = _read(relative_path)

    assert '$(PYTEST) tests/parity -m "parity and not parity_diagnostic" -s' in body
    assert relative_path.as_posix().startswith("tests/parity/")
    assert "pytest.mark.parity" in source
