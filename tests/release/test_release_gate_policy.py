from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pytest

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
