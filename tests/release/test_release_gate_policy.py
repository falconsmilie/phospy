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

AUTHORITATIVE_RELEASE_COMMAND = "make release-check"
PUBLIC_RELEASE_INSTRUCTION_DOCS = (
    Path("README.md"),
    Path("docs/maintenance.md"),
    Path("docs/testing/README.md"),
    Path("docs/testing/pytest_markers.md"),
    Path("docs/contributing.md"),
    Path(".github/CONTRIBUTING.md"),
)
RELEASE_CHECK_DEPENDENCIES = (
    "lint",
    "type-check",
    "test-unit",
    "test-parity",
    "test-performance",
    "validate-reference-bundles",
    "test-release-gates",
    "build",
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


def _make_target_line(target_name: str) -> str:
    lines = _read("Makefile").splitlines()
    return next(line for line in lines if line.startswith(f"{target_name}:"))


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
        if line.strip() or body:
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
    list_match = re.search(
        r"(?ms)^\s+needs:\s*\n(?P<items>(?:\s+- .+\n)+)",
        job_block,
    )
    if list_match is not None:
        return (
            re.search(
                rf"(?m)^\s+- {re.escape(dependency)}\s*$",
                list_match.group("items"),
            )
            is not None
        )
    return (
        re.search(rf"(?m)^\s+needs:\s*{re.escape(dependency)}\s*$", job_block)
        is not None
    )


def _assert_supported_python_matrix(job_block: str) -> None:
    assert "python-version: ['3.10', '3.11', '3.12']" in job_block


def test_default_pytest_keeps_parity_out_of_fast_local_loop() -> None:
    assert _pytest_config()["addopts"] == '-m "not parity"'


def test_release_check_command_is_maintained_project_authority() -> None:
    target_line = _make_target_line("release-check")
    target_dependencies = tuple(target_line.split(":", maxsplit=1)[1].split())

    assert AUTHORITATIVE_RELEASE_COMMAND == "make release-check"
    assert target_dependencies == RELEASE_CHECK_DEPENDENCIES
    assert _make_target_body("release-check") == ""


@pytest.mark.parametrize(
    "relative_path",
    PUBLIC_RELEASE_INSTRUCTION_DOCS,
    ids=lambda path: path.as_posix(),
)
def test_maintained_release_instructions_point_to_lightweight_release_check(
    relative_path: Path,
) -> None:
    text = _read(relative_path)
    normalized = re.sub(r"\s+", " ", text.lower())

    assert AUTHORITATIVE_RELEASE_COMMAND in text
    assert "normal ci/build confidence" in normalized
    assert "formal exact-source/exact-artifact attestation" in normalized


def test_make_release_check_covers_declared_blocking_marker_policy() -> None:
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

    makefile = _read("Makefile")
    assert re.search(r"(?ms)^\.PHONY:.*test-release-gates", makefile) is not None
    assert "make test-release-gates" in _make_target_body("help")

    assert '$(PYTEST) -m "not parity"' in _make_target_body("test-unit")
    assert (
        '$(PYTEST) tests/parity -m "parity and not parity_diagnostic" -s'
        in _make_target_body("test-parity")
    )
    assert 'tests/performance -m "performance or release_gate"' in _make_target_body(
        "test-performance"
    )
    release_gates = _make_target_body("test-release-gates")
    assert '$(MKDIR_P) "$(PYTEST_REPORT_DIR)"' in release_gates
    assert "$(PYTEST) -o addopts= tests/release tests/golden" in release_gates
    assert '-m "release_gate or golden or reproducibility"' in release_gates
    assert '--junitxml "$(PYTEST_REPORT_DIR)/release-gates.xml"' in release_gates
    assert (
        "$(PYTHON) scripts/validate_reference_bundle_index.py --repo-root ."
        in _make_target_body("validate-reference-bundles")
    )


def test_make_build_is_conventional_and_git_independent() -> None:
    build = _make_target_body("build")

    assert "$(RM_RF) dist" in build
    assert "$(BUILD)" in build
    assert "dist/*.whl" in build
    assert "dist/*.tar.gz" in build
    assert "$(TWINE) check dist/*" in build
    assert (
        "$(PYTHON) scripts/validate_reference_bundle_distribution.py "
        "--no-git-index-compare dist/*"
    ) in build
    assert build.index("$(RM_RF) dist") < build.index("$(BUILD)")
    assert build.index("$(BUILD)") < build.index("$(TWINE) check dist/*")
    assert build.index("$(TWINE) check dist/*") < build.index(
        "scripts/validate_reference_bundle_distribution.py"
    )


def test_publish_workflow_builds_once_and_publishes_uploaded_dist() -> None:
    workflow = _read(".github/workflows/publish.yml")
    build = _workflow_job_block(workflow, "build")
    testpypi = _workflow_job_block(workflow, "publish-to-testpypi")
    pypi = _workflow_job_block(workflow, "publish-to-pypi")

    assert "run: make release-check" in build
    assert "name: python-package-distributions" in build
    assert _has_needs(testpypi, "build")
    assert _has_needs(pypi, "build")
    assert "packages-dir: dist/" in testpypi
    assert "packages-dir: dist/" in pypi
    assert "id-token: write" in testpypi
    assert "id-token: write" in pypi


def test_ci_keeps_supported_python_source_tests_and_single_build_smoke() -> None:
    workflow = _read(".github/workflows/ci.yml")
    clean_install = _workflow_job_block(workflow, "clean-constrained-install")
    default_suite = _workflow_job_block(workflow, "default-suite")
    activity_parity = _workflow_job_block(workflow, "activity-parity-gate")
    hard_parity = _workflow_job_block(workflow, "parity-tests")
    diagnostics = _workflow_job_block(workflow, "parity-diagnostics")
    performance = _workflow_job_block(workflow, "performance-contracts")
    reference_bundles = _workflow_job_block(workflow, "reference-bundles")
    fixture_integrity = _workflow_job_block(workflow, "fixture-integrity")
    release_gates = _workflow_job_block(workflow, "release-gates")
    build = _workflow_job_block(workflow, "build-distributions")

    _assert_supported_python_matrix(clean_install)
    _assert_supported_python_matrix(default_suite)
    _assert_supported_python_matrix(activity_parity)
    _assert_supported_python_matrix(hard_parity)
    _assert_supported_python_matrix(performance)
    _assert_supported_python_matrix(release_gates)
    assert "python-version: '3.10'" in diagnostics
    assert "make test-performance" in performance
    assert "make validate-reference-bundles" in reference_bundles
    assert "runs-on: ${{ matrix.os }}" in fixture_integrity
    assert "os: [ubuntu-latest, windows-latest]" in fixture_integrity
    assert "python-version: '3.12'" in fixture_integrity
    assert (
        "test_manifest_fixture_byte_reproducibility.py::"
        "test_manifest_governed_fixtures_use_canonical_lf_bytes_and_valid_hashes"
    ) in fixture_integrity
    assert (
        "test_large_limma_trend_fixture_manifest_hashes_match_checked_in_files"
        in fixture_integrity
    )
    assert "make test-release-gates" in release_gates
    assert '-m "parity and activity_parity"' in activity_parity
    assert '-m "parity and not parity_diagnostic"' in hard_parity
    assert "continue-on-error: true" not in hard_parity
    assert "continue-on-error: true" in diagnostics
    assert '-m "parity_diagnostic"' in diagnostics
    assert "make build" in build
    assert "python -m venv" in build
    assert "-m pip check" in build
    assert "import pathlib, phospy" in build
