from __future__ import annotations

import json
import platform
import re
from datetime import datetime
from pathlib import Path
from typing import Any

import pytest

from phospy.release.metadata import (
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

AUTHORITATIVE_RELEASE_GATE_COMMAND = "make test-release-gate"
PUBLIC_RELEASE_INSTRUCTION_DOCS = (
    Path("README.md"),
    Path("docs/maintenance.md"),
    Path("docs/testing/README.md"),
    Path("docs/testing/pytest_markers.md"),
    Path("docs/contributing.md"),
    Path(".github/CONTRIBUTING.md"),
    Path("docs/scientific-coverage.md"),
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


def test_release_gate_command_is_maintained_project_authority() -> None:
    body = _make_target_body("test-release-gate")

    assert DEFAULT_TEST_COMMAND == AUTHORITATIVE_RELEASE_GATE_COMMAND
    assert body


def test_pyright_requirement_matches_ci_constraint() -> None:
    pyproject = _load_pyproject()
    dev_requirements = pyproject["project"]["optional-dependencies"]["dev"]
    pyright_requirement = next(
        requirement
        for requirement in dev_requirements
        if requirement.startswith("pyright>=")
    )
    constraint_text = _read("constraints/ci.txt")

    required_version = re.search(r"pyright>=(\d+\.\d+\.\d+)", pyright_requirement)
    constrained_version = re.search(
        r"(?m)^pyright==(\d+\.\d+\.\d+)$",
        constraint_text,
    )

    assert required_version is not None
    assert constrained_version is not None
    assert constrained_version.group(1) == required_version.group(1)


@pytest.mark.parametrize(
    "relative_path",
    PUBLIC_RELEASE_INSTRUCTION_DOCS,
    ids=lambda path: path.as_posix(),
)
def test_maintained_release_instructions_require_full_gate(
    relative_path: Path,
) -> None:
    text = _read(relative_path)
    normalized = re.sub(r"\s+", " ", text.lower())
    normalized_without_markdown = normalized.replace("`", "")

    assert AUTHORITATIVE_RELEASE_GATE_COMMAND in text
    assert "default pytest" in normalized_without_markdown
    assert "not sufficient for release" in normalized
    assert "release tests" in normalized
    assert "parity" in normalized
    assert "golden" in normalized or "reproducibility" in normalized
    assert "performance" in normalized


def test_publish_workflow_cannot_publish_without_scientific_release_gate() -> None:
    workflow = _read(".github/workflows/publish.yml")
    release_gate = _workflow_job_block(workflow, "release-gate")
    build = _workflow_job_block(workflow, "build")
    distribution_install = _workflow_job_block(workflow, "distribution-install-tests")
    release_attestation = _workflow_job_block(workflow, "release-attestation")
    testpypi = _workflow_job_block(workflow, "publish-to-testpypi")
    pypi = _workflow_job_block(workflow, "publish-to-pypi")

    assert "run: make test-release-gate" in release_gate
    _assert_supported_python_matrix(release_gate)
    assert "release-gate-py${{ matrix.python-version }}" in release_gate
    assert _has_needs(build, "release-gate")
    assert "make build" in build
    assert build.index("make build") < build.index("twine check dist/*")
    assert build.index("twine check dist/*") < build.index(
        "uses: actions/upload-artifact@v5"
    )
    assert _has_needs(distribution_install, "build")
    _assert_supported_python_matrix(distribution_install)
    assert '"${wheel[0]}"' in distribution_install
    assert '"${wheel[0]}[test]"' not in distribution_install
    assert '"${sdist[0]}"' in distribution_install
    assert '"${sdist[0]}[test]"' not in distribution_install
    assert 'python -I "$GITHUB_WORKSPACE/scripts/verify_distribution_artifact.py"' in (
        distribution_install
    )
    assert "--artifact-kind wheel" in distribution_install
    assert "--artifact-kind sdist" in distribution_install
    assert "--artifact-path" in distribution_install
    assert "--build-manifest" in distribution_install
    assert "python-package-build-manifest" in distribution_install
    assert "$GITHUB_WORKSPACE/build/reports/build-manifest.json" in (
        distribution_install
    )
    assert "tests/distribution" not in distribution_install
    assert "python -m pytest" not in distribution_install
    assert _has_needs(release_attestation, "distribution-install-tests")
    assert "scripts/release/build_release_attestation.py" in release_attestation
    assert "scripts/release/verify_release_attestation.py" in release_attestation
    assert "release-attestation.json" in release_attestation
    assert _has_needs(testpypi, "release-attestation")
    assert _has_needs(pypi, "release-attestation")
    assert "scripts/release/prepare_publish_directory.py" in testpypi
    assert "packages-dir: build/publish/" in testpypi
    assert "scripts/release/prepare_publish_directory.py" in pypi
    assert "packages-dir: build/publish/" in pypi


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
        "$(PYTHON) scripts/validate_reference_bundle_index.py",
        (
            '$(PYTEST) $(PYTEST_DURATION_ARGS) -m "not parity and not '
            'performance and not release_gate"'
        ),
        '--junitxml "$(PYTEST_REPORT_DIR)/release-default.xml"',
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
        '$(PYTEST) $(PYTEST_DURATION_ARGS) tests/release -m "release_gate"',
        (
            '$(PYTEST) $(PYTEST_DURATION_ARGS) tests/parity -m "parity and '
            'not parity_diagnostic" -s'
        ),
        (
            "$(PYTEST) $(PYTEST_DURATION_ARGS) tests/performance -m "
            '"performance or release_gate" -q'
        ),
        '--junitxml "$(PYTEST_REPORT_DIR)/release-performance.xml"',
    )
    for fragment in required_fragments:
        assert fragment in body


def test_make_release_gate_writes_source_identity_before_source_checks() -> None:
    makefile = _read("Makefile")
    body = _make_target_body("test-release-gate")
    body_lines = body.splitlines()

    source_identity_fragment = "$(PYTHON) scripts/release/create_source_identity.py"

    assert "SOURCE_IDENTITY_PATH ?= build/reports/source-identity.json" in makefile
    assert source_identity_fragment in body
    assert '--output "$(SOURCE_IDENTITY_PATH)"' in body
    assert "scripts/release/write_source_check_report.py" in body
    source_identity_index = next(
        index
        for index, line in enumerate(body_lines)
        if source_identity_fragment in line
    )
    first_pytest_index = next(
        index
        for index, line in enumerate(body_lines)
        if line.startswith(
            '$(PYTEST) $(PYTEST_DURATION_ARGS) -m "not parity and not '
            'performance and not release_gate"'
        )
    )
    first_source_report_index = next(
        index
        for index, line in enumerate(body_lines)
        if "scripts/release/write_source_check_report.py" in line
    )
    assert source_identity_index < first_pytest_index
    assert first_pytest_index < first_source_report_index


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
        if "pytest.mark.release_gate" in path.read_text(encoding="utf-8")
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


def test_ci_release_verdict_requires_supported_version_verification() -> None:
    workflow = _read(".github/workflows/ci.yml")
    clean_install = _workflow_job_block(workflow, "clean-constrained-install")
    default_suite = _workflow_job_block(workflow, "default-suite")
    release_gate = _workflow_job_block(workflow, "release-gate")
    distribution_install = _workflow_job_block(workflow, "distribution-install-tests")
    verdict = _workflow_job_block(workflow, "release-verdict")

    _assert_supported_python_matrix(clean_install)
    _assert_supported_python_matrix(default_suite)
    _assert_supported_python_matrix(release_gate)
    _assert_supported_python_matrix(distribution_install)

    assert 'python -m pip install -e ".[dev,test]"' in clean_install
    assert 'pytest --durations=25 --durations-min=0.01 -m "not parity"' in (
        default_suite
    )
    assert "run: make test-release-gate" in release_gate
    assert "release-gate-py${{ matrix.python-version }}" in release_gate
    assert '"${wheel[0]}"' in distribution_install
    assert '"${wheel[0]}[test]"' not in distribution_install
    assert '"${sdist[0]}"' in distribution_install
    assert '"${sdist[0]}[test]"' not in distribution_install
    assert 'python -I "$GITHUB_WORKSPACE/scripts/verify_distribution_artifact.py"' in (
        distribution_install
    )
    assert "--artifact-path" in distribution_install
    assert "--build-manifest" in distribution_install
    assert "python-package-build-manifest" in distribution_install
    assert "$GITHUB_WORKSPACE/build/reports/build-manifest.json" in (
        distribution_install
    )
    assert "tests/distribution" not in distribution_install
    assert "python -m pytest" not in distribution_install
    for dependency in (
        "clean-constrained-install",
        "default-suite",
        "performance-contracts",
        "release-gate",
        "build-distributions",
        "distribution-install-tests",
    ):
        assert _has_needs(verdict, dependency)


def test_ci_distribution_build_validates_reference_bundle_archives() -> None:
    workflow = _read(".github/workflows/ci.yml")
    build = _workflow_job_block(workflow, "build-distributions")
    makefile = _read("Makefile")
    make_build = _make_target_body("build")
    validation_command = (
        "python scripts/validate_reference_bundle_distribution.py dist/*"
    )

    assert "make build" in build
    assert build.index("make build") < build.index("twine check dist/*")
    assert "build/reports/build-manifest.json" in build
    assert "build: check-tools" in makefile
    assert "$(PYTHON) scripts/release/create_source_identity.py" in make_build
    assert "$(PYTHON) scripts/validate_reference_bundle_index.py" in make_build
    assert "$(BUILD)" in make_build
    assert "$(PYTHON) scripts/write_build_manifest.py" in make_build
    assert "$(PYTHON) scripts/validate_reference_bundle_distribution.py dist/*" in (
        make_build
    )
    assert validation_command.replace("python", "$(PYTHON)") in make_build


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

    assert (
        '$(PYTEST) $(PYTEST_DURATION_ARGS) tests/parity -m "parity and '
        'not parity_diagnostic" -s'
    ) in body
    assert relative_path.as_posix().startswith("tests/parity/")
    assert "pytest.mark.parity" in source
