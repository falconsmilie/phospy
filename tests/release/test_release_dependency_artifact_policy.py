from __future__ import annotations

import re
import shlex
from pathlib import Path
from typing import Any

import pytest

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 fallback
    import tomli as tomllib


ROOT = Path(__file__).resolve().parents[2]
SUPPORTED_PYTHON_VERSIONS = ("3.10", "3.11", "3.12")
WORKFLOW_PATHS = (
    Path(".github/workflows/ci.yml"),
    Path(".github/workflows/publish.yml"),
)
ARTIFACT_VERIFIER = Path("scripts/verify_distribution_artifact.py")
BUILD_MANIFEST_WRITER = Path("scripts/write_build_manifest.py")
EXPLICITLY_JUSTIFIED_UNPINNED_DIRECT_INSTALLS = {
    "pip": "bootstrap installer upgraded before constraint-enforced installs",
}

pytestmark = pytest.mark.release_gate


def _read(relative_path: str | Path) -> str:
    return (ROOT / relative_path).read_text(encoding="utf-8").replace("\r\n", "\n")


def _load_pyproject() -> dict[str, Any]:
    with (ROOT / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)


def _normalise_distribution_name(value: str) -> str:
    return value.strip().lower().replace("_", "-")


def _requirement_name(requirement: str) -> str:
    requirement = requirement.split(";", maxsplit=1)[0].strip()
    requirement = requirement.split("[", maxsplit=1)[0].strip()
    return _normalise_distribution_name(
        re.split(r"\s*(?:===|==|~=|!=|<=|>=|<|>)", requirement, maxsplit=1)[0]
    )


def _pinned_constraint_names() -> set[str]:
    names: set[str] = set()
    for raw_line in _read("constraints/ci.txt").splitlines():
        line = raw_line.split("#", maxsplit=1)[0].strip()
        if not line:
            continue
        requirement = line.split(";", maxsplit=1)[0].strip()
        if "==" not in requirement:
            continue
        names.add(_requirement_name(requirement))
    return names


def _pyproject_release_dependency_names() -> set[str]:
    pyproject = _load_pyproject()
    project = pyproject["project"]
    names = {_requirement_name(item) for item in project["dependencies"]}
    extras = project["optional-dependencies"]
    for extra_name in ("dev", "test", "parquet"):
        names.update(_requirement_name(item) for item in extras[extra_name])
    names.update(
        _requirement_name(item) for item in pyproject["build-system"]["requires"]
    )
    return names


def _workflow_direct_install_names(workflow_text: str) -> set[str]:
    names: set[str] = set()
    for match in re.finditer(
        r"(?m)^\s*(?:python -m )?pip install (?P<arguments>.+)$",
        workflow_text,
    ):
        tokens = shlex.split(match.group("arguments"))
        skip_next = False
        for token in tokens:
            if skip_next:
                skip_next = False
                continue
            if token in {"-c", "--constraint"}:
                skip_next = True
                continue
            if token.startswith("-"):
                continue
            if token.startswith(".") or token.startswith("$") or "/" in token:
                continue
            names.add(_requirement_name(token))
    return names


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


def _expected_matrix_literal() -> str:
    quoted = ", ".join(f"'{version}'" for version in SUPPORTED_PYTHON_VERSIONS)
    return f"python-version: [{quoted}]"


def test_ci_constraints_pin_every_direct_release_dependency() -> None:
    direct_dependencies = _pyproject_release_dependency_names()
    for workflow_path in WORKFLOW_PATHS:
        direct_dependencies.update(_workflow_direct_install_names(_read(workflow_path)))

    unpinned = sorted(
        direct_dependencies
        - _pinned_constraint_names()
        - set(EXPLICITLY_JUSTIFIED_UNPINNED_DIRECT_INSTALLS)
    )

    assert unpinned == []


def test_ci_and_publish_install_commands_share_the_ci_constraints() -> None:
    for workflow_path in WORKFLOW_PATHS:
        workflow = _read(workflow_path)
        assert "PIP_CONSTRAINT: ${{ github.workspace }}/constraints/ci.txt" in workflow
        assert (
            "PIP_BUILD_CONSTRAINT: ${{ github.workspace }}/constraints/ci.txt"
            in workflow
        )


def test_distribution_build_uses_preinstalled_constrained_pep517_dependencies() -> None:
    makefile = _read("Makefile")
    assert "BUILD ?= $(PYTHON) -m build --no-isolation" in makefile
    assert "$(PYTHON) scripts/write_build_manifest.py" in makefile

    for workflow_path in WORKFLOW_PATHS:
        workflow = _read(workflow_path)
        build_job = _workflow_job_block(
            workflow,
            "build" if workflow_path.name == "publish.yml" else "build-distributions",
        )
        assert (
            "python -m pip install -c constraints/ci.txt build twine setuptools wheel"
            in build_job
        )
        assert "make build" in build_job
        assert "build/reports/build-manifest.json" in build_job


def test_python_support_declaration_matches_release_workflow_matrices() -> None:
    pyproject = _load_pyproject()
    classifiers = set(pyproject["project"]["classifiers"])
    assert pyproject["project"]["requires-python"] == ">=3.10,<3.13"
    assert {
        classifier.rsplit("::", maxsplit=1)[1].strip()
        for classifier in classifiers
        if classifier.startswith("Programming Language :: Python :: 3.")
    } == set(SUPPORTED_PYTHON_VERSIONS)

    matrix_literal = _expected_matrix_literal()
    for workflow_path in WORKFLOW_PATHS:
        workflow = _read(workflow_path)
        release_gate = _workflow_job_block(workflow, "release-gate")
        distribution_install = _workflow_job_block(
            workflow,
            "distribution-install-tests",
        )
        assert matrix_literal in release_gate
        assert matrix_literal in distribution_install


def test_installed_artifact_verification_is_standalone_and_outside_tests() -> None:
    verifier_path = ROOT / ARTIFACT_VERIFIER
    verifier_source = verifier_path.read_text(encoding="utf-8")
    manifest_writer_path = ROOT / BUILD_MANIFEST_WRITER

    assert verifier_path.is_file()
    assert manifest_writer_path.is_file()
    assert ARTIFACT_VERIFIER.parts[0] != "tests"
    assert BUILD_MANIFEST_WRITER.parts[0] != "tests"
    assert not (ROOT / "tests" / "distribution").exists()
    assert not re.search(
        r"(?m)^\s*(?:import pytest|from pytest import )", verifier_source
    )
    assert not re.search(
        r"(?m)^\s*(?:import tests|from tests(?:\.|\s))", verifier_source
    )
    assert "sys.path" not in verifier_source


def test_artifact_install_jobs_run_standalone_verifier_for_wheel_and_sdist() -> None:
    verifier_command = (
        'python -I "$GITHUB_WORKSPACE/scripts/verify_distribution_artifact.py"'
    )
    for workflow_path in WORKFLOW_PATHS:
        workflow = _read(workflow_path)
        distribution_install = _workflow_job_block(
            workflow,
            "distribution-install-tests",
        )

        assert "mktemp -d" in distribution_install
        assert 'cd "$verify_dir"' in distribution_install
        assert distribution_install.count(verifier_command) == 2
        assert "--artifact-kind wheel" in distribution_install
        assert "--artifact-kind sdist" in distribution_install
        assert "--artifact-path" in distribution_install
        assert "--build-manifest" in distribution_install
        assert '--repository-root "$GITHUB_WORKSPACE"' in distribution_install
        assert "python-package-build-manifest" in distribution_install
        assert "$GITHUB_WORKSPACE/build/reports/build-manifest.json" in (
            distribution_install
        )
        assert (
            "wheel-artifact-py${{ matrix.python-version }}.json" in distribution_install
        )
        assert (
            "sdist-artifact-py${{ matrix.python-version }}.json" in distribution_install
        )
        assert "tests/distribution" not in distribution_install
        assert "python -m pytest" not in distribution_install
        assert '"${wheel[0]}"' in distribution_install
        assert '"${wheel[0]}[test]"' not in distribution_install
        assert '"${sdist[0]}"' in distribution_install
        assert '"${sdist[0]}[test]"' not in distribution_install
        assert (
            'python -m pip install -c "$GITHUB_WORKSPACE/constraints/ci.txt" '
            "setuptools wheel"
        ) in distribution_install
        assert (
            "python -m pip install --no-build-isolation -c "
            '"$GITHUB_WORKSPACE/constraints/ci.txt" "${sdist[0]}"'
        ) in distribution_install
        assert distribution_install.index('"${wheel[0]}"') < (
            distribution_install.index("--artifact-kind wheel")
        )
        assert distribution_install.index('"${sdist[0]}"') < (
            distribution_install.index("--artifact-kind sdist")
        )


def test_distribution_artifact_failures_block_publish_jobs() -> None:
    workflow = _read(".github/workflows/publish.yml")
    distribution_install = _workflow_job_block(workflow, "distribution-install-tests")
    release_attestation = _workflow_job_block(workflow, "release-attestation")
    assert "continue-on-error" not in distribution_install
    assert 'python -I "$GITHUB_WORKSPACE/scripts/verify_distribution_artifact.py"' in (
        distribution_install
    )
    assert "tests/distribution" not in distribution_install

    assert _has_needs(release_attestation, "distribution-install-tests")
    assert "scripts/release/build_release_attestation.py" in release_attestation
    assert "scripts/release/verify_release_attestation.py" in release_attestation

    for publish_job in ("publish-to-testpypi", "publish-to-pypi"):
        publish_block = _workflow_job_block(workflow, publish_job)
        assert _has_needs(publish_block, "release-attestation")
        assert "scripts/release/prepare_publish_directory.py" in publish_block
        assert "packages-dir: build/publish/" in publish_block
