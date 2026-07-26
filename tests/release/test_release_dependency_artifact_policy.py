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
MINIMUM_CONSTRAINT_LOWER_BOUND_PINS = {
    "numpy": "1.24.0",
    "pandas": "2.0.0",
    "scipy": "1.10.0",
    "scikit-learn": "1.6.0",
    "tomli": "2.0.0",
    "hypothesis": "6.0.0",
    "pytest": "8.0.0",
    "setuptools": "77.0.3",
    "wheel": "0.45.1",
}
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


def _pinned_constraint_names(
    constraint_path: str | Path = "constraints/ci.txt",
) -> set[str]:
    names: set[str] = set()
    for raw_line in _read(constraint_path).splitlines():
        line = raw_line.split("#", maxsplit=1)[0].strip()
        if not line:
            continue
        requirement = line.split(";", maxsplit=1)[0].strip()
        if "==" not in requirement:
            continue
        names.add(_requirement_name(requirement))
    return names


def _pinned_constraint_versions(
    constraint_path: str | Path,
) -> dict[str, str]:
    versions: dict[str, str] = {}
    for raw_line in _read(constraint_path).splitlines():
        line = raw_line.split("#", maxsplit=1)[0].strip()
        if not line:
            continue
        requirement = line.split(";", maxsplit=1)[0].strip()
        if "==" not in requirement:
            continue
        name, version = requirement.split("==", maxsplit=1)
        versions[_requirement_name(name)] = version.strip()
    return versions


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
        - _pinned_constraint_names("constraints/ci.txt")
        - set(EXPLICITLY_JUSTIFIED_UNPINNED_DIRECT_INSTALLS)
    )

    assert unpinned == []


def test_minimum_constraints_pin_declared_runtime_and_test_lower_bounds() -> None:
    assert (
        _pinned_constraint_versions("constraints/minimum.txt")
        == MINIMUM_CONSTRAINT_LOWER_BOUND_PINS
    )


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
    assert "TWINE ?= $(PYTHON) -m twine" in makefile
    assert "$(TWINE) check dist/*" in makefile

    ci_build = _workflow_job_block(
        _read(".github/workflows/ci.yml"), "build-distributions"
    )
    publish_build = _workflow_job_block(_read(".github/workflows/publish.yml"), "build")
    for build_job in (ci_build, publish_build):
        assert (
            "python -m pip install -c constraints/ci.txt build twine setuptools wheel"
            in build_job
        )
    assert "make build" in ci_build
    assert "make release-check" in publish_build


def test_python_support_declaration_matches_source_test_matrices() -> None:
    pyproject = _load_pyproject()
    classifiers = set(pyproject["project"]["classifiers"])
    assert pyproject["project"]["requires-python"] == ">=3.10,<3.13"
    assert {
        classifier.rsplit("::", maxsplit=1)[1].strip()
        for classifier in classifiers
        if classifier.startswith("Programming Language :: Python :: 3.")
    } == set(SUPPORTED_PYTHON_VERSIONS)

    workflow = _read(".github/workflows/ci.yml")
    matrix_literal = _expected_matrix_literal()
    for job_name in (
        "clean-constrained-install",
        "default-suite",
        "performance-contracts",
        "activity-parity-gate",
        "parity-tests",
        "release-gates",
    ):
        assert matrix_literal in _workflow_job_block(workflow, job_name)


def test_minimum_dependency_ci_lane_uses_minimum_constraints_and_release_science_selectors() -> (
    None
):
    job = _workflow_job_block(
        _read(".github/workflows/ci.yml"), "minimum-dependency-suite"
    )

    assert "PIP_CONSTRAINT: ${{ github.workspace }}/constraints/minimum.txt" in job
    assert (
        "PIP_BUILD_CONSTRAINT: ${{ github.workspace }}/constraints/minimum.txt" in job
    )
    assert "constraints/ci.txt" not in job
    assert "python-version: '3.10'" in job
    assert 'python -m pip install -c constraints/minimum.txt -e ".[test]"' in job
    assert "python -m pip check" in job
    assert 'pytest --durations=25 --durations-min=0.01 -m "not parity"' in job
    assert (
        "pytest -o addopts= tests/release tests/golden "
        '-m "release_gate or golden or reproducibility"'
    ) in job


def test_ci_build_smoke_installs_one_wheel_outside_checkout() -> None:
    build_job = _workflow_job_block(
        _read(".github/workflows/ci.yml"), "build-distributions"
    )

    assert "make build" in build_job
    assert "dist/*.whl" in build_job
    assert "dist/*.tar.gz" in build_job
    assert "Expected exactly one wheel" in build_job
    assert "Expected exactly one sdist" in build_job
    assert "mktemp -d" in build_job
    assert "python -m venv" in build_job
    assert "-m pip check" in build_job
    assert 'cd "$smoke_dir"' in build_job
    assert 'python" -I -c' in build_job
    assert "import pathlib, phospy" in build_job


def test_publish_jobs_wait_for_single_release_check_build() -> None:
    workflow = _read(".github/workflows/publish.yml")
    build = _workflow_job_block(workflow, "build")

    assert "run: make release-check" in build
    assert "path: dist/" in build
    for publish_job in ("publish-to-testpypi", "publish-to-pypi"):
        publish_block = _workflow_job_block(workflow, publish_job)
        assert _has_needs(publish_block, "build")
        assert "uses: actions/download-artifact@v6" in publish_block
        assert "name: python-package-distributions" in publish_block
        assert "packages-dir: dist/" in publish_block
        assert "id-token: write" in publish_block
