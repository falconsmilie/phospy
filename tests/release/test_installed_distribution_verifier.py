from __future__ import annotations

import ast
import copy
import io
import os
import shutil
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

import pytest

import scripts.verify_installed_distributions as verifier
from scripts.verify_installed_distributions import (
    INSTALLED_PROBE_SOURCE,
    DistributionArtifacts,
    InstalledDistributionReport,
    InstalledDistributionVerificationError,
    _probe_payload,
    _require_path_inside,
    _require_path_outside,
    _sha256_file,
    _verify_one_artifact,
    find_distribution_artifacts,
)

pytestmark = pytest.mark.release_gate

ROOT = Path(__file__).resolve().parents[2]
CONSTRAINT = ROOT / "constraints" / "ci.txt"
WHEEL_DECLARED_RESOURCE = (
    "phospy/data/reference_bundles/rat/l6_native/substrate_map.csv"
)
SDIST_DECLARED_RESOURCE_SUFFIX = (
    "/src/phospy/data/reference_bundles/rat/l6_native/substrate_map.csv"
)
DAMAGED_RESOURCE_SUFFIX = b"\n# deliberately damaged by release verifier test\n"
ARTIFACT_REGRESSION_INSTALL_KWARGS = {}
DAMAGED_ARTIFACT_INSTALL_KWARGS = {
    "install_packaging_tools": False,
}
SUPPORTED_RELEASE_PYTHON_VERSIONS = frozenset(
    tuple(int(part) for part in version.split("."))
    for version in verifier.SUPPORTED_PYTHON_VERSIONS
)
RUNNING_SUPPORTED_RELEASE_PYTHON = (
    sys.version_info.major,
    sys.version_info.minor,
) in SUPPORTED_RELEASE_PYTHON_VERSIONS
requires_supported_release_python = pytest.mark.skipif(
    not RUNNING_SUPPORTED_RELEASE_PYTHON,
    reason=(
        "installed distribution verification requires a supported release Python "
        f"({', '.join(verifier.SUPPORTED_PYTHON_VERSIONS)}); running "
        f"{sys.version_info.major}.{sys.version_info.minor}"
    ),
)


@pytest.fixture(scope="session")
def built_distribution_artifacts(
    tmp_path_factory: pytest.TempPathFactory,
) -> DistributionArtifacts:
    dist = tmp_path_factory.mktemp("phospy-built-dist") / "dist"
    dist.mkdir()
    build_script = """
import pathlib
import sys

import setuptools.build_meta as build_meta

dist = pathlib.Path(sys.argv[1])
build_meta.build_sdist(str(dist))
build_meta.build_wheel(str(dist))
"""
    command = [
        sys.executable,
        "-c",
        build_script,
        str(dist),
    ]
    result = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        pytest.fail(
            "failed to build clean wheel and sdist for installed verifier tests\n"
            f"command: {command!r}\n\nstdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"
        )
    return find_distribution_artifacts(dist)


def _copy_clean_dist(
    artifacts: DistributionArtifacts,
    target_dist: Path,
) -> DistributionArtifacts:
    target_dist.mkdir()
    wheel = target_dist / artifacts.wheel.name
    sdist = target_dist / artifacts.sdist.name
    shutil.copy2(artifacts.wheel, wheel)
    shutil.copy2(artifacts.sdist, sdist)
    return DistributionArtifacts(wheel=wheel.resolve(), sdist=sdist.resolve())


def _damaged_wheel(source: Path, destination: Path) -> Path:
    damaged = False
    with zipfile.ZipFile(source, "r") as source_zip:
        with zipfile.ZipFile(destination, "w") as destination_zip:
            for item in source_zip.infolist():
                data = source_zip.read(item.filename)
                if item.filename == WHEEL_DECLARED_RESOURCE:
                    data = data + DAMAGED_RESOURCE_SUFFIX
                    damaged = True
                destination_zip.writestr(item, data)
    if not damaged:
        raise AssertionError(
            f"wheel did not contain expected manifest resource: {WHEEL_DECLARED_RESOURCE}"
        )
    return destination.resolve()


def _damaged_sdist(source: Path, destination: Path) -> Path:
    damaged = False
    with tarfile.open(source, "r:gz") as source_tar:
        with tarfile.open(destination, "w:gz") as destination_tar:
            for member in source_tar.getmembers():
                member_copy = copy.copy(member)
                if member.isfile():
                    extracted = source_tar.extractfile(member)
                    if extracted is None:
                        raise AssertionError(
                            f"sdist member could not be read: {member.name}"
                        )
                    data = extracted.read()
                    if member.name.endswith(SDIST_DECLARED_RESOURCE_SUFFIX):
                        data = data + DAMAGED_RESOURCE_SUFFIX
                        member_copy.size = len(data)
                        damaged = True
                    destination_tar.addfile(member_copy, io.BytesIO(data))
                else:
                    destination_tar.addfile(member_copy)
    if not damaged:
        raise AssertionError(
            "sdist did not contain expected manifest resource suffix: "
            f"{SDIST_DECLARED_RESOURCE_SUFFIX}"
        )
    return destination.resolve()


def _verify_single_artifact(
    *,
    artifact_kind: str,
    artifact_path: Path,
    tmp_path: Path,
) -> InstalledDistributionReport:
    work_root = tmp_path / f"{artifact_kind}-verification-work"
    work_root.mkdir()
    return _verify_one_artifact(
        artifact_kind=artifact_kind,
        artifact_path=artifact_path,
        work_root=work_root,
        repo_root=ROOT.resolve(),
        python_executable=sys.executable,
        constraint=CONSTRAINT,
        **DAMAGED_ARTIFACT_INSTALL_KWARGS,
    )


def test_distribution_artifact_discovery_requires_exactly_one_wheel_and_sdist(
    tmp_path: Path,
) -> None:
    dist = tmp_path / "dist"
    dist.mkdir()
    wheel = dist / "phospy-1.6.0-py3-none-any.whl"
    sdist = dist / "phospy-1.6.0.tar.gz"
    wheel.write_bytes(b"placeholder wheel")
    sdist.write_bytes(b"placeholder sdist")

    artifacts = find_distribution_artifacts(dist)

    assert artifacts.wheel == wheel.resolve()
    assert artifacts.sdist == sdist.resolve()


def test_distribution_artifact_discovery_rejects_duplicates(tmp_path: Path) -> None:
    dist = tmp_path / "dist"
    dist.mkdir()
    (dist / "phospy-1.6.0-py3-none-any.whl").write_bytes(b"one")
    (dist / "phospy-1.6.0-1-py3-none-any.whl").write_bytes(b"two")
    (dist / "phospy-1.6.0.tar.gz").write_bytes(b"sdist")

    with pytest.raises(
        InstalledDistributionVerificationError,
        match="expected exactly one wheel",
    ):
        find_distribution_artifacts(dist)


def test_verifier_dispatches_wheel_and_sdist_to_isolated_artifact_checks(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dist = tmp_path / "dist"
    repo_root = tmp_path / "repo"
    constraint = tmp_path / "constraints.txt"
    dist.mkdir()
    repo_root.mkdir()
    constraint.write_text("", encoding="utf-8")
    wheel = dist / "phospy-1.6.0-py3-none-any.whl"
    sdist = dist / "phospy-1.6.0.tar.gz"
    wheel.write_bytes(b"placeholder wheel")
    sdist.write_bytes(b"placeholder sdist")
    calls: list[tuple[str, Path, Path, Path, str, Path | None]] = []

    def fake_verify_one_artifact(
        *,
        artifact_kind: str,
        artifact_path: Path,
        work_root: Path,
        repo_root: Path,
        python_executable: str,
        constraint: Path | None,
        use_system_site_packages: bool = False,
        install_dependencies: bool = True,
        build_isolation: bool = True,
        install_packaging_tools: bool = True,
        ignore_requires_python: bool = False,
    ) -> InstalledDistributionReport:
        del (
            use_system_site_packages,
            install_dependencies,
            build_isolation,
            install_packaging_tools,
            ignore_requires_python,
        )
        _require_path_outside(work_root, repo_root, label="verification tempdir")
        calls.append(
            (
                artifact_kind,
                artifact_path,
                work_root,
                repo_root,
                python_executable,
                constraint,
            )
        )
        environment_root = work_root / f"{artifact_kind}-venv"
        return InstalledDistributionReport(
            artifact_kind=artifact_kind,
            artifact_path=artifact_path,
            artifact_sha256="0" * 64,
            environment_root=environment_root,
            run_directory=work_root / f"{artifact_kind}-run",
            phospy_file=(
                environment_root / "Lib" / "site-packages" / "phospy" / "__init__.py"
            ),
            python_version="3.12.10",
            requires_python=verifier.EXPECTED_REQUIRES_PYTHON,
            resource_count=5,
            ticket_1_boundary_status="withdrawn_asserted",
            constraint_path=None if constraint is None else constraint.resolve(),
            constraint_sha256=None if constraint is None else _sha256_file(constraint),
        )

    monkeypatch.setattr(verifier, "_verify_one_artifact", fake_verify_one_artifact)

    reports = verifier.verify_installed_distributions(
        dist_dir=dist,
        repo_root=repo_root,
        python_executable="python-under-test",
        constraint=constraint,
    )

    assert [call[0] for call in calls] == ["wheel", "sdist"]
    assert [call[1] for call in calls] == [wheel.resolve(), sdist.resolve()]
    assert {call[2] for call in calls}
    assert all(call[3] == repo_root.resolve() for call in calls)
    assert all(call[4] == "python-under-test" for call in calls)
    assert all(call[5] == constraint.resolve() for call in calls)
    assert [report.artifact_kind for report in reports] == ["wheel", "sdist"]


def test_installed_artifact_regression_tests_keep_requires_python_enforced() -> None:
    assert "ignore_requires_python" not in ARTIFACT_REGRESSION_INSTALL_KWARGS
    assert "ignore_requires_python" not in DAMAGED_ARTIFACT_INSTALL_KWARGS


def test_verifier_subprocesses_ignore_ambient_pip_constraints(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repo_root = tmp_path / "repo"
    cwd = tmp_path / "work"
    repo_root.mkdir()
    cwd.mkdir()
    minimum_constraint = ROOT / "constraints" / "minimum.txt"
    captured_env: dict[str, str] = {}

    monkeypatch.setenv("PYTHONPATH", str(ROOT / "src"))
    monkeypatch.setenv("PYTHONHOME", "poison-python-home")
    monkeypatch.setenv("PIP_CONSTRAINT", str(minimum_constraint))
    monkeypatch.setenv("PIP_BUILD_CONSTRAINT", str(minimum_constraint))

    def fake_run(
        command: list[str],
        *,
        cwd: Path,
        env: dict[str, str],
        text: bool,
        capture_output: bool,
        check: bool,
    ) -> subprocess.CompletedProcess[str]:
        del cwd, text, capture_output, check
        captured_env.update(env)
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(verifier.subprocess, "run", fake_run)

    verifier._run(["probe"], cwd=cwd, repo_root=repo_root, context="probe")

    assert "PYTHONPATH" not in captured_env
    assert "PYTHONHOME" not in captured_env
    assert "PIP_CONSTRAINT" not in captured_env
    assert "PIP_BUILD_CONSTRAINT" not in captured_env
    assert captured_env["PIP_DISABLE_PIP_VERSION_CHECK"] == "1"


def test_pip_installs_apply_explicit_constraint_to_build_dependencies(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repo_root = tmp_path / "repo"
    cwd = tmp_path / "work"
    constraint = tmp_path / "constraints" / "ci.txt"
    repo_root.mkdir()
    cwd.mkdir()
    constraint.parent.mkdir()
    constraint.write_text("setuptools==77.0.3\n", encoding="utf-8")
    captured: dict[str, object] = {}

    def fake_run(
        command: list[str],
        *,
        cwd: Path,
        repo_root: Path,
        context: str,
        environment_overrides: dict[str, str] | None = None,
    ) -> subprocess.CompletedProcess[str]:
        captured.update(
            {
                "command": command,
                "cwd": cwd,
                "repo_root": repo_root,
                "context": context,
                "environment_overrides": environment_overrides,
            }
        )
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(verifier, "_run", fake_run)

    verifier._run_pip(
        tmp_path / "venv" / "bin" / "python",
        "install sdist",
        repo_root=repo_root,
        cwd=cwd,
        constraint=constraint,
        arguments=["install", "phospy-1.7.0.tar.gz"],
    )

    assert captured["command"] == [
        str(tmp_path / "venv" / "bin" / "python"),
        "-m",
        "pip",
        "--disable-pip-version-check",
        "install",
        "-c",
        str(constraint),
        "phospy-1.7.0.tar.gz",
    ]
    assert captured["environment_overrides"] == {
        "PIP_BUILD_CONSTRAINT": str(constraint)
    }


def test_built_artifact_metadata_declares_supported_python_contract(
    built_distribution_artifacts: DistributionArtifacts,
) -> None:
    assert (
        verifier._requires_python_specifiers(
            verifier._artifact_requires_python(
                artifact_kind="wheel",
                artifact_path=built_distribution_artifacts.wheel,
            )
        )
        == verifier.EXPECTED_REQUIRES_PYTHON_SPECIFIERS
    )
    assert (
        verifier._requires_python_specifiers(
            verifier._artifact_requires_python(
                artifact_kind="sdist",
                artifact_path=built_distribution_artifacts.sdist,
            )
        )
        == verifier.EXPECTED_REQUIRES_PYTHON_SPECIFIERS
    )


@requires_supported_release_python
def test_clean_wheel_and_sdist_install_and_execute_standalone_probe(
    tmp_path: Path,
    built_distribution_artifacts: DistributionArtifacts,
) -> None:
    dist = tmp_path / "clean-dist"
    _copy_clean_dist(built_distribution_artifacts, dist)

    reports = verifier.verify_installed_distributions(
        dist_dir=dist,
        repo_root=ROOT,
        python_executable=sys.executable,
        constraint=CONSTRAINT,
        **ARTIFACT_REGRESSION_INSTALL_KWARGS,
    )

    assert [report.artifact_kind for report in reports] == ["wheel", "sdist"]
    for report in reports:
        _require_path_inside(
            report.phospy_file,
            report.environment_root,
            label="phospy.__file__",
        )
        _require_path_outside(report.phospy_file, ROOT, label="phospy.__file__")
        assert (
            verifier._requires_python_specifiers(report.requires_python)
            == verifier.EXPECTED_REQUIRES_PYTHON_SPECIFIERS
        )
        assert report.resource_count >= 5
        assert report.ticket_1_boundary_status == "withdrawn_asserted"
        assert report.artifact_sha256 == _sha256_file(report.artifact_path)
        assert report.constraint_path == CONSTRAINT.resolve()
        assert report.constraint_sha256 == _sha256_file(CONSTRAINT)


@requires_supported_release_python
def test_damaged_wheel_manifest_resource_fails_installed_probe(
    tmp_path: Path,
    built_distribution_artifacts: DistributionArtifacts,
) -> None:
    wheel = _damaged_wheel(
        built_distribution_artifacts.wheel,
        tmp_path / built_distribution_artifacts.wheel.name,
    )

    with pytest.raises(
        InstalledDistributionVerificationError,
        match="manifest-declared resource hash mismatch",
    ):
        _verify_single_artifact(
            artifact_kind="wheel",
            artifact_path=wheel,
            tmp_path=tmp_path,
        )


@requires_supported_release_python
def test_damaged_sdist_manifest_resource_fails_installed_probe(
    tmp_path: Path,
    built_distribution_artifacts: DistributionArtifacts,
) -> None:
    sdist = _damaged_sdist(
        built_distribution_artifacts.sdist,
        tmp_path / built_distribution_artifacts.sdist.name,
    )

    with pytest.raises(
        InstalledDistributionVerificationError,
        match="manifest-declared resource hash mismatch",
    ):
        _verify_single_artifact(
            artifact_kind="sdist",
            artifact_path=sdist,
            tmp_path=tmp_path,
        )


def test_installed_probe_rejects_source_checkout_import_origin() -> None:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            INSTALLED_PROBE_SOURCE,
            str(ROOT.resolve()),
            str(ROOT.resolve()),
            "source-checkout",
        ],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode != 0
    assert "phospy.__file__ unexpectedly resolves inside checkout" in result.stderr


def test_probe_payload_rejects_failed_status_report() -> None:
    with pytest.raises(
        InstalledDistributionVerificationError,
        match="did not report ok status",
    ):
        _probe_payload(
            '{"status": "failed", "reason": "manifest-declared resource hash mismatch"}',
            artifact_kind="wheel",
        )


def test_probe_path_guards_enforce_installed_origin_outside_checkout(
    tmp_path: Path,
) -> None:
    repo_root = tmp_path / "repo"
    environment_root = tmp_path / "venv"
    installed_file = (
        environment_root / "Lib" / "site-packages" / "phospy" / "__init__.py"
    )
    checkout_file = repo_root / "src" / "phospy" / "__init__.py"
    installed_file.parent.mkdir(parents=True)
    checkout_file.parent.mkdir(parents=True)
    installed_file.write_text("", encoding="utf-8")
    checkout_file.write_text("", encoding="utf-8")

    _require_path_inside(installed_file, environment_root, label="phospy.__file__")
    _require_path_outside(installed_file, repo_root, label="phospy.__file__")
    with pytest.raises(
        InstalledDistributionVerificationError,
        match="unexpectedly resolves inside checkout",
    ):
        _require_path_outside(checkout_file, repo_root, label="phospy.__file__")


def test_installed_probe_source_avoids_repository_tests_and_fixtures() -> None:
    tree = ast.parse(INSTALLED_PROBE_SOURCE)
    imported_modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imported_modules.add(node.module)

    assert not any(
        module == "tests" or module.startswith("tests.") for module in imported_modules
    )
    assert "conftest" not in imported_modules
    assert "PYTHONPATH" not in INSTALLED_PROBE_SOURCE
    assert "manifest-declared resource is missing" in INSTALLED_PROBE_SOURCE
    assert "manifest-declared resource hash mismatch" in INSTALLED_PROBE_SOURCE
    assert "load_bundled_reference_manifest" in INSTALLED_PROBE_SOURCE
    assert "DifferentialAnalysisWorkflow" in INSTALLED_PROBE_SOURCE
    assert "_exercise_duplicate_correlation_differential_contract" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION" in INSTALLED_PROBE_SOURCE
    assert "KinaseWorkflow" in INSTALLED_PROBE_SOURCE
    assert "publish_dataset" in INSTALLED_PROBE_SOURCE
    assert "advanced_table_publisher" in INSTALLED_PROBE_SOURCE
    assert "save_kinase_workflow_bundle" in INSTALLED_PROBE_SOURCE
    assert "load_kinase_workflow_bundle" in INSTALLED_PROBE_SOURCE
    assert "ticket_1_posthoc_peptide_to_site_boundary" in INSTALLED_PROBE_SOURCE
    assert "withdrawn_asserted" in INSTALLED_PROBE_SOURCE
