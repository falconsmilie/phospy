from __future__ import annotations

import ast
import copy
import io
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import tomllib
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
FIXED_DISTRIBUTION_ARTIFACTS_ENV = "PHOSPY_TEST_DISTRIBUTION_ARTIFACTS"
CURRENT_RELEASE_VERSION = "1.7.6"
EXPECTED_FIXED_ARTIFACTS = {
    "wheel": (
        "phospy-1.7.6-py3-none-any.whl",
        1_731_940,
        "be16e6dba5945d1f6a8d96bef6060eebf983ebf76bb16d562f963b496f9271fa",
    ),
    "sdist": (
        "phospy-1.7.6.tar.gz",
        1_397_979,
        "f734950dfc5eb23b7751f3ab873d583920ca6efa510fd6ebb531b0c4458b803d",
    ),
}
WHEEL_DECLARED_RESOURCE = (
    "phospy/data/reference_bundles/rat/l6_native/substrate_map.csv"
)
WHEEL_ENRICHMENT_PROVENANCE_ASSEMBLER = "phospy/enrichment/provenance.py"
WHEEL_SPS_DISCOVERY_CONTRACTS = (
    "phospy/preprocessing/control_sites/discovery_contracts.py"
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
    fixed_dist = os.environ.get(FIXED_DISTRIBUTION_ARTIFACTS_ENV)
    if fixed_dist is not None:
        artifacts = find_distribution_artifacts(Path(fixed_dist))
        for artifact_kind, artifact in (
            ("wheel", artifacts.wheel),
            ("sdist", artifacts.sdist),
        ):
            expected_name, expected_size, expected_sha256 = EXPECTED_FIXED_ARTIFACTS[
                artifact_kind
            ]
            assert artifact.name == expected_name
            assert artifact.stat().st_size == expected_size
            assert _sha256_file(artifact) == expected_sha256
        return artifacts

    fixture_root = tmp_path_factory.mktemp("phospy-built-dist")
    dist = fixture_root / "dist"
    build_base = fixture_root / "build"
    dist.mkdir()
    build_script = """
import pathlib
import sys

import setuptools.build_meta as build_meta

dist = pathlib.Path(sys.argv[1])
build_base = pathlib.Path(sys.argv[2])
build_meta.build_sdist(str(dist))
build_meta.build_wheel(
    str(dist),
    {"--global-option": ["build", "--build-base", str(build_base)]},
)
"""
    command = [
        sys.executable,
        "-c",
        build_script,
        str(dist),
        str(build_base),
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


def _wheel_omitting_workflow_identity(source: Path, destination: Path) -> Path:
    original = b'workflow_name="enrichment",'
    replacement = b"workflow_name=None,"
    changed = False
    with zipfile.ZipFile(source, "r") as source_zip:
        with zipfile.ZipFile(destination, "w") as destination_zip:
            for item in source_zip.infolist():
                data = source_zip.read(item.filename)
                if item.filename == WHEEL_ENRICHMENT_PROVENANCE_ASSEMBLER:
                    if data.count(original) != 1:
                        raise AssertionError(
                            "wheel Enrichment provenance assembler did not contain "
                            "the expected workflow_name serialization"
                        )
                    data = data.replace(original, replacement)
                    changed = True
                destination_zip.writestr(item, data)
    if not changed:
        raise AssertionError(
            "wheel did not contain the Enrichment provenance assembler: "
            f"{WHEEL_ENRICHMENT_PROVENANCE_ASSEMBLER}"
        )
    return destination.resolve()


def _wheel_mutating_sps_provenance(source: Path, destination: Path) -> Path:
    original = b'SPS_DISCOVERY_ALGORITHM_ID = "phospy_sps_consensus_stability"'
    replacement = b'SPS_DISCOVERY_ALGORITHM_ID = "phospy_sps_consensus_changed__"'
    if len(original) != len(replacement):
        raise AssertionError("SPS provenance mutation must preserve archive size")
    changed = False
    with zipfile.ZipFile(source, "r") as source_zip:
        with zipfile.ZipFile(destination, "w") as destination_zip:
            for item in source_zip.infolist():
                data = source_zip.read(item.filename)
                if item.filename == WHEEL_SPS_DISCOVERY_CONTRACTS:
                    if data.count(original) != 1:
                        raise AssertionError(
                            "wheel SPS owner did not contain the expected algorithm ID"
                        )
                    data = data.replace(original, replacement)
                    changed = True
                destination_zip.writestr(item, data)
    if not changed:
        raise AssertionError(
            "wheel did not contain the SPS discovery contracts: "
            f"{WHEEL_SPS_DISCOVERY_CONTRACTS}"
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
    del tmp_path
    # Keep the isolated environment close to the runner-provided temporary root.
    # Deep pytest node directories otherwise exceed the legacy Windows path limit
    # when dependencies such as joblib install long fixture filenames.
    with tempfile.TemporaryDirectory(prefix=f"phospy-{artifact_kind}-damage-") as root:
        return _verify_one_artifact(
            artifact_kind=artifact_kind,
            artifact_path=artifact_path,
            work_root=Path(root),
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
            distribution_version=CURRENT_RELEASE_VERSION,
            runtime_version=CURRENT_RELEASE_VERSION,
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

    sdist_name = f"phospy-{CURRENT_RELEASE_VERSION}.tar.gz"

    verifier._run_pip(
        tmp_path / "venv" / "bin" / "python",
        "install sdist",
        repo_root=repo_root,
        cwd=cwd,
        constraint=constraint,
        arguments=["install", sdist_name],
    )

    assert captured["command"] == [
        str(tmp_path / "venv" / "bin" / "python"),
        "-m",
        "pip",
        "--disable-pip-version-check",
        "--no-cache-dir",
        "install",
        "-c",
        str(constraint),
        sdist_name,
    ]
    assert captured["environment_overrides"] == {
        "PIP_BUILD_CONSTRAINT": str(constraint),
        "PIP_NO_CACHE_DIR": "1",
    }


def test_project_version_authority_reports_current_release() -> None:
    with (ROOT / "pyproject.toml").open("rb") as handle:
        pyproject = tomllib.load(handle)

    assert pyproject["project"]["version"] == CURRENT_RELEASE_VERSION


def test_built_artifact_names_and_metadata_report_current_release_version(
    built_distribution_artifacts: DistributionArtifacts,
) -> None:
    assert (
        built_distribution_artifacts.wheel.name
        == f"phospy-{CURRENT_RELEASE_VERSION}-py3-none-any.whl"
    )
    assert (
        built_distribution_artifacts.sdist.name
        == f"phospy-{CURRENT_RELEASE_VERSION}.tar.gz"
    )
    assert (
        verifier._artifact_version(
            artifact_kind="wheel",
            artifact_path=built_distribution_artifacts.wheel,
        )
        == CURRENT_RELEASE_VERSION
    )
    assert (
        verifier._artifact_version(
            artifact_kind="sdist",
            artifact_path=built_distribution_artifacts.sdist,
        )
        == CURRENT_RELEASE_VERSION
    )


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
        assert report.distribution_version == CURRENT_RELEASE_VERSION
        assert report.runtime_version == CURRENT_RELEASE_VERSION
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
def test_installed_probe_rejects_omitted_run_provenance_workflow_identity(
    tmp_path: Path,
    built_distribution_artifacts: DistributionArtifacts,
) -> None:
    wheel = _wheel_omitting_workflow_identity(
        built_distribution_artifacts.wheel,
        tmp_path / built_distribution_artifacts.wheel.name,
    )

    with pytest.raises(
        InstalledDistributionVerificationError,
        match="complete RunProvenance payload changed",
    ):
        _verify_single_artifact(
            artifact_kind="wheel",
            artifact_path=wheel,
            tmp_path=tmp_path,
        )


@requires_supported_release_python
def test_installed_probe_rejects_changed_sps_provenance(
    tmp_path: Path,
    built_distribution_artifacts: DistributionArtifacts,
) -> None:
    wheel = _wheel_mutating_sps_provenance(
        built_distribution_artifacts.wheel,
        tmp_path / built_distribution_artifacts.wheel.name,
    )

    with pytest.raises(
        InstalledDistributionVerificationError,
        match="complete SPS discovery/control/correction payload changed",
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


def test_installed_probe_rejects_source_checkout_import_origin(tmp_path: Path) -> None:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    probe_path = tmp_path / "installed_distribution_probe.py"
    probe_path.write_text(INSTALLED_PROBE_SOURCE, encoding="utf-8", newline="\n")
    result = subprocess.run(
        [
            sys.executable,
            str(probe_path),
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
    assert "from phospy.signalome.requests import SignalomeWorkflowRequest" in (
        INSTALLED_PROBE_SOURCE
    )
    assert 'SignalomeWorkflowRequest.__module__ != "phospy.signalome.requests"' in (
        INSTALLED_PROBE_SOURCE
    )
    for retired_route in verifier.RETIRED_HORIZONTAL_ROUTES:
        assert retired_route in INSTALLED_PROBE_SOURCE
    assert '"exception_type": "ModuleNotFoundError"' in INSTALLED_PROBE_SOURCE
    assert 'route_result["module_not_found_error_name"]' in INSTALLED_PROBE_SOURCE
    assert 'module_name.startswith(missing_name + ".")' in INSTALLED_PROBE_SOURCE
    assert "manifest-declared resource is missing" in INSTALLED_PROBE_SOURCE
    assert "api.ReferenceBundleBuilder is not OwnerReferenceBundleBuilder" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "api.ReferenceBundleBuildRequest" in INSTALLED_PROBE_SOURCE
    assert "_exercise_reference_bundle_builder" in INSTALLED_PROBE_SOURCE
    assert 'bundle.site_sequences.to_dict(orient="index")' in INSTALLED_PROBE_SOURCE
    assert "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA" in INSTALLED_PROBE_SOURCE
    assert "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA" in INSTALLED_PROBE_SOURCE
    assert "expected_manifest_metadata" in INSTALLED_PROBE_SOURCE
    assert "expected_provenance_metadata" in INSTALLED_PROBE_SOURCE
    assert "expected_source_file_metadata" in INSTALLED_PROBE_SOURCE
    assert "expected_table_fingerprints" in INSTALLED_PROBE_SOURCE
    assert '"rows": 2' in INSTALLED_PROBE_SOURCE
    assert '"columns": 4' in INSTALLED_PROBE_SOURCE
    assert "sha256-stable-json-v1" in INSTALLED_PROBE_SOURCE
    assert (
        "a490fb6499d5c9deb87a7ce28284db5e916be44527efc9657827274520e696c6"
        in INSTALLED_PROBE_SOURCE
    )
    assert (
        "d6ca2a89d15e7a9882e23e7e2ff55325b84e818125b7ae903a218f5fd88e021a"
        in INSTALLED_PROBE_SOURCE
    )
    assert "_verify_references_evidence_import_orders" in INSTALLED_PROBE_SOURCE
    assert '("references", "evidence")' in INSTALLED_PROBE_SOURCE
    assert '("evidence", "references")' in INSTALLED_PROBE_SOURCE
    assert "aggregate_inventory_sha256" in INSTALLED_PROBE_SOURCE
    assert '"environment_versions"' in INSTALLED_PROBE_SOURCE
    assert '"runtime_platform"' in INSTALLED_PROBE_SOURCE
    assert '"dependency_identities"' in INSTALLED_PROBE_SOURCE
    assert "_exercise_enrichment_capability" in INSTALLED_PROBE_SOURCE
    assert "_verify_enrichment_import_orders" in INSTALLED_PROBE_SOURCE
    assert "_verify_control_site_import_orders" in INSTALLED_PROBE_SOURCE
    assert "_verify_batch_correction_boundary_import_orders" in INSTALLED_PROBE_SOURCE
    assert "_verify_protein_aware_boundary_import_orders" in INSTALLED_PROBE_SOURCE
    assert "ProteinAwarePreparationStageOutput" in INSTALLED_PROBE_SOURCE
    assert "compose_protein_aware_preparation_result" in INSTALLED_PROBE_SOURCE
    assert "phospy.datasets.models" in INSTALLED_PROBE_SOURCE
    assert "phospy.datasets.internal_view" in INSTALLED_PROBE_SOURCE
    assert "phospy.datasets.protein_aware" in INSTALLED_PROBE_SOURCE
    assert "phospy.science.datasets.builders.protein_aware_binding" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "phospy.science.datasets.preprocessing.protein_aware_models" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "phospy.contracts.results.preprocessing" in INSTALLED_PROBE_SOURCE
    assert "name.startswith('phospy.run_lineage.')" in INSTALLED_PROBE_SOURCE
    assert "name.startswith('phospy.datasets.')" in INSTALLED_PROBE_SOURCE
    assert "name.startswith('phospy.science.datasets.')" in INSTALLED_PROBE_SOURCE
    assert "phospy.provenance.serialization.batch_correction" in INSTALLED_PROBE_SOURCE
    assert (
        "phospy.preprocessing.batch_correction.adapters.preprocessing"
        in INSTALLED_PROBE_SOURCE
    )
    assert "synthetic_gene_sets" in INSTALLED_PROBE_SOURCE
    assert "synthetic_ptm_sets" in INSTALLED_PROBE_SOURCE
    assert "correction_owner" in INSTALLED_PROBE_SOURCE
    assert "run_provenance_to_payload" in INSTALLED_PROBE_SOURCE
    assert "_fixed_enrichment_environment" in INSTALLED_PROBE_SOURCE
    assert "1a20b51c3ae39867e2d005c12dafde31077dc89f88e290b5a4bba40b23e0727d" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "d34915364a6b40ef0d80a59bddd1acaeace302c314f15a33246dbc4a7a2cee56" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "ccc7cd07bd3a37df872906a6664cf8e820a5e24609c875baa6adaba59445b6d9" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "1e68f999766e008ce26433937343a064f0a88e9ce5f49f14958cf53a6540788f" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "manifest-declared resource hash mismatch" in INSTALLED_PROBE_SOURCE
    assert "load_bundled_reference_manifest" in INSTALLED_PROBE_SOURCE
    assert "DifferentialAnalysisWorkflow" in INSTALLED_PROBE_SOURCE
    assert "_exercise_duplicate_correlation_differential_contract" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION" in INSTALLED_PROBE_SOURCE
    assert "DatasetProteinAwarePreparationConfig" in INSTALLED_PROBE_SOURCE
    assert "DifferentialProteinAwareModelConfig" in INSTALLED_PROBE_SOURCE
    assert "protein_aware_configuration" in INSTALLED_PROBE_SOURCE
    assert (
        "DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1"
        in INSTALLED_PROBE_SOURCE
    )
    assert "KinaseWorkflow" in INSTALLED_PROBE_SOURCE
    assert "phospy.kinase.result_boundary" in INSTALLED_PROBE_SOURCE
    assert "installed Kinase boundary type ownership changed" in INSTALLED_PROBE_SOURCE
    assert (
        "installed Kinase selection policy identity changed" in INSTALLED_PROBE_SOURCE
    )
    assert "installed Kinase boundary retained a former route" in INSTALLED_PROBE_SOURCE
    assert "publish_dataset" in INSTALLED_PROBE_SOURCE
    assert "advanced_table_publisher" in INSTALLED_PROBE_SOURCE
    assert "save_kinase_workflow_bundle" in INSTALLED_PROBE_SOURCE
    assert "load_kinase_workflow_bundle" in INSTALLED_PROBE_SOURCE
    assert "publish_kinase_workflow" in INSTALLED_PROBE_SOURCE
    assert "_exercise_signalome_round_trip" in INSTALLED_PROBE_SOURCE
    assert "save_signalome_workflow_bundle" in INSTALLED_PROBE_SOURCE
    assert "load_signalome_workflow_bundle" in INSTALLED_PROBE_SOURCE
    assert "publish_signalome_workflow" in INSTALLED_PROBE_SOURCE
    assert "processing_state_round_trip" in INSTALLED_PROBE_SOURCE
    assert "_assert_frame_snapshot_isolation" in INSTALLED_PROBE_SOURCE
    assert "_assert_published_table_isolation" in INSTALLED_PROBE_SOURCE
    assert "ticket_1_posthoc_peptide_to_site_boundary" in INSTALLED_PROBE_SOURCE
    assert "withdrawn_asserted" in INSTALLED_PROBE_SOURCE
    assert "SpsDiscoveryWorkflow" in INSTALLED_PROBE_SOURCE
    assert "SpsDiscoveryRequest" in INSTALLED_PROBE_SOURCE
    assert "SpsReferenceDataset.from_condition_relative_log2" in INSTALLED_PROBE_SOURCE
    assert "RuvIIIReplicateStructure.from_assignments" in INSTALLED_PROBE_SOURCE
    assert "run_ruv_iii" in INSTALLED_PROBE_SOURCE
    assert "installed RUV-III changed the phospho matrix shape" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "installed RUV-III reported an unexpected method" in INSTALLED_PROBE_SOURCE
    assert '"ruv_iii_method"' in INSTALLED_PROBE_SOURCE
    assert '"ruv_iii_output_shape"' in INSTALLED_PROBE_SOURCE
    assert "BatchCorrectionWorkflow" in INSTALLED_PROBE_SOURCE
    assert "BatchCorrectionWorkflowRequest" in INSTALLED_PROBE_SOURCE
    assert "LinearResidualizeBatchCorrectionEngine" in INSTALLED_PROBE_SOURCE
    assert "default dataset build lost no-op correction execution" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "collect_batch_correction_environment_provenance" in (INSTALLED_PROBE_SOURCE)
    assert "installed RUV-III did not restore the governed cell" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "SpsRuvBatchCorrectionConfig" in INSTALLED_PROBE_SOURCE
    assert "batch_correction_provenance_to_payload" in INSTALLED_PROBE_SOURCE
    assert "batch_correction_provenance_from_payload" in INSTALLED_PROBE_SOURCE
    assert "run_provenance_from_payload" in INSTALLED_PROBE_SOURCE
    assert "installed direct correction provenance round trip changed" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "installed nested correction provenance round trip changed" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "installed malformed correction payload was accepted" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "installed correction recorded the wrong Python version" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "installed correction environment recorded the wrong Python version" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "installed correction recorded the wrong dependency versions" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "compared-with-rtol-1e-10-atol-1e-12" in INSTALLED_PROBE_SOURCE
    assert (
        "published processing-state payload unexpectedly contains batch_correction"
        in (INSTALLED_PROBE_SOURCE)
    )
    assert "_exercise_sps_ruv_contracts" in INSTALLED_PROBE_SOURCE
    assert '"sps_ruv": sps_ruv_report' in INSTALLED_PROBE_SOURCE
    assert "SpsDiscoveryResult.from_payload" in INSTALLED_PROBE_SOURCE
    assert "ControlSiteSet.from_payload" in INSTALLED_PROBE_SOURCE
    assert "corrected_preprocessing_output.output_observation_mask" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "corrected_preprocessing_output.corrected_cell_status" in (
        INSTALLED_PROBE_SOURCE
    )
    assert '"complete_correction"' in INSTALLED_PROBE_SOURCE
    assert '"sps_discovery_identity"' in INSTALLED_PROBE_SOURCE
    assert "inspect.signature(SpsDiscoveryRequest)" in INSTALLED_PROBE_SOURCE
    assert "inspect.signature(SpsDiscoveryWorkflow)" in INSTALLED_PROBE_SOURCE
    assert "OwnerSpsDiscoveryRequest" in INSTALLED_PROBE_SOURCE
    assert "OwnerSpsDiscoveryResult" in INSTALLED_PROBE_SOURCE
    assert "phospy-sps-discovery-result-v1" in INSTALLED_PROBE_SOURCE
    assert "legacy_restored.control_site_set.to_payload()" in INSTALLED_PROBE_SOURCE
    assert "legacy_controls = ControlSiteSet.from_payload(" in INSTALLED_PROBE_SOURCE
    assert 'legacy_controls.source_metadata.source_version != "1.0.0"' in (
        INSTALLED_PROBE_SOURCE
    )
    assert "legacy_controls.source_metadata.sps_discovery_identity" in (
        INSTALLED_PROBE_SOURCE
    )
    assert '"legacy_restored_controls"' in INSTALLED_PROBE_SOURCE
    assert "6f43b82b9ebaec9779fe1f9cf4fb02751369c715e797212d630a278a5e15fe23" in (
        INSTALLED_PROBE_SOURCE
    )
    assert "complete SPS discovery/control/correction payload changed" in (
        INSTALLED_PROBE_SOURCE
    )
