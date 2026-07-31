from __future__ import annotations

import ast
from pathlib import Path

import pytest

import scripts.verify_installed_distributions as verifier
from scripts.verify_installed_distributions import (
    INSTALLED_PROBE_SOURCE,
    InstalledDistributionReport,
    InstalledDistributionVerificationError,
    _probe_payload,
    _require_path_inside,
    _require_path_outside,
    find_distribution_artifacts,
)

pytestmark = pytest.mark.release_gate


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
    ) -> InstalledDistributionReport:
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
            environment_root=environment_root,
            run_directory=work_root / f"{artifact_kind}-run",
            phospy_file=(
                environment_root / "Lib" / "site-packages" / "phospy" / "__init__.py"
            ),
            python_version="3.12.10",
            resource_count=5,
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


def test_probe_payload_rejects_damaged_artifact_report() -> None:
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
    assert "KinaseWorkflow" in INSTALLED_PROBE_SOURCE
