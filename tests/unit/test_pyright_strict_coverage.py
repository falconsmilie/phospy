from __future__ import annotations

import json
import subprocess
import sys
from collections.abc import Iterable
from pathlib import Path

from tools.testing.pyright_strict_coverage import (
    DEFAULT_STRICT_COVERAGE_FLOOR_PERCENT,
    ConfiguredStrictPathIssue,
    StrictSuppressionPolicyIssue,
    pyright_strict_coverage,
)

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "tools" / "testing" / "pyright_strict_coverage.py"


def test_repository_strict_coverage_policy_is_clean_and_above_floor() -> None:
    coverage = pyright_strict_coverage(ROOT)

    assert coverage.included_source_file_count > 0
    assert coverage.declared_strict_source_file_count >= 95
    assert coverage.strict_source_files_inside_include_count == (
        coverage.declared_strict_source_file_count
    )
    assert coverage.strict_source_files_outside_include == ()
    assert coverage.missing_or_invalid_configured_strict_paths == ()
    assert coverage.strict_suppression_policy_issues == ()
    assert coverage.strict_source_percent >= DEFAULT_STRICT_COVERAGE_FLOOR_PERCENT
    assert coverage.passes_policy


def test_contract_results_contribute_to_strict_coverage() -> None:
    coverage = pyright_strict_coverage(ROOT)
    contract_result_files = {
        path.relative_to(ROOT).as_posix()
        for path in (ROOT / "src" / "phospy" / "contracts" / "results").glob("*.py")
    }

    assert contract_result_files
    assert contract_result_files <= set(coverage.declared_strict_source_files)
    assert contract_result_files <= set(coverage.strict_source_files_inside_include)


def test_strict_file_outside_include_is_reported_as_policy_failure(
    tmp_path: Path,
) -> None:
    _write_pyright_config(
        tmp_path,
        include_paths=("pkg/included",),
        strict_paths=("pkg/outside.py",),
    )
    _write_source(tmp_path, "pkg/included/module.py")
    _write_source(tmp_path, "pkg/outside.py")

    coverage = pyright_strict_coverage(tmp_path, strict_coverage_floor_percent=0.0)

    assert coverage.included_source_files == ("pkg/included/module.py",)
    assert coverage.declared_strict_source_files == ("pkg/outside.py",)
    assert coverage.strict_source_files_inside_include == ()
    assert coverage.strict_source_files_outside_include == ("pkg/outside.py",)
    assert any(
        "outside [tool.pyright].include" in issue
        for issue in (coverage.policy_failure_messages())
    )
    assert not coverage.passes_policy


def test_nonexistent_strict_path_is_reported_as_policy_failure(
    tmp_path: Path,
) -> None:
    _write_pyright_config(
        tmp_path,
        include_paths=("pkg",),
        strict_paths=("pkg/missing.py",),
    )
    _write_source(tmp_path, "pkg/module.py")

    coverage = pyright_strict_coverage(tmp_path, strict_coverage_floor_percent=0.0)

    assert coverage.declared_strict_source_files == ()
    assert coverage.missing_or_invalid_configured_strict_paths == (
        ConfiguredStrictPathIssue(
            path="pkg/missing.py",
            reason="does_not_exist",
        ),
    )
    assert not coverage.passes_policy


def test_non_python_strict_path_is_reported_as_invalid_configuration(
    tmp_path: Path,
) -> None:
    _write_pyright_config(
        tmp_path,
        include_paths=("pkg",),
        strict_paths=("pkg/notes.txt",),
    )
    _write_source(tmp_path, "pkg/module.py")
    notes_path = tmp_path / "pkg" / "notes.txt"
    notes_path.write_text("not Python source\n", encoding="utf-8")

    coverage = pyright_strict_coverage(tmp_path, strict_coverage_floor_percent=0.0)

    assert coverage.declared_strict_source_files == ()
    assert coverage.missing_or_invalid_configured_strict_paths == (
        ConfiguredStrictPathIssue(
            path="pkg/notes.txt",
            reason="not_python_source_file",
        ),
    )
    assert not coverage.passes_policy


def test_directory_and_individual_file_strict_entries_are_reported(
    tmp_path: Path,
) -> None:
    _write_pyright_config(
        tmp_path,
        include_paths=("pkg",),
        strict_paths=("pkg/strict_dir", "pkg/strict_file.py"),
    )
    _write_source(tmp_path, "pkg/strict_dir/a.py")
    _write_source(tmp_path, "pkg/strict_dir/b.py")
    _write_source(tmp_path, "pkg/strict_file.py")

    coverage = pyright_strict_coverage(tmp_path, strict_coverage_floor_percent=0.0)

    assert coverage.declared_strict_source_files == (
        "pkg/strict_dir/a.py",
        "pkg/strict_dir/b.py",
        "pkg/strict_file.py",
    )
    assert coverage.strict_source_files_inside_include == (
        "pkg/strict_dir/a.py",
        "pkg/strict_dir/b.py",
        "pkg/strict_file.py",
    )
    assert coverage.strict_source_files_outside_include == ()
    assert coverage.passes_policy


def test_reasoned_rule_specific_pyright_ignore_passes(tmp_path: Path) -> None:
    _write_pyright_config(
        tmp_path,
        include_paths=("pkg",),
        strict_paths=("pkg/strict.py",),
    )
    _write_source(
        tmp_path,
        "pkg/strict.py",
        (
            "VALUE = object()  # pyright: ignore[reportUnknownMemberType] - "
            "third-party descriptor returns a runtime int but its stub exposes object.\n"
        ),
    )

    coverage = pyright_strict_coverage(tmp_path, strict_coverage_floor_percent=0.0)

    assert coverage.strict_suppression_policy_issues == ()
    assert coverage.passes_policy


def test_rule_specific_pyright_ignore_without_reason_fails(tmp_path: Path) -> None:
    _write_pyright_config(
        tmp_path,
        include_paths=("pkg",),
        strict_paths=("pkg/strict.py",),
    )
    _write_source(
        tmp_path,
        "pkg/strict.py",
        "VALUE = object()  # pyright: ignore[reportUnknownMemberType]\n",
    )

    coverage = pyright_strict_coverage(tmp_path, strict_coverage_floor_percent=0.0)

    assert coverage.strict_suppression_policy_issues == (
        StrictSuppressionPolicyIssue(
            path="pkg/strict.py",
            line=1,
            suppression="reportUnknownMemberType",
            reason=(
                "rule-specific Pyright ignore must include ' - ' followed "
                "by a technical rationale"
            ),
        ),
    )
    assert not coverage.passes_policy


def test_blanket_inline_pyright_ignore_fails(tmp_path: Path) -> None:
    _write_pyright_config(
        tmp_path,
        include_paths=("pkg",),
        strict_paths=("pkg/strict.py",),
    )
    _write_source(
        tmp_path,
        "pkg/strict.py",
        "VALUE = object()  # pyright: ignore\n",
    )

    coverage = pyright_strict_coverage(tmp_path, strict_coverage_floor_percent=0.0)

    assert coverage.strict_suppression_policy_issues == (
        StrictSuppressionPolicyIssue(
            path="pkg/strict.py",
            line=1,
            suppression="pyright: ignore",
            reason=(
                "blanket Pyright ignores are not allowed in effective strict files"
            ),
        ),
    )
    assert not coverage.passes_policy


def test_file_wide_strictness_downgrade_fails(tmp_path: Path) -> None:
    _write_pyright_config(
        tmp_path,
        include_paths=("pkg",),
        strict_paths=("pkg/strict.py",),
    )
    _write_source(
        tmp_path,
        "pkg/strict.py",
        "# pyright: basic\nVALUE: int = 1\n",
    )

    coverage = pyright_strict_coverage(tmp_path, strict_coverage_floor_percent=0.0)

    assert coverage.strict_suppression_policy_issues == (
        StrictSuppressionPolicyIssue(
            path="pkg/strict.py",
            line=1,
            suppression="pyright: basic",
            reason=(
                "effective strict files may not downgrade file-wide Pyright "
                "type-checking mode"
            ),
        ),
    )
    assert not coverage.passes_policy


def test_file_wide_diagnostic_suppression_fails(tmp_path: Path) -> None:
    _write_pyright_config(
        tmp_path,
        include_paths=("pkg",),
        strict_paths=("pkg/strict.py",),
    )
    _write_source(
        tmp_path,
        "pkg/strict.py",
        "# pyright: reportUnknownMemberType=false\nVALUE: int = 1\n",
    )

    coverage = pyright_strict_coverage(tmp_path, strict_coverage_floor_percent=0.0)

    assert coverage.strict_suppression_policy_issues == (
        StrictSuppressionPolicyIssue(
            path="pkg/strict.py",
            line=1,
            suppression="reportUnknownMemberType=false",
            reason=(
                "effective strict files may not use file-wide diagnostic "
                "severity overrides; use a line-local rule-specific "
                "suppression with a rationale if unavoidable"
            ),
        ),
    )
    assert not coverage.passes_policy


def test_placeholder_pyright_ignore_reason_fails(tmp_path: Path) -> None:
    _write_pyright_config(
        tmp_path,
        include_paths=("pkg",),
        strict_paths=("pkg/strict.py",),
    )
    _write_source(
        tmp_path,
        "pkg/strict.py",
        "VALUE = object()  # pyright: ignore[reportUnknownMemberType] - TODO\n",
    )

    coverage = pyright_strict_coverage(tmp_path, strict_coverage_floor_percent=0.0)

    assert coverage.strict_suppression_policy_issues == (
        StrictSuppressionPolicyIssue(
            path="pkg/strict.py",
            line=1,
            suppression="reportUnknownMemberType",
            reason=(
                "rule-specific Pyright ignore rationale must describe the "
                "concrete typing limitation and runtime safety"
            ),
        ),
    )
    assert not coverage.passes_policy


def test_unexplained_suppression_outside_effective_strict_set_is_ignored(
    tmp_path: Path,
) -> None:
    _write_pyright_config(
        tmp_path,
        include_paths=("pkg",),
        strict_paths=("pkg/strict.py",),
    )
    _write_source(tmp_path, "pkg/strict.py")
    _write_source(
        tmp_path,
        "pkg/loose.py",
        "VALUE = object()  # pyright: ignore\n",
    )

    coverage = pyright_strict_coverage(tmp_path, strict_coverage_floor_percent=0.0)

    assert coverage.strict_suppression_policy_issues == ()
    assert coverage.passes_policy


def test_package_wide_ignore_intersecting_strict_path_fails(tmp_path: Path) -> None:
    _write_pyright_config(
        tmp_path,
        include_paths=("pkg",),
        strict_paths=("pkg/strict",),
        ignore_paths=("pkg/strict",),
    )
    _write_source(tmp_path, "pkg/strict/a.py")

    coverage = pyright_strict_coverage(tmp_path, strict_coverage_floor_percent=0.0)

    assert coverage.strict_suppression_policy_issues == (
        StrictSuppressionPolicyIssue(
            path="pyproject.toml",
            line=4,
            suppression="[tool.pyright].ignore=pkg/strict",
            reason=(
                "package-wide Pyright ignore intersects effective strict "
                "source file pkg/strict/a.py"
            ),
        ),
    )
    assert not coverage.passes_policy


def test_duplicate_and_overlapping_strict_entries_are_deduplicated(
    tmp_path: Path,
) -> None:
    _write_pyright_config(
        tmp_path,
        include_paths=("pkg",),
        strict_paths=("pkg/strict", "pkg/strict/a.py", "pkg/strict"),
    )
    _write_source(tmp_path, "pkg/strict/a.py")
    _write_source(tmp_path, "pkg/strict/b.py")

    coverage = pyright_strict_coverage(tmp_path, strict_coverage_floor_percent=0.0)

    assert coverage.declared_strict_source_files == (
        "pkg/strict/a.py",
        "pkg/strict/b.py",
    )
    assert coverage.strict_source_files_inside_include == (
        "pkg/strict/a.py",
        "pkg/strict/b.py",
    )
    assert coverage.missing_or_invalid_configured_strict_paths == ()


def test_check_mode_exits_nonzero_for_invalid_configuration(tmp_path: Path) -> None:
    _write_pyright_config(
        tmp_path,
        include_paths=("pkg",),
        strict_paths=("pkg/missing.py",),
    )
    _write_source(tmp_path, "pkg/module.py")

    result = _run_coverage_cli(
        tmp_path,
        "--check",
        "--floor-percent",
        "0",
    )

    assert result.returncode == 1
    assert "missing or invalid configured strict paths: 1" in result.stderr


def test_check_mode_exits_nonzero_when_floor_is_not_met(tmp_path: Path) -> None:
    _write_pyright_config(
        tmp_path,
        include_paths=("pkg",),
        strict_paths=("pkg/a.py",),
    )
    _write_source(tmp_path, "pkg/a.py")
    _write_source(tmp_path, "pkg/b.py")
    _write_source(tmp_path, "pkg/c.py")

    result = _run_coverage_cli(
        tmp_path,
        "--check",
        "--floor-percent",
        "50",
    )

    assert result.returncode == 1
    assert "strict coverage floor not met" in result.stderr


def test_json_output_is_stable_and_machine_readable(tmp_path: Path) -> None:
    _write_pyright_config(
        tmp_path,
        include_paths=("pkg",),
        strict_paths=("pkg/a.py",),
    )
    _write_source(tmp_path, "pkg/a.py")

    result = _run_coverage_cli(
        tmp_path,
        "--json",
        "--floor-percent",
        "0",
    )

    assert result.returncode == 0
    assert json.loads(result.stdout) == {
        "declared_strict_source_file_count": 1,
        "declared_strict_source_files": ["pkg/a.py"],
        "included_source_file_count": 1,
        "included_source_files": ["pkg/a.py"],
        "missing_or_invalid_configured_strict_paths": [],
        "policy_failures": [],
        "strict_suppression_policy_issues": [],
        "strict_coverage_floor_percent": 0.0,
        "strict_source_files_inside_include": ["pkg/a.py"],
        "strict_source_files_inside_include_count": 1,
        "strict_source_files_outside_include": [],
        "strict_source_files_outside_include_count": 0,
        "strict_source_fraction": 1.0,
        "strict_source_percent": 100.0,
    }


def _write_pyright_config(
    repo_root: Path,
    *,
    include_paths: Iterable[str],
    strict_paths: Iterable[str],
    ignore_paths: Iterable[str] = (),
) -> None:
    include_payload = _toml_string_list(include_paths)
    strict_payload = _toml_string_list(strict_paths)
    ignore_payload = _toml_string_list(ignore_paths)
    ignore_line = "" if not ignore_payload else f"ignore = [{ignore_payload}]\n"
    (repo_root / "pyproject.toml").write_text(
        "[tool.pyright]\n"
        f"include = [{include_payload}]\n"
        f"strict = [{strict_payload}]\n"
        f"{ignore_line}"
        'exclude = ["**/__pycache__", ".venv", "build", "dist"]\n',
        encoding="utf-8",
    )


def _toml_string_list(values: Iterable[str]) -> str:
    return ", ".join(json.dumps(value) for value in values)


def _write_source(
    repo_root: Path,
    relative_path: str,
    content: str = "VALUE: int = 1\n",
) -> None:
    path = repo_root / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _run_coverage_cli(repo_root: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--repo-root",
            str(repo_root),
            *args,
        ],
        check=False,
        capture_output=True,
        text=True,
    )
