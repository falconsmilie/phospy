from __future__ import annotations

import argparse
import fnmatch
import json
import re
import sys
import tomllib
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

DEFAULT_STRICT_COVERAGE_FLOOR_PERCENT = 15.0

StrictPathIssueReason = Literal[
    "outside_repository",
    "does_not_exist",
    "not_python_source_file",
    "no_python_source_files",
]

_EXCLUDED_PATH_PARTS = frozenset(
    {
        "__pycache__",
        ".pytest_cache",
        ".venv",
        "build",
        "dist",
        "legacy_archive",
        "src/phospy.egg-info",
    }
)
_PYRIGHT_DIRECTIVE_RE = re.compile(r"#\s*pyright:\s*(?P<body>.*)$")
_PYRIGHT_IGNORE_RE = re.compile(
    r"^ignore(?P<rules>\[(?P<rule_list>[^\]]*)\])?(?P<trailing>.*)$"
)
_PYRIGHT_RULE_ID_RE = re.compile(r"^report[A-Za-z][A-Za-z0-9]*$")
_FILE_WIDE_DIAGNOSTIC_SUPPRESSION_RE = re.compile(
    r"\b(?P<rule>report[A-Za-z][A-Za-z0-9]*)\s*=\s*"
    r"(?P<value>false|none|warning|information|hint)\b",
    re.IGNORECASE,
)
_FILE_WIDE_TYPE_CHECKING_DOWNGRADE_RE = re.compile(
    r"^(?P<mode>basic|standard|off)\b"
    r"|(?:\btypeCheckingMode\s*=\s*(?P<assigned_mode>basic|standard|off)\b)"
    r"|(?:\bstrict\s*=\s*false\b)",
    re.IGNORECASE,
)
_TYPE_IGNORE_RE = re.compile(r"#\s*type:\s*ignore(?:\[[^\]]*\])?")
_PLACEHOLDER_RATIONALE_TOKENS = frozenset(
    {
        "fixme",
        "todo",
        "tbd",
        "xxx",
    }
)
_GENERIC_RATIONALES = frozenset(
    {
        "because pyright",
        "needed for pandas",
        "none",
        "n/a",
        "na",
        "pyright bug",
        "pyright is wrong",
        "typing",
        "typing issue",
        "type checker issue",
    }
)


@dataclass(frozen=True, slots=True)
class ConfiguredStrictPathIssue:
    path: str
    reason: StrictPathIssueReason

    def to_payload(self) -> dict[str, str]:
        return {
            "path": self.path,
            "reason": self.reason,
        }


@dataclass(frozen=True, slots=True)
class StrictSuppressionPolicyIssue:
    path: str
    line: int
    suppression: str
    reason: str

    def to_payload(self) -> dict[str, object]:
        return {
            "path": self.path,
            "line": self.line,
            "suppression": self.suppression,
            "reason": self.reason,
        }


@dataclass(frozen=True, slots=True)
class _ConfiguredPathSourceFiles:
    source_files: set[str]
    issues: tuple[ConfiguredStrictPathIssue, ...]


@dataclass(frozen=True, slots=True)
class PyrightStrictCoverage:
    included_source_files: tuple[str, ...]
    declared_strict_source_files: tuple[str, ...]
    strict_source_files_inside_include: tuple[str, ...]
    strict_source_files_outside_include: tuple[str, ...]
    missing_or_invalid_configured_strict_paths: tuple[ConfiguredStrictPathIssue, ...]
    strict_suppression_policy_issues: tuple[StrictSuppressionPolicyIssue, ...]
    strict_coverage_floor_percent: float = DEFAULT_STRICT_COVERAGE_FLOOR_PERCENT

    @property
    def included_source_file_count(self) -> int:
        return len(self.included_source_files)

    @property
    def declared_strict_source_file_count(self) -> int:
        return len(self.declared_strict_source_files)

    @property
    def strict_source_files_inside_include_count(self) -> int:
        return len(self.strict_source_files_inside_include)

    @property
    def strict_source_files_outside_include_count(self) -> int:
        return len(self.strict_source_files_outside_include)

    @property
    def strict_source_fraction(self) -> float:
        if self.included_source_file_count == 0:
            return 0.0
        return self.strict_source_files_inside_include_count / (
            self.included_source_file_count
        )

    @property
    def strict_source_percent(self) -> float:
        return self.strict_source_fraction * 100.0

    def policy_failure_messages(self) -> tuple[str, ...]:
        messages: list[str] = []
        if self.strict_source_files_outside_include_count:
            messages.append(
                "declared strict source files outside [tool.pyright].include: "
                f"{self.strict_source_files_outside_include_count}"
            )
        if self.missing_or_invalid_configured_strict_paths:
            messages.append(
                "missing or invalid configured strict paths: "
                f"{len(self.missing_or_invalid_configured_strict_paths)}"
            )
        if self.strict_suppression_policy_issues:
            messages.append(
                "strict suppression policy violations: "
                f"{len(self.strict_suppression_policy_issues)}"
            )
        if self.strict_source_percent < self.strict_coverage_floor_percent:
            messages.append(
                "strict coverage floor not met: "
                f"{self.strict_source_percent:.1f}% < "
                f"{self.strict_coverage_floor_percent:.1f}%"
            )
        return tuple(messages)

    @property
    def passes_policy(self) -> bool:
        return not self.policy_failure_messages()

    def to_payload(self) -> dict[str, object]:
        return {
            "declared_strict_source_file_count": (
                self.declared_strict_source_file_count
            ),
            "declared_strict_source_files": list(self.declared_strict_source_files),
            "included_source_file_count": self.included_source_file_count,
            "included_source_files": list(self.included_source_files),
            "missing_or_invalid_configured_strict_paths": [
                issue.to_payload()
                for issue in self.missing_or_invalid_configured_strict_paths
            ],
            "policy_failures": list(self.policy_failure_messages()),
            "strict_suppression_policy_issues": [
                issue.to_payload() for issue in self.strict_suppression_policy_issues
            ],
            "strict_coverage_floor_percent": self.strict_coverage_floor_percent,
            "strict_source_files_inside_include": list(
                self.strict_source_files_inside_include
            ),
            "strict_source_files_inside_include_count": (
                self.strict_source_files_inside_include_count
            ),
            "strict_source_files_outside_include": list(
                self.strict_source_files_outside_include
            ),
            "strict_source_files_outside_include_count": (
                self.strict_source_files_outside_include_count
            ),
            "strict_source_fraction": self.strict_source_fraction,
            "strict_source_percent": self.strict_source_percent,
        }


def pyright_strict_coverage(
    repo_root: Path,
    *,
    strict_coverage_floor_percent: float = DEFAULT_STRICT_COVERAGE_FLOOR_PERCENT,
) -> PyrightStrictCoverage:
    repo_root = repo_root.resolve()
    pyright_config = _pyright_config(repo_root)
    exclude_paths = _optional_string_list(pyright_config, "exclude")
    included = _source_files_for_paths(
        repo_root,
        _string_list(pyright_config, "include"),
        exclude_paths=exclude_paths,
    )
    strict = _declared_strict_source_files_for_paths(
        repo_root, _string_list(pyright_config, "strict")
    )
    strict_inside_include = strict.source_files & included
    strict_outside_include = strict.source_files - included
    suppression_issues = _strict_suppression_policy_issues(
        repo_root,
        pyright_config=pyright_config,
        strict_source_files=strict_inside_include,
    )
    return PyrightStrictCoverage(
        included_source_files=tuple(sorted(included)),
        declared_strict_source_files=tuple(sorted(strict.source_files)),
        strict_source_files_inside_include=tuple(sorted(strict_inside_include)),
        strict_source_files_outside_include=tuple(sorted(strict_outside_include)),
        missing_or_invalid_configured_strict_paths=strict.issues,
        strict_suppression_policy_issues=suppression_issues,
        strict_coverage_floor_percent=strict_coverage_floor_percent,
    )


def _pyright_config(repo_root: Path) -> Mapping[str, object]:
    pyproject = repo_root / "pyproject.toml"
    payload = tomllib.loads(pyproject.read_text(encoding="utf-8"))
    tool = payload.get("tool")
    if not isinstance(tool, Mapping):
        raise ValueError("pyproject.toml is missing [tool]")
    pyright = tool.get("pyright")
    if not isinstance(pyright, Mapping):
        raise ValueError("pyproject.toml is missing [tool.pyright]")
    return pyright


def _string_list(payload: Mapping[str, object], key: str) -> tuple[str, ...]:
    value = payload.get(key)
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ValueError(f"[tool.pyright].{key} must be a list of strings")
    return tuple(value)


def _optional_string_list(payload: Mapping[str, object], key: str) -> tuple[str, ...]:
    value = payload.get(key, [])
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ValueError(f"[tool.pyright].{key} must be a list of strings")
    return tuple(value)


def _source_files_for_paths(
    repo_root: Path,
    paths: Iterable[str],
    *,
    exclude_paths: Iterable[str],
) -> set[str]:
    files: set[str] = set()
    for path_text in paths:
        path = _resolve_configured_path_under_repo(repo_root, path_text)
        if path is None:
            continue
        if path.is_file():
            relative_path = _normalized_relative_path(repo_root, path)
            if path.suffix == ".py" and not _is_excluded_source_path(
                relative_path, exclude_paths
            ):
                files.add(relative_path)
            continue
        if path.is_dir():
            for source_file in path.rglob("*.py"):
                relative_path = _normalized_relative_path(repo_root, source_file)
                if not _is_excluded_source_path(relative_path, exclude_paths):
                    files.add(relative_path)
    return files


def _declared_strict_source_files_for_paths(
    repo_root: Path,
    paths: Iterable[str],
) -> _ConfiguredPathSourceFiles:
    files: set[str] = set()
    issues: list[ConfiguredStrictPathIssue] = []
    for path_text in paths:
        normalized_path_text = _normalized_configured_path_text(path_text)
        path = _resolve_configured_path_under_repo(repo_root, path_text)
        if path is None:
            issues.append(
                ConfiguredStrictPathIssue(
                    path=normalized_path_text,
                    reason="outside_repository",
                )
            )
            continue
        if not path.exists():
            issues.append(
                ConfiguredStrictPathIssue(
                    path=normalized_path_text,
                    reason="does_not_exist",
                )
            )
            continue
        if path.is_file():
            if path.suffix == ".py":
                files.add(_normalized_relative_path(repo_root, path))
            else:
                issues.append(
                    ConfiguredStrictPathIssue(
                        path=normalized_path_text,
                        reason="not_python_source_file",
                    )
                )
            continue
        source_files = {
            _normalized_relative_path(repo_root, source_file)
            for source_file in path.rglob("*.py")
        }
        if source_files:
            files.update(source_files)
        else:
            issues.append(
                ConfiguredStrictPathIssue(
                    path=normalized_path_text,
                    reason="no_python_source_files",
                )
            )
    return _ConfiguredPathSourceFiles(
        source_files=files,
        issues=tuple(sorted(issues, key=lambda issue: (issue.path, issue.reason))),
    )


def _resolve_configured_path_under_repo(repo_root: Path, path_text: str) -> Path | None:
    path = Path(path_text)
    resolved = (path if path.is_absolute() else repo_root / path).resolve()
    try:
        resolved.relative_to(repo_root)
    except ValueError:
        return None
    return resolved


def _normalized_relative_path(repo_root: Path, path: Path) -> str:
    return path.relative_to(repo_root).as_posix()


def _normalized_configured_path_text(path_text: str) -> str:
    normalized = path_text.replace("\\", "/")
    while normalized.startswith("./"):
        normalized = normalized[2:]
    return normalized.rstrip("/") or "."


def _is_excluded_source_path(
    relative_path: str,
    exclude_paths: Iterable[str],
) -> bool:
    configured_excludes = tuple(
        _normalized_configured_path_text(path) for path in exclude_paths
    )
    path_parts = relative_path.split("/")
    if any(part in _EXCLUDED_PATH_PARTS for part in path_parts):
        return True
    return any(
        _matches_exclude_pattern(relative_path, pattern)
        for pattern in configured_excludes
    )


def _matches_exclude_pattern(relative_path: str, pattern: str) -> bool:
    if not pattern:
        return False
    if _is_relative_path_at_or_below(relative_path, pattern):
        return True
    candidate_paths = _path_and_parent_paths(relative_path)
    return any(
        Path(candidate).as_posix() != "."
        and (fnmatch.fnmatchcase(candidate, pattern) or Path(candidate).match(pattern))
        for candidate in candidate_paths
    )


def _is_relative_path_at_or_below(relative_path: str, base_path: str) -> bool:
    base = base_path.rstrip("/")
    return relative_path == base or relative_path.startswith(f"{base}/")


def _path_and_parent_paths(relative_path: str) -> tuple[str, ...]:
    path = Path(relative_path)
    parents = tuple(parent.as_posix() for parent in path.parents)
    return (relative_path, *parents)


def _strict_suppression_policy_issues(
    repo_root: Path,
    *,
    pyright_config: Mapping[str, object],
    strict_source_files: Iterable[str],
) -> tuple[StrictSuppressionPolicyIssue, ...]:
    strict_paths = tuple(sorted(strict_source_files))
    issues: list[StrictSuppressionPolicyIssue] = []
    for relative_path in strict_paths:
        issues.extend(_strict_file_suppression_policy_issues(repo_root, relative_path))
    issues.extend(
        _strict_package_ignore_policy_issues(
            repo_root,
            pyright_config=pyright_config,
            strict_source_files=strict_paths,
        )
    )
    return tuple(
        sorted(
            issues,
            key=lambda issue: (issue.path, issue.line, issue.suppression),
        )
    )


def _strict_file_suppression_policy_issues(
    repo_root: Path,
    relative_path: str,
) -> tuple[StrictSuppressionPolicyIssue, ...]:
    path = repo_root / relative_path
    issues: list[StrictSuppressionPolicyIssue] = []
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), 1
    ):
        issues.extend(
            _pyright_directive_policy_issues(
                path=relative_path,
                line_number=line_number,
                line=line,
            )
        )
        type_ignore = _TYPE_IGNORE_RE.search(line)
        if type_ignore is not None:
            issues.append(
                StrictSuppressionPolicyIssue(
                    path=relative_path,
                    line=line_number,
                    suppression=type_ignore.group(0).lstrip("# "),
                    reason=(
                        "strict files must use the machine-checkable "
                        "'# pyright: ignore[reportRule] - reason' format"
                    ),
                )
            )
    return tuple(issues)


def _pyright_directive_policy_issues(
    *,
    path: str,
    line_number: int,
    line: str,
) -> tuple[StrictSuppressionPolicyIssue, ...]:
    directive = _PYRIGHT_DIRECTIVE_RE.search(line)
    if directive is None:
        return ()
    body = directive.group("body").strip()
    ignore = _PYRIGHT_IGNORE_RE.match(body)
    if ignore is not None:
        return _pyright_ignore_policy_issues(
            path=path,
            line_number=line_number,
            rule_list=ignore.group("rule_list"),
            trailing=ignore.group("trailing"),
        )
    issues: list[StrictSuppressionPolicyIssue] = []
    type_checking_downgrade = _FILE_WIDE_TYPE_CHECKING_DOWNGRADE_RE.search(body)
    if type_checking_downgrade is not None:
        issues.append(
            StrictSuppressionPolicyIssue(
                path=path,
                line=line_number,
                suppression=f"pyright: {type_checking_downgrade.group(0).strip()}",
                reason=(
                    "effective strict files may not downgrade file-wide Pyright "
                    "type-checking mode"
                ),
            )
        )
    for suppression in _FILE_WIDE_DIAGNOSTIC_SUPPRESSION_RE.finditer(body):
        issues.append(
            StrictSuppressionPolicyIssue(
                path=path,
                line=line_number,
                suppression=suppression.group(0).strip(),
                reason=(
                    "effective strict files may not use file-wide diagnostic "
                    "severity overrides; use a line-local rule-specific "
                    "suppression with a rationale if unavoidable"
                ),
            )
        )
    return tuple(issues)


def _pyright_ignore_policy_issues(
    *,
    path: str,
    line_number: int,
    rule_list: str | None,
    trailing: str,
) -> tuple[StrictSuppressionPolicyIssue, ...]:
    if rule_list is None:
        return (
            StrictSuppressionPolicyIssue(
                path=path,
                line=line_number,
                suppression="pyright: ignore",
                reason=(
                    "blanket Pyright ignores are not allowed in effective strict files"
                ),
            ),
        )
    rules = tuple(rule.strip() for rule in rule_list.split(","))
    if not rules or any(not rule for rule in rules):
        return (
            StrictSuppressionPolicyIssue(
                path=path,
                line=line_number,
                suppression="pyright: ignore[]",
                reason="Pyright ignore must name at least one diagnostic rule",
            ),
        )
    invalid_rules = tuple(
        rule for rule in rules if _PYRIGHT_RULE_ID_RE.fullmatch(rule) is None
    )
    if invalid_rules:
        return tuple(
            StrictSuppressionPolicyIssue(
                path=path,
                line=line_number,
                suppression=rule,
                reason="Pyright ignore must use explicit report* diagnostic rule identifiers",
            )
            for rule in invalid_rules
        )
    if not trailing.startswith(" - "):
        return tuple(
            StrictSuppressionPolicyIssue(
                path=path,
                line=line_number,
                suppression=rule,
                reason=(
                    "rule-specific Pyright ignore must include ' - ' followed "
                    "by a technical rationale"
                ),
            )
            for rule in rules
        )
    rationale = trailing[3:].strip()
    if not rationale:
        return tuple(
            StrictSuppressionPolicyIssue(
                path=path,
                line=line_number,
                suppression=rule,
                reason="rule-specific Pyright ignore rationale must be non-empty",
            )
            for rule in rules
        )
    if _is_placeholder_rationale(rationale):
        return tuple(
            StrictSuppressionPolicyIssue(
                path=path,
                line=line_number,
                suppression=rule,
                reason=(
                    "rule-specific Pyright ignore rationale must describe the "
                    "concrete typing limitation and runtime safety"
                ),
            )
            for rule in rules
        )
    return ()


def _is_placeholder_rationale(rationale: str) -> bool:
    normalized = re.sub(r"\s+", " ", rationale.strip().lower()).strip(" .:-")
    if not normalized:
        return True
    if normalized in _GENERIC_RATIONALES:
        return True
    tokens = set(re.findall(r"[a-z]+", normalized))
    return bool(tokens & _PLACEHOLDER_RATIONALE_TOKENS)


def _strict_package_ignore_policy_issues(
    repo_root: Path,
    *,
    pyright_config: Mapping[str, object],
    strict_source_files: Iterable[str],
) -> tuple[StrictSuppressionPolicyIssue, ...]:
    strict_paths = tuple(strict_source_files)
    if not strict_paths:
        return ()
    ignore_paths = _optional_string_list(pyright_config, "ignore")
    if not ignore_paths:
        return ()
    line_number = _pyproject_tool_pyright_key_line(repo_root, "ignore")
    issues: list[StrictSuppressionPolicyIssue] = []
    for ignore_path in ignore_paths:
        normalized_ignore = _normalized_configured_path_text(ignore_path)
        matching_strict_file = next(
            (
                strict_path
                for strict_path in strict_paths
                if _configured_path_matches_source_file(
                    strict_path,
                    normalized_ignore,
                )
            ),
            None,
        )
        if matching_strict_file is None:
            continue
        issues.append(
            StrictSuppressionPolicyIssue(
                path="pyproject.toml",
                line=line_number,
                suppression=f"[tool.pyright].ignore={normalized_ignore}",
                reason=(
                    "package-wide Pyright ignore intersects effective strict "
                    f"source file {matching_strict_file}"
                ),
            )
        )
    return tuple(issues)


def _configured_path_matches_source_file(
    relative_path: str,
    configured_path: str,
) -> bool:
    if configured_path in {"", "."}:
        return True
    if _is_relative_path_at_or_below(relative_path, configured_path):
        return True
    return _matches_exclude_pattern(relative_path, configured_path)


def _pyproject_tool_pyright_key_line(repo_root: Path, key: str) -> int:
    in_pyright_table = False
    key_re = re.compile(rf"^\s*{re.escape(key)}\s*=")
    for line_number, line in enumerate(
        (repo_root / "pyproject.toml").read_text(encoding="utf-8").splitlines(),
        1,
    ):
        stripped = line.strip()
        if stripped == "[tool.pyright]":
            in_pyright_table = True
            continue
        if in_pyright_table and stripped.startswith("["):
            break
        if in_pyright_table and key_re.match(line):
            return line_number
    return 1


def _format_human_report(coverage: PyrightStrictCoverage) -> str:
    lines = [
        "Pyright strict source coverage:",
        f"  included source files: {coverage.included_source_file_count}",
        (
            "  declared strict source files: "
            f"{coverage.declared_strict_source_file_count}"
        ),
        (
            "  strict files inside include scope: "
            f"{coverage.strict_source_files_inside_include_count}"
        ),
        (
            "  strict files outside include scope: "
            f"{coverage.strict_source_files_outside_include_count}"
        ),
        (
            "  missing/invalid configured strict paths: "
            f"{len(coverage.missing_or_invalid_configured_strict_paths)}"
        ),
        (
            "  strict coverage: "
            f"{coverage.strict_source_files_inside_include_count}/"
            f"{coverage.included_source_file_count} "
            f"({coverage.strict_source_percent:.1f}%; "
            f"floor {coverage.strict_coverage_floor_percent:.1f}%)"
        ),
    ]
    if coverage.strict_source_files_outside_include:
        lines.append("  outside include:")
        lines.extend(
            f"    - {path}" for path in coverage.strict_source_files_outside_include
        )
    if coverage.missing_or_invalid_configured_strict_paths:
        lines.append("  missing/invalid strict paths:")
        lines.extend(
            f"    - {issue.path} ({issue.reason})"
            for issue in coverage.missing_or_invalid_configured_strict_paths
        )
    if coverage.strict_suppression_policy_issues:
        lines.append("  strict suppression policy issues:")
        lines.extend(
            (f"    - {issue.path}:{issue.line}: {issue.suppression} ({issue.reason})")
            for issue in coverage.strict_suppression_policy_issues
        )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Report Pyright strict source coverage"
    )
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root containing pyproject.toml",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit machine-readable JSON",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help=(
            "Exit non-zero when strict paths are outside include, configured "
            "strict paths are missing/invalid, or the strict floor is unmet"
        ),
    )
    parser.add_argument(
        "--floor-percent",
        type=float,
        default=DEFAULT_STRICT_COVERAGE_FLOOR_PERCENT,
        help="Minimum strict source coverage percentage required by --check",
    )
    args = parser.parse_args()

    coverage = pyright_strict_coverage(
        Path(args.repo_root).resolve(),
        strict_coverage_floor_percent=args.floor_percent,
    )
    if args.json:
        print(json.dumps(coverage.to_payload(), indent=2, sort_keys=True))
    else:
        print(_format_human_report(coverage))
    if args.check and not coverage.passes_policy:
        if not args.json:
            sys.stderr.write("Pyright strict coverage policy failed:\n")
            for message in coverage.policy_failure_messages():
                sys.stderr.write(f"- {message}\n")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
