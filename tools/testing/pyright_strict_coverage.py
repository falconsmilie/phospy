from __future__ import annotations

import argparse
import fnmatch
import json
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
    return PyrightStrictCoverage(
        included_source_files=tuple(sorted(included)),
        declared_strict_source_files=tuple(sorted(strict.source_files)),
        strict_source_files_inside_include=tuple(sorted(strict_inside_include)),
        strict_source_files_outside_include=tuple(sorted(strict_outside_include)),
        missing_or_invalid_configured_strict_paths=strict.issues,
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
