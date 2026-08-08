from __future__ import annotations

import argparse
import json
import tomllib
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass
from pathlib import Path

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
class PyrightStrictCoverage:
    included_source_files: int
    strict_source_files: int

    @property
    def strict_source_fraction(self) -> float:
        if self.included_source_files == 0:
            return 0.0
        return self.strict_source_files / self.included_source_files

    @property
    def strict_source_percent(self) -> float:
        return self.strict_source_fraction * 100.0

    def to_payload(self) -> dict[str, object]:
        payload = asdict(self)
        payload["strict_source_fraction"] = self.strict_source_fraction
        payload["strict_source_percent"] = self.strict_source_percent
        return payload


def pyright_strict_coverage(repo_root: Path) -> PyrightStrictCoverage:
    pyright_config = _pyright_config(repo_root)
    included = _source_files_for_paths(
        repo_root, _string_list(pyright_config, "include")
    )
    strict = _source_files_for_paths(repo_root, _string_list(pyright_config, "strict"))
    return PyrightStrictCoverage(
        included_source_files=len(included),
        strict_source_files=len(strict & included),
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


def _source_files_for_paths(repo_root: Path, paths: Iterable[str]) -> set[str]:
    files: set[str] = set()
    for path_text in paths:
        path = repo_root / path_text
        if path.is_file():
            if path.suffix == ".py":
                files.add(_normalized_relative_path(repo_root, path))
            continue
        if path.is_dir():
            for source_file in path.rglob("*.py"):
                relative_path = _normalized_relative_path(repo_root, source_file)
                if not _is_excluded_source_path(relative_path):
                    files.add(relative_path)
    return files


def _normalized_relative_path(repo_root: Path, path: Path) -> str:
    return path.relative_to(repo_root).as_posix()


def _is_excluded_source_path(relative_path: str) -> bool:
    path_parts = relative_path.split("/")
    return any(part in _EXCLUDED_PATH_PARTS for part in path_parts)


def main() -> None:
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
    args = parser.parse_args()

    coverage = pyright_strict_coverage(Path(args.repo_root).resolve())
    if args.json:
        print(json.dumps(coverage.to_payload(), indent=2, sort_keys=True))
        return
    print(
        "Pyright strict source coverage: "
        f"{coverage.strict_source_files}/{coverage.included_source_files} "
        f"({coverage.strict_source_percent:.1f}%)"
    )


if __name__ == "__main__":
    main()
