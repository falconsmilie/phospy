"""PhosPy package metadata."""

from __future__ import annotations

import tomllib as _tomllib
from collections.abc import Mapping as _Mapping
from importlib import metadata as _metadata
from pathlib import Path as _Path


def _source_project_version() -> str | None:
    """Return the version from this checkout's fixed ``src`` layout."""

    package_directory = _Path(__file__).resolve().parent
    source_directory = package_directory.parent
    if package_directory.name != "phospy" or source_directory.name != "src":
        return None

    pyproject = source_directory.parent / "pyproject.toml"
    try:
        payload = _tomllib.loads(pyproject.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, _tomllib.TOMLDecodeError):
        return None
    project = payload.get("project")
    if not isinstance(project, _Mapping):
        return None
    if _normalized_string(project.get("name")) != "phospy":
        return None
    return _normalized_string(project.get("version"))


def _installed_distribution_version() -> str | None:
    try:
        return _metadata.version("phospy")
    except _metadata.PackageNotFoundError:  # pragma: no cover - source tree fallback
        return None


def _normalized_string(value: object) -> str | None:
    if not isinstance(value, str):
        return None
    normalized = value.strip()
    if normalized:
        return normalized
    return None


__version__ = (
    _installed_distribution_version() or _source_project_version() or "0+unknown"
)

__all__: list[str] = []
