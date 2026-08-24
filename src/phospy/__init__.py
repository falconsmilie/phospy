"""Curated convenience entrypoints for the PhosPy package root."""

from __future__ import annotations

import tomllib as _tomllib
from collections.abc import Mapping as _Mapping
from importlib import metadata as _metadata
from pathlib import Path as _Path

from phospy.api import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset
from phospy.workflows.differential.public import DifferentialAnalysisWorkflow
from phospy.workflows.kinase.public import KinaseWorkflow
from phospy.workflows.signalome.public import SignalomeWorkflow


def _source_project_version() -> str | None:
    """Return the authoritative source-tree project version when available."""

    for parent in _Path(__file__).resolve().parents:
        pyproject = parent / "pyproject.toml"
        if not pyproject.is_file():
            continue
        try:
            payload = _tomllib.loads(pyproject.read_text(encoding="utf-8"))
        except Exception:  # pragma: no cover - defensive source-tree fallback
            return None
        project = payload.get("project")
        if not isinstance(project, _Mapping):
            return None
        if _normalized(project.get("name")) != "phospy":
            return None
        return _normalized(project.get("version"))
    return None


def _installed_distribution_version() -> str | None:
    try:
        return _metadata.version("phospy")
    except _metadata.PackageNotFoundError:  # pragma: no cover - source tree fallback
        return None


def _normalized(value: object) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip()
    if normalized:
        return normalized
    return None


__version__ = (
    _source_project_version() or _installed_distribution_version() or "0+unknown"
)

__all__ = [
    "AnalysisReadyDatasetBuilder",
    "AnalysisReadyPhosphoDataset",
    "DifferentialAnalysisWorkflow",
    "KinaseWorkflow",
    "SignalomeWorkflow",
]
