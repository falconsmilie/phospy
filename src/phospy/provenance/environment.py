"""Environment provenance construction."""

from __future__ import annotations

import platform
from importlib import metadata

from phospy.provenance.models import EnvironmentProvenance

CORE_ENVIRONMENT_DEPENDENCIES = ("numpy", "pandas", "scipy", "scikit-learn")
OPTIONAL_ENVIRONMENT_DEPENDENCIES = ("pyarrow", "openpyxl")
DEFAULT_ENVIRONMENT_DEPENDENCIES = (
    *CORE_ENVIRONMENT_DEPENDENCIES,
    *OPTIONAL_ENVIRONMENT_DEPENDENCIES,
)


def collect_environment_provenance(
    *,
    package_name: str = "phospy",
    dependency_names: tuple[str, ...] = DEFAULT_ENVIRONMENT_DEPENDENCIES,
) -> EnvironmentProvenance:
    """Collect package/python/dependency versions for run provenance."""

    return EnvironmentProvenance(
        package_name=package_name,
        package_version=_distribution_version(package_name) or "unknown",
        python_version=platform.python_version(),
        dependency_versions={
            dependency: _distribution_version(dependency)
            for dependency in dependency_names
        },
        platform=_platform_provenance(),
    )


def _distribution_version(distribution_name: str) -> str | None:
    try:
        return metadata.version(distribution_name)
    except metadata.PackageNotFoundError:
        return None


def _platform_provenance() -> dict[str, str]:
    return {
        "platform": _normalize_platform_value(platform.platform()),
        "system": _normalize_platform_value(platform.system()),
        "release": _normalize_platform_value(platform.release()),
        "version": _normalize_platform_value(platform.version()),
        "machine": _normalize_platform_value(platform.machine()),
        "processor": _normalize_platform_value(platform.processor()),
        "python_implementation": _normalize_platform_value(
            platform.python_implementation()
        ),
    }


def _normalize_platform_value(value: str) -> str:
    normalized = str(value).strip()
    if normalized:
        return normalized
    return "unknown"


__all__ = [
    "CORE_ENVIRONMENT_DEPENDENCIES",
    "DEFAULT_ENVIRONMENT_DEPENDENCIES",
    "OPTIONAL_ENVIRONMENT_DEPENDENCIES",
    "collect_environment_provenance",
]
