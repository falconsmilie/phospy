"""Environment provenance construction."""

from __future__ import annotations

import platform
from importlib import metadata

from phospy.provenance.models import EnvironmentProvenance

DEFAULT_ENVIRONMENT_DEPENDENCIES = ("numpy", "pandas", "scikit-learn")


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
    )


def _distribution_version(distribution_name: str) -> str | None:
    try:
        return metadata.version(distribution_name)
    except metadata.PackageNotFoundError:
        return None


__all__ = [
    "DEFAULT_ENVIRONMENT_DEPENDENCIES",
    "collect_environment_provenance",
]
