"""Compatibility helpers for retired ``phospy.api`` export routes."""

from __future__ import annotations

import importlib
from functools import lru_cache
from types import ModuleType
from typing import Any, NamedTuple

from phospy._deprecations import (
    RetainedDeprecation,
    api_compatibility_deprecation_id,
    compatibility_deprecation_record,
    compatibility_deprecation_records,
    warn_deprecated,
)


class CompatibilityExport(NamedTuple):
    """Documented metadata for one deprecated compatibility route."""

    old_module: str
    name: str
    owner_module: str
    replacement_module: str
    introduced_version: str
    planned_removal_version: str
    stability: str


def deprecated_advanced_export(name: str, *, old_module: str) -> Any:
    """Resolve an advanced symbol from ``phospy.advanced`` with a warning."""

    if name.startswith("__"):
        raise AttributeError(f"module {old_module!r} has no attribute {name!r}")
    compat = compatibility_export(old_module=old_module, name=name)
    advanced = importlib.import_module("phospy.advanced")
    if compat is None or compat.owner_module != "phospy.advanced":
        raise AttributeError(f"module {old_module!r} has no attribute {name!r}")
    _warn(compat)
    return getattr(advanced, name)


def deprecated_advanced_result_export(name: str, *, old_module: str) -> Any:
    """Resolve a result inspection symbol from ``phospy.advanced.results``."""

    if name.startswith("__"):
        raise AttributeError(f"module {old_module!r} has no attribute {name!r}")
    compat = compatibility_export(old_module=old_module, name=name)
    advanced_results = importlib.import_module("phospy.advanced.results")
    if compat is None or compat.owner_module != "phospy.advanced.results":
        raise AttributeError(f"module {old_module!r} has no attribute {name!r}")
    _warn(compat)
    return getattr(advanced_results, name)


def deprecated_request_export(name: str, *, old_module: str) -> Any:
    """Resolve a retired request-wrapper export with replacement guidance."""

    if name.startswith("__"):
        raise AttributeError(f"module {old_module!r} has no attribute {name!r}")
    compat = compatibility_export(old_module=old_module, name=name)
    if compat is not None:
        owner = importlib.import_module(compat.owner_module)
        if _module_exports(owner, name):
            _warn(compat)
            return getattr(owner, name)
    raise AttributeError(f"module {old_module!r} has no attribute {name!r}")


def deprecated_config_export(name: str, *, old_module: str) -> Any:
    """Resolve a retired config-wrapper export with replacement guidance."""

    if name.startswith("__"):
        raise AttributeError(f"module {old_module!r} has no attribute {name!r}")
    compat = compatibility_export(old_module=old_module, name=name)
    if compat is not None:
        owner = importlib.import_module(compat.owner_module)
        if _module_exports(owner, name):
            _warn(compat)
            return getattr(owner, name)
    raise AttributeError(f"module {old_module!r} has no attribute {name!r}")


def compatibility_export(*, old_module: str, name: str) -> CompatibilityExport | None:
    """Return metadata for a retained deprecated export route."""

    record = compatibility_deprecation_record(old_module=old_module, name=name)
    if record is None:
        return None
    return _compatibility_export_from_record(record)


def compatibility_exports() -> tuple[CompatibilityExport, ...]:
    """Return all documented deprecated export routes."""

    return tuple(_compatibility_export_map().values())


@lru_cache(maxsize=1)
def _compatibility_export_map() -> dict[tuple[str, str], CompatibilityExport]:
    return {
        (compat.old_module, compat.name): compat
        for compat in (
            _compatibility_export_from_record(record)
            for record in compatibility_deprecation_records()
        )
    }


def _module_exports(module: ModuleType, name: str) -> bool:
    exported = getattr(module, "__all__", ())
    return name in exported or hasattr(module, name)


def _compatibility_export_from_record(
    record: RetainedDeprecation,
) -> CompatibilityExport:
    if record.deprecated_module is None or record.deprecated_name is None:
        raise RuntimeError(
            "API compatibility records must carry deprecated module/name metadata"
        )
    return CompatibilityExport(
        old_module=record.deprecated_module,
        name=record.deprecated_name,
        owner_module=record.owner_module,
        replacement_module=record.replacement_module,
        introduced_version=record.introduced_version,
        planned_removal_version=record.planned_removal_version,
        stability=record.stability,
    )


def _warn(
    compat: CompatibilityExport,
) -> None:
    warn_deprecated(
        api_compatibility_deprecation_id(
            old_module=compat.old_module,
            name=compat.name,
        ),
        stacklevel=4,
    )


__all__ = [
    "CompatibilityExport",
    "compatibility_export",
    "compatibility_exports",
    "deprecated_advanced_export",
    "deprecated_advanced_result_export",
    "deprecated_config_export",
    "deprecated_request_export",
]
