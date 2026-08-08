"""Compatibility helpers for retired ``phospy.api`` export routes."""

from __future__ import annotations

import importlib
import warnings
from functools import lru_cache
from types import ModuleType
from typing import Any, NamedTuple

from phospy._api_inventory import (
    ADVANCED_CONFIG_API,
    ADVANCED_PUBLIC_API,
    ADVANCED_RESULT_API,
    API_COMPATIBILITY_INTRODUCED_VERSION,
    API_COMPATIBILITY_PLANNED_REMOVAL_VERSION,
    COMPATIBILITY_CONFIG_MODULES,
    CONFIG_COMPATIBILITY_ADVANCED_ROUTE_OVERRIDES,
    REQUEST_COMPATIBILITY_ADVANCED_API,
    STABLE_CONFIG_API,
    STABLE_REQUEST_API,
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
    if compat is None or name not in ADVANCED_PUBLIC_API:
        raise AttributeError(f"module {old_module!r} has no attribute {name!r}")
    _warn(compat)
    return getattr(advanced, name)


def deprecated_advanced_result_export(name: str, *, old_module: str) -> Any:
    """Resolve a result inspection symbol from ``phospy.advanced.results``."""

    if name.startswith("__"):
        raise AttributeError(f"module {old_module!r} has no attribute {name!r}")
    compat = compatibility_export(old_module=old_module, name=name)
    advanced_results = importlib.import_module("phospy.advanced.results")
    if compat is None or name not in getattr(advanced_results, "__all__", ()):
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

    return _compatibility_export_map().get((old_module, name))


def compatibility_exports() -> tuple[CompatibilityExport, ...]:
    """Return all documented deprecated export routes."""

    return tuple(_compatibility_export_map().values())


@lru_cache(maxsize=1)
def _compatibility_export_map() -> dict[tuple[str, str], CompatibilityExport]:
    entries: dict[tuple[str, str], CompatibilityExport] = {}

    for name in ADVANCED_PUBLIC_API:
        _register(
            entries,
            old_module="phospy.api",
            name=name,
            owner_module="phospy.advanced",
            replacement_module="phospy.advanced",
            stability="advanced",
        )

    for name in REQUEST_COMPATIBILITY_ADVANCED_API:
        _register(
            entries,
            old_module="phospy.api.requests",
            name=name,
            owner_module="phospy.advanced",
            replacement_module="phospy.advanced",
            stability="advanced",
        )

    advanced_config_routes = _advanced_config_routes()
    for name in ADVANCED_CONFIG_API:
        for old_module in advanced_config_routes.get(name, ()):
            _register(
                entries,
                old_module=old_module,
                name=name,
                owner_module="phospy.advanced.configs",
                replacement_module="phospy.advanced.configs",
                stability="advanced",
            )

    for name in ADVANCED_RESULT_API:
        _register(
            entries,
            old_module="phospy.api.results",
            name=name,
            owner_module="phospy.advanced.results",
            replacement_module="phospy.advanced.results",
            stability="advanced",
        )

    _register_contract_request_exports(entries)
    _register_contract_config_exports(entries)
    return entries


def _register_contract_request_exports(
    entries: dict[tuple[str, str], CompatibilityExport],
) -> None:
    owner_module = "phospy.contracts.requests"
    owner = importlib.import_module(owner_module)
    for name in getattr(owner, "__all__", ()):
        if name in STABLE_REQUEST_API:
            continue
        _register(
            entries,
            old_module="phospy.api.requests",
            name=name,
            owner_module=owner_module,
            replacement_module=owner_module,
            stability="unsupported",
        )


def _register_contract_config_exports(
    entries: dict[tuple[str, str], CompatibilityExport],
) -> None:
    for old_module in COMPATIBILITY_CONFIG_MODULES:
        owner_module = _config_owner_module(old_module)
        owner = importlib.import_module(owner_module)
        for name in getattr(owner, "__all__", ()):
            if name in STABLE_CONFIG_API or name in ADVANCED_CONFIG_API:
                continue
            _register(
                entries,
                old_module=old_module,
                name=name,
                owner_module=owner_module,
                replacement_module=owner_module,
                stability="unsupported",
            )


def _advanced_config_routes() -> dict[str, tuple[str, ...]]:
    routes: dict[str, list[str]] = {name: [] for name in ADVANCED_CONFIG_API}
    for old_module in COMPATIBILITY_CONFIG_MODULES:
        owner_module = _config_owner_module(old_module)
        owner = importlib.import_module(owner_module)
        owner_exports = set(getattr(owner, "__all__", ()))
        for name in ADVANCED_CONFIG_API:
            if old_module == "phospy.api.configs" or name in owner_exports:
                routes[name].append(old_module)
    for name, route_overrides in CONFIG_COMPATIBILITY_ADVANCED_ROUTE_OVERRIDES.items():
        routes.setdefault(name, [])
        for old_module in route_overrides:
            if old_module not in routes[name]:
                routes[name].append(old_module)
    return {name: tuple(route_names) for name, route_names in routes.items()}


def _config_owner_module(old_module: str) -> str:
    return old_module.replace(
        "phospy.api.configs",
        "phospy.contracts.configs",
        1,
    )


def _register(
    entries: dict[tuple[str, str], CompatibilityExport],
    *,
    old_module: str,
    name: str,
    owner_module: str,
    replacement_module: str,
    stability: str,
) -> None:
    entries.setdefault(
        (old_module, name),
        CompatibilityExport(
            old_module=old_module,
            name=name,
            owner_module=owner_module,
            replacement_module=replacement_module,
            introduced_version=API_COMPATIBILITY_INTRODUCED_VERSION,
            planned_removal_version=API_COMPATIBILITY_PLANNED_REMOVAL_VERSION,
            stability=stability,
        ),
    )


def _module_exports(module: ModuleType, name: str) -> bool:
    exported = getattr(module, "__all__", ())
    return name in exported or hasattr(module, name)


def _warn(
    compat: CompatibilityExport,
) -> None:
    stability_note = (
        "an advanced API route"
        if compat.stability == "advanced"
        else "an unsupported compatibility route"
    )
    warnings.warn(
        f"{compat.old_module}.{compat.name} is deprecated as {stability_note}; "
        f"use `from {compat.replacement_module} import {compat.name}`. "
        f"This compatibility route was introduced in PhosPy "
        f"{compat.introduced_version} and is planned for removal in PhosPy "
        f"{compat.planned_removal_version}.",
        DeprecationWarning,
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
