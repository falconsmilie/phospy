"""Compatibility helpers for retired ``phospy.api`` export routes."""

from __future__ import annotations

import importlib
import warnings
from types import ModuleType
from typing import Any

from phospy._api_inventory import ADVANCED_PUBLIC_API


def deprecated_advanced_export(name: str, *, old_module: str) -> Any:
    """Resolve an advanced symbol from ``phospy.advanced`` with a warning."""

    if name.startswith("__"):
        raise AttributeError(f"module {old_module!r} has no attribute {name!r}")
    advanced = importlib.import_module("phospy.advanced")
    if name not in ADVANCED_PUBLIC_API:
        raise AttributeError(f"module {old_module!r} has no attribute {name!r}")
    _warn(
        old_module=old_module,
        name=name,
        replacement_module="phospy.advanced",
    )
    return getattr(advanced, name)


def deprecated_advanced_result_export(name: str, *, old_module: str) -> Any:
    """Resolve a result inspection symbol from ``phospy.advanced.results``."""

    if name.startswith("__"):
        raise AttributeError(f"module {old_module!r} has no attribute {name!r}")
    advanced_results = importlib.import_module("phospy.advanced.results")
    if name not in getattr(advanced_results, "__all__", ()):
        raise AttributeError(f"module {old_module!r} has no attribute {name!r}")
    _warn(
        old_module=old_module,
        name=name,
        replacement_module="phospy.advanced.results",
    )
    return getattr(advanced_results, name)


def deprecated_request_export(name: str, *, old_module: str) -> Any:
    """Resolve a retired request-wrapper export with replacement guidance."""

    if name.startswith("__"):
        raise AttributeError(f"module {old_module!r} has no attribute {name!r}")
    if name in ADVANCED_PUBLIC_API:
        return deprecated_advanced_export(name, old_module=old_module)
    owner = importlib.import_module("phospy.contracts.requests")
    if _module_exports(owner, name):
        _warn(
            old_module=old_module,
            name=name,
            replacement_module="phospy.contracts.requests",
            supported=False,
        )
        return getattr(owner, name)
    raise AttributeError(f"module {old_module!r} has no attribute {name!r}")


def deprecated_config_export(name: str, *, old_module: str) -> Any:
    """Resolve a retired config-wrapper export with replacement guidance."""

    if name.startswith("__"):
        raise AttributeError(f"module {old_module!r} has no attribute {name!r}")
    advanced_configs = importlib.import_module("phospy.advanced.configs")
    if name in getattr(advanced_configs, "__all__", ()):
        _warn(
            old_module=old_module,
            name=name,
            replacement_module="phospy.advanced.configs",
        )
        return getattr(advanced_configs, name)

    owner_module_name = old_module.replace(
        "phospy.api.configs",
        "phospy.contracts.configs",
        1,
    )
    owner = importlib.import_module(owner_module_name)
    if _module_exports(owner, name):
        _warn(
            old_module=old_module,
            name=name,
            replacement_module=owner_module_name,
            supported=False,
        )
        return getattr(owner, name)
    aggregate_owner = importlib.import_module("phospy.contracts.configs")
    if _module_exports(aggregate_owner, name):
        _warn(
            old_module=old_module,
            name=name,
            replacement_module="phospy.contracts.configs",
            supported=False,
        )
        return getattr(aggregate_owner, name)
    raise AttributeError(f"module {old_module!r} has no attribute {name!r}")


def _module_exports(module: ModuleType, name: str) -> bool:
    exported = getattr(module, "__all__", ())
    return name in exported or hasattr(module, name)


def _warn(
    *,
    old_module: str,
    name: str,
    replacement_module: str,
    supported: bool = True,
) -> None:
    stability_note = (
        "an advanced API route" if supported else "an unsupported compatibility route"
    )
    warnings.warn(
        f"{old_module}.{name} is deprecated as {stability_note}; "
        f"use `from {replacement_module} import {name}`.",
        DeprecationWarning,
        stacklevel=4,
    )


__all__ = [
    "deprecated_advanced_export",
    "deprecated_advanced_result_export",
    "deprecated_config_export",
    "deprecated_request_export",
]
