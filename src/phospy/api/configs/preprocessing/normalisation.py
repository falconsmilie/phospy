"""Deprecated compatibility adapter for advanced normalisation configs."""

from __future__ import annotations

from phospy.api._compat import deprecated_config_export

__all__: list[str] = []


def __getattr__(name: str) -> object:
    return deprecated_config_export(name, old_module=__name__)
