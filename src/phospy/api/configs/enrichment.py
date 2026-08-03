"""Stable enrichment config facade with compatibility adapters."""

from __future__ import annotations

from phospy.api._compat import deprecated_config_export
from phospy.contracts.configs.enrichment import EnrichmentConfig

__all__ = ["EnrichmentConfig"]


def __getattr__(name: str) -> object:
    return deprecated_config_export(name, old_module=__name__)
