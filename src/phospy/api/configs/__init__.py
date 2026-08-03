# pyright: reportUnsupportedDunderAll=false
# ruff: noqa: F401
"""Stable public config facade with compatibility adapters.

Only stable configuration objects are exported from this module. Advanced
configuration and policy objects are supported from ``phospy.advanced.configs``.
Historical imports from ``phospy.api.configs`` remain available during
migration through module ``__getattr__`` and emit ``DeprecationWarning``.
"""

from __future__ import annotations

from phospy._api_inventory import STABLE_CONFIG_API
from phospy.api._compat import deprecated_config_export
from phospy.contracts.configs import (
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    EnrichmentConfig,
)

__all__ = STABLE_CONFIG_API


def __getattr__(name: str) -> object:
    return deprecated_config_export(name, old_module=__name__)
