"""Stable preprocessing config facade with compatibility adapters."""

from __future__ import annotations

from phospy.api._compat import deprecated_config_export
from phospy.contracts.configs.preprocessing.localisation import (
    DatasetLocalisationConfig,
)

__all__ = ["DatasetLocalisationConfig"]


def __getattr__(name: str) -> object:
    return deprecated_config_export(name, old_module=__name__)
