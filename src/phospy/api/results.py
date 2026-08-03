# pyright: reportUnsupportedDunderAll=false
# ruff: noqa: F401
"""Stable public result models."""

from __future__ import annotations

from phospy._api_inventory import STABLE_RESULT_API
from phospy.api._compat import deprecated_advanced_result_export
from phospy.contracts.results import (
    DifferentialAnalysisResult,
    EnrichmentResultRecord,
    EnrichmentWorkflowResult,
    KinaseActivityResult,
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
    PhosphositeImportResult,
    ResultCaveat,
    SignalomeWorkflowResult,
)

__all__ = STABLE_RESULT_API


def __getattr__(name: str) -> object:
    return deprecated_advanced_result_export(name, old_module=__name__)
