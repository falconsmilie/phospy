"""Prediction domain.

This package owns kinase prediction engines, scoring components, prediction
execution, and prediction result models.

Only the stable public prediction API is exported from this package module.
Advanced and internal helpers must be imported from concrete submodules.
"""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS: dict[str, tuple[str, str]] = {
    "KinaseMotifScorer": (".motif_scoring", "KinaseMotifScorer"),
    "KinasePredictor": (".engines", "KinasePredictor"),
    "KinasePredictionResult": (".results", "KinasePredictionResult"),
    "KinaseProfilePolicy": (".profiles", "KinaseProfilePolicy"),
    "KinaseScorer": (".scoring", "KinaseScorer"),
    "KinaseScoringResult": (".scoring", "KinaseScoringResult"),
    "MotifScoringResult": (".motif_scoring", "MotifScoringResult"),
    "PredMatResult": (".results", "PredMatResult"),
}

__all__ = [
    "KinaseMotifScorer",
    "KinasePredictor",
    "KinasePredictionResult",
    "KinaseProfilePolicy",
    "KinaseScorer",
    "KinaseScoringResult",
    "MotifScoringResult",
    "PredMatResult",
]


def __getattr__(name: str) -> Any:
    if name not in _EXPORTS:
        msg = f"module {__name__!r} has no attribute {name!r}"
        raise AttributeError(msg)
    module_name, attribute_name = _EXPORTS[name]
    module = import_module(module_name, __name__)
    value = getattr(module, attribute_name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
