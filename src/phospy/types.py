from __future__ import annotations

from typing import Literal

DegenerateProbabilityPolicy = Literal["uniform", "error"]
PredictionSvmMode = Literal["default", "r_parity"]

__all__ = ["DegenerateProbabilityPolicy", "PredictionSvmMode"]
