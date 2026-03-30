from __future__ import annotations

from typing import Literal

PredictionSvmMode = Literal["default", "r_parity"]
PredictionTraceLevel = Literal["none", "summary", "full"]
PredictionTraceFormat = Literal["csv", "parquet"]

__all__ = [
    "PredictionSvmMode",
    "PredictionTraceFormat",
    "PredictionTraceLevel",
]
