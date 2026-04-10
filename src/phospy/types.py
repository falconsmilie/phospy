from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Literal

PredictionSvmMode = Literal["default", "r_parity"]
PredictionTraceLevel = Literal["none", "summary", "full"]
PredictionTraceFormat = Literal["csv", "parquet"]

KinaseSubstrateMap = Mapping[str, Sequence[str]]
KinaseMotifSequenceMap = Mapping[str, Sequence[str]]

__all__ = [
    "KinaseMotifSequenceMap",
    "KinaseSubstrateMap",
    "PredictionSvmMode",
    "PredictionTraceFormat",
    "PredictionTraceLevel",
]
