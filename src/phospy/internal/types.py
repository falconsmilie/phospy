from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Literal

PredictionSvmMode = Literal["default", "r_parity"]
PredictionTraceLevel = Literal["none", "summary", "full"]
PredictionTraceFormat = Literal["csv", "parquet"]
DuplicateSiteStrategy = Literal[
    "max_mean_signal",
    "first",
    "aggregate_mean",
    "aggregate_median",
    "error",
]
KinaseProfileMissingValueStrategy = Literal[
    "propagate_any_missing",
    "median_skipna",
]
SignalomeModuleSelectionStrategy = Literal[
    "correlation_thresholds",
    "single_module",
]

KinaseSubstrateMap = Mapping[str, Sequence[str]]
KinaseMotifSequenceMap = Mapping[str, Sequence[str]]

__all__ = [
    "DuplicateSiteStrategy",
    "KinaseMotifSequenceMap",
    "KinaseProfileMissingValueStrategy",
    "KinaseSubstrateMap",
    "PredictionSvmMode",
    "PredictionTraceFormat",
    "PredictionTraceLevel",
    "SignalomeModuleSelectionStrategy",
]
