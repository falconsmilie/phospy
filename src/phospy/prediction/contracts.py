from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from ..validation.requests import PredictionRequest
    from .execution import (
        KinasePredictionBatch,
        PredictionSamplingSession,
        PredictionTraceState,
    )


class EnsemblePredictorContract(ABC):
    """Explicit contract for one kinase prediction executor."""

    @abstractmethod
    def predict_kinase(
        self,
        *,
        kinase: str,
        substrates: list[str],
        feature_mat: pd.DataFrame,
        request: PredictionRequest,
        trace_state: PredictionTraceState,
        sampling_session: PredictionSamplingSession,
    ) -> KinasePredictionBatch:
        raise NotImplementedError


__all__ = ["EnsemblePredictorContract"]
