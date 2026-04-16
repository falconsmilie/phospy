from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from ..validation.requests.prediction import PredictionRequest
    from .execution import (
        KinasePredictionBatch,
        PredictionSamplingSession,
        PredictionTraceState,
    )


CUSTOM_PREDICTOR_SCORE_INDEX_CONTRACT: str = "label_aligned"
"""Score-index contract for custom predictor batch outputs.

`label_aligned` means:
- `score_index` must contain the same phosphosite labels as `feature_mat.index`.
- Reordered labels are allowed and are aligned by label before aggregation.
- Missing, unexpected, duplicate, or length-mismatched labels are rejected.
"""


class EnsemblePredictorContract(ABC):
    """Explicit contract for one kinase prediction executor.

    Custom implementations are extension points. `predict_kinase(...)` must
    return a `KinasePredictionBatch` that satisfies all the following:

    - `batch.kinase` equals the requested `kinase`.
    - `batch.score_values` is a 1D numeric array with one value per phosphosite
      row in `feature_mat`.
    - all score values are finite (`NaN` and `+/-inf` are rejected).
    - `batch.score_index` follows `CUSTOM_PREDICTOR_SCORE_INDEX_CONTRACT`.
    """

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

    def clear_cache(self) -> None:
        """Release any optional temporary caches held across kinase runs."""

        return None


__all__ = [
    "CUSTOM_PREDICTOR_SCORE_INDEX_CONTRACT",
    "EnsemblePredictorContract",
]
