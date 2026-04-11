from __future__ import annotations

from ...errors import InputCompatibilityError
from ...prediction.contracts import EnsemblePredictorContract


def validate_ensemble_predictor(
    ensemble_predictor: EnsemblePredictorContract | None,
) -> EnsemblePredictorContract | None:
    """Validate an injected ensemble predictor against the prediction contract."""

    if ensemble_predictor is None:
        return None
    if not isinstance(ensemble_predictor, EnsemblePredictorContract):
        msg = (
            "ensemble_predictor must implement EnsemblePredictorContract with "
            "predict_kinase(..., sampling_session=...)"
        )
        raise InputCompatibilityError(msg)
    return ensemble_predictor


__all__ = ["validate_ensemble_predictor"]
