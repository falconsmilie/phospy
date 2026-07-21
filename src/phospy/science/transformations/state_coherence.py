"""Science-owned intensity-scale-state coherence checks."""

from __future__ import annotations

from typing import cast

from phospy.errors.validation import TransformationValidationError
from phospy.science.transformations.models import IntensityScaleState


def require_intensity_scale_state_coherence(
    intensity_scale_state: IntensityScaleState,
    *,
    has_total_matrix: bool,
    require_established: bool = False,
) -> IntensityScaleState:
    """Return a coherent intensity-scale state or raise a validation error."""

    if not isinstance(cast(object, intensity_scale_state), IntensityScaleState):
        raise TransformationValidationError(
            "dataset.intensity_scale_state must be an IntensityScaleState instance"
        )
    if require_established and not intensity_scale_state.is_established:
        raise TransformationValidationError(
            "dataset.intensity_scale_state must be established through a "
            "supported PhosPy path; use AnalysisReadyDatasetBuilder or a "
            "supported transformer/bundle reconstruction path"
        )
    if require_established and intensity_scale_state.quantity is None:
        raise TransformationValidationError(
            "dataset.intensity_scale_state must have established quantitative "
            "meaning provenance"
        )
    if (
        require_established
        and intensity_scale_state.quantitative_meaning_provenance is None
    ):
        raise TransformationValidationError(
            "dataset.intensity_scale_state must carry quantitative meaning "
            "provenance separate from intensity-scale establishment provenance"
        )
    if has_total_matrix and intensity_scale_state.total is None:
        raise TransformationValidationError(
            "intensity_scale_state.total is required when dataset.total is provided"
        )
    if not has_total_matrix and intensity_scale_state.total is not None:
        raise TransformationValidationError(
            "intensity_scale_state.total must be None when dataset.total is absent"
        )
    if (
        intensity_scale_state.total is not None
        and intensity_scale_state.total.kind is not intensity_scale_state.phospho.kind
    ):
        raise TransformationValidationError(
            "phospho and total matrices must share one intensity scale kind"
        )
    return intensity_scale_state


__all__ = ["require_intensity_scale_state_coherence"]
