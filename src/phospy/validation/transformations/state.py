"""Intensity-scale-state validator."""

from __future__ import annotations

from phospy.errors.validation import TransformationValidationError
from phospy.transformations.models import IntensityScaleState


class IntensityScaleStateValidator:
    """Validate explicit intensity-scale state coherence."""

    def run(
        self,
        intensity_scale_state: IntensityScaleState,
        *,
        has_total_matrix: bool,
        require_established: bool = False,
    ) -> IntensityScaleState:
        if not isinstance(intensity_scale_state, IntensityScaleState):
            raise TransformationValidationError(
                "dataset.intensity_scale_state must be an IntensityScaleState instance"
            )
        if require_established and not intensity_scale_state.is_established:
            raise TransformationValidationError(
                "dataset.intensity_scale_state must be established through a "
                "supported PhosPy path; use AnalysisReadyDatasetBuilder or a "
                "supported transformer/bundle reconstruction path"
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
            and intensity_scale_state.total.kind
            is not intensity_scale_state.phospho.kind
        ):
            raise TransformationValidationError(
                "phospho and total matrices must share one intensity scale kind"
            )
        return intensity_scale_state
