"""Intensity-scale-state validator."""

from __future__ import annotations

from phospy.science.transformations.models import IntensityScaleState
from phospy.science.transformations.state_coherence import (
    require_intensity_scale_state_coherence,
)


class IntensityScaleStateValidator:
    """Validate explicit intensity-scale state coherence."""

    def run(
        self,
        intensity_scale_state: IntensityScaleState,
        *,
        has_total_matrix: bool,
        require_established: bool = False,
    ) -> IntensityScaleState:
        return require_intensity_scale_state_coherence(
            intensity_scale_state,
            has_total_matrix=has_total_matrix,
            require_established=require_established,
        )
