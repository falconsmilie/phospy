"""Identity transformer for the supported analysis-ready builder lane."""

from __future__ import annotations

import pandas as pd

from phospy.science.transformations.contracts import TransformationResult
from phospy.science.transformations.models import IntensityScaleState


class IdentityTransformer:
    """Pass through matrix values without establishing intensity scale state."""

    preserves_input_scale_state = True
    changes_numeric_values = False
    requires_established_input_state = False

    def run(
        self,
        phospho: pd.DataFrame,
        total: pd.DataFrame | None = None,
    ) -> TransformationResult:
        return TransformationResult(
            phospho=phospho,
            total=total,
            state=IntensityScaleState.raw(has_total_matrix=total is not None),
        )
