"""Identity transformer for the supported analysis-ready builder lane."""

from __future__ import annotations

import pandas as pd

from phospy.science.transformations._authority import (
    _identity_transformer_establishment_authority,
)
from phospy.science.transformations.contracts import TransformationResult
from phospy.science.transformations.models import IntensityScaleState


class IdentityTransformer:
    """Establish linear state without changing matrix values.

    This is a strict pass-through establisher used when inputs are already
    donor-aligned analysis-ready matrices.
    """

    def run(
        self,
        phospho: pd.DataFrame,
        total: pd.DataFrame | None = None,
    ) -> TransformationResult:
        return TransformationResult(
            phospho=phospho,
            total=total,
            state=IntensityScaleState.established_raw(
                has_total_matrix=total is not None,
                _authority=_identity_transformer_establishment_authority(),
            ),
        )
