"""Identity transformer used for the initial rewrite boundary."""

from __future__ import annotations

import pandas as pd

from phospy.transformations.contracts import TransformationResult
from phospy.transformations.models import TransformationState


class IdentityTransformer:
    """Establish a linear transformation state without changing matrix values."""

    def run(
        self,
        phospho: pd.DataFrame,
        total: pd.DataFrame | None = None,
    ) -> TransformationResult:
        return TransformationResult(
            phospho=phospho,
            total=total,
            state=TransformationState.raw(has_total_matrix=total is not None),
        )
