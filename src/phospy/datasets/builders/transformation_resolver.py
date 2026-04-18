"""Internal transformation-state establishment for dataset builder execution."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.transformations.contracts import Transformer
from phospy.transformations.models import TransformationState
from phospy.transformations.transformers.identity import IdentityTransformer


@dataclass(frozen=True, slots=True)
class ResolvedTransformation:
    """Quantitative matrices paired with an established transformation state."""

    phospho: pd.DataFrame
    total: pd.DataFrame | None
    transformation_state: TransformationState


class DatasetTransformationResolver:
    """Resolve transformation state for a dataset build request."""

    def __init__(self, *, transformer: Transformer | None = None) -> None:
        self._transformer = transformer or IdentityTransformer()

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        total: pd.DataFrame | None,
        transformation_state: TransformationState | None,
    ) -> ResolvedTransformation:
        if transformation_state is not None:
            return ResolvedTransformation(
                phospho=phospho,
                total=total,
                transformation_state=transformation_state,
            )
        transformed = self._transformer.run(phospho=phospho, total=total)
        return ResolvedTransformation(
            phospho=transformed.phospho,
            total=transformed.total,
            transformation_state=transformed.state,
        )
