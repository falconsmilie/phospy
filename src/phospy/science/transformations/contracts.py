"""Transformer contracts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import pandas as pd

from phospy.science.transformations.models import IntensityScaleState


@dataclass(frozen=True, slots=True)
class TransformationResult:
    """Result of applying a transformer to dataset quantitative matrices."""

    phospho: pd.DataFrame
    total: pd.DataFrame | None
    state: IntensityScaleState


class Transformer(Protocol):
    """Contract for components that establish dataset intensity scale state."""

    def run(
        self,
        phospho: pd.DataFrame,
        total: pd.DataFrame | None = None,
    ) -> TransformationResult:
        """Apply the supported transformation policy and return explicit state."""
        ...
