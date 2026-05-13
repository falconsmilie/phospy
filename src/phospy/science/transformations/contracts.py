"""Transformer contracts."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Protocol

import pandas as pd

from phospy.science.transformations.models import IntensityScaleState


def _default_transformation_provenance() -> dict[str, object]:
    return {}


@dataclass(frozen=True, slots=True)
class TransformationResult:
    """Result of applying a transformer to dataset quantitative matrices."""

    phospho: pd.DataFrame
    total: pd.DataFrame | None
    state: IntensityScaleState
    provenance: Mapping[str, object] = field(
        default_factory=_default_transformation_provenance
    )


class Transformer(Protocol):
    """Contract for components that establish dataset intensity scale state."""

    @property
    def preserves_input_scale_state(self) -> bool:
        """Return whether declared input scale state may be preserved."""
        ...

    @property
    def changes_numeric_values(self) -> bool:
        """Return whether transformer run changes matrix numeric values."""
        ...

    @property
    def requires_established_input_state(self) -> bool:
        """Return whether transformer requires explicit established input scale state."""
        ...

    def run(
        self,
        phospho: pd.DataFrame,
        total: pd.DataFrame | None = None,
    ) -> TransformationResult:
        """Apply the supported transformation policy and return explicit state."""
        ...
