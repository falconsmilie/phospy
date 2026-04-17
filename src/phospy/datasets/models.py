"""Dataset domain models."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from phospy.references.models import Organism
from phospy.transformations.models import TransformationState


@dataclass(frozen=True, slots=True)
class AnalysisReadyPhosphoDataset:
    """Public analysis-ready dataset contract."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None = None
    total: pd.DataFrame | None = None
    organism: Organism | None = None
    transformation_state: TransformationState = field(
        default_factory=TransformationState
    )
