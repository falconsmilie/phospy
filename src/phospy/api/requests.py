"""Public request models."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from phospy.api.configs import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    SignalomeConfig,
)
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.references.models import Organism, ReferenceBundle, ReferencePreset
from phospy.transformations.models import TransformationState

if TYPE_CHECKING:
    from phospy.api.results import SimpleKinaseWorkflowResult


@dataclass(frozen=True, slots=True)
class DatasetBuildRequest:
    """Request for building an ``AnalysisReadyPhosphoDataset``."""

    phospho: pd.DataFrame | str | Path
    site_metadata: pd.DataFrame | str | Path | None = None
    sample_metadata: pd.DataFrame | str | Path | None = None
    total: pd.DataFrame | str | Path | None = None
    organism: Organism | None = None
    transformation_state: TransformationState | None = None


@dataclass(frozen=True, slots=True)
class SimpleKinaseWorkflowRequest:
    """Request for the public kinase workflow."""

    dataset: AnalysisReadyPhosphoDataset
    references: ReferencePreset | ReferenceBundle = ReferencePreset.AUTO
    scoring_config: KinaseScoringConfig = field(default_factory=KinaseScoringConfig)
    prediction_config: KinasePredictionConfig = field(
        default_factory=KinasePredictionConfig
    )
    activity_config: KinaseActivityConfig | None = field(
        default_factory=KinaseActivityConfig
    )


@dataclass(frozen=True, slots=True)
class SignalomeWorkflowRequest:
    """Request for the public signalome workflow."""

    kinase_result: SimpleKinaseWorkflowResult
    config: SignalomeConfig = field(default_factory=SignalomeConfig)


__all__ = [
    "DatasetBuildRequest",
    "SignalomeWorkflowRequest",
    "SimpleKinaseWorkflowRequest",
]
