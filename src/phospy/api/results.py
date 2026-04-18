"""Public result models."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.activities.models import KinaseActivityResult
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.references.models import ReferenceBundle
from phospy.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
)


@dataclass(frozen=True, slots=True)
class KinaseWorkflowResult:
    """Top-level public kinase workflow result."""

    dataset: AnalysisReadyPhosphoDataset
    references: ReferenceBundle
    scoring_result: KinaseScoringResult
    prediction_result: KinasePredictionResult
    activity_result: KinaseActivityResult | None = None


@dataclass(frozen=True, slots=True)
class SignalomeWorkflowResult:
    """Top-level public signalome workflow result."""

    dataset: AnalysisReadyPhosphoDataset
    kinase_result: KinaseWorkflowResult
    module_assignments: SignalomeAssignments
    signalome_modules: SignalomeModules
    kinase_network: KinaseNetwork
    expanded_signalome: pd.DataFrame | None = None


__all__ = [
    "KinaseActivityResult",
    "KinasePredictionResult",
    "KinaseScoringResult",
    "KinaseWorkflowResult",
    "SignalomeWorkflowResult",
]
