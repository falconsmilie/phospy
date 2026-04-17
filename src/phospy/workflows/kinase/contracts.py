"""Kinase workflow internal contracts."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.api.configs import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
)
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.references.models import ReferenceBundle


@dataclass(frozen=True, slots=True)
class ResolvedKinaseWorkflowRequest:
    """Interpreter output for kinase workflow execution."""

    dataset: AnalysisReadyPhosphoDataset
    references: ReferenceBundle
    scoring_config: KinaseScoringConfig
    prediction_config: KinasePredictionConfig
    activity_config: KinaseActivityConfig | None
