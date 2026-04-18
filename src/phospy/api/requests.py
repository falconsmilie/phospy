"""Public request models."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from phospy.api.configs import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    SignalomeConfig,
)
from phospy.datasets.builders.contracts import DatasetInput
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.references.models import Organism, ReferenceBundle, ReferencePreset

if TYPE_CHECKING:
    from phospy.api.results import KinaseWorkflowResult


@dataclass(frozen=True, slots=True)
class DatasetBuildRequest:
    """Request for building an ``AnalysisReadyPhosphoDataset``.

    Supported public inputs are pandas ``DataFrame`` values or file paths.
    """

    phospho: DatasetInput
    site_metadata: DatasetInput
    sample_metadata: DatasetInput | None = None
    total: DatasetInput | None = None
    organism: Organism | None = None


@dataclass(frozen=True, slots=True)
class KinaseWorkflowRequest:
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
    """Request for the public signalome workflow.

    Signalome execution requires resolvable protein identity per interpreted site:
    ``dataset.site_metadata.protein_id`` when present, otherwise non-empty protein
    prefixes in the interpreted site identifiers.
    """

    kinase_result: KinaseWorkflowResult
    config: SignalomeConfig = field(default_factory=SignalomeConfig)


__all__ = [
    "DatasetBuildRequest",
    "KinaseWorkflowRequest",
    "SignalomeWorkflowRequest",
]
