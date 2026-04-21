"""Public request models."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from phospy.api.configs import (
    DatasetPreprocessingConfig,
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
    Preprocessing policy remains builder-owned via ``preprocessing_config``,
    including explicit site-matrix duplicate-policy controls under a strict
    complete-case analysis-ready boundary.
    """

    phospho: DatasetInput
    site_metadata: DatasetInput
    sample_metadata: DatasetInput | None = None
    total: DatasetInput | None = None
    organism: Organism | None = None
    preprocessing_config: DatasetPreprocessingConfig = field(
        default_factory=DatasetPreprocessingConfig
    )


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

    Signalome execution requires explicit protein identity per interpreted site via
    ``dataset.site_metadata.protein_id``.
    """

    kinase_result: KinaseWorkflowResult
    config: SignalomeConfig = field(default_factory=SignalomeConfig)


__all__ = [
    "DatasetBuildRequest",
    "KinaseWorkflowRequest",
    "SignalomeWorkflowRequest",
]
