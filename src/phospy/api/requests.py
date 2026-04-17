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
from phospy.datasets.builders.contracts import DatasetInput
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import WorkflowValidationError
from phospy.references.models import Organism, ReferenceBundle, ReferencePreset
from phospy.transformations.models import TransformationState

if TYPE_CHECKING:
    from phospy.api.results import SimpleKinaseWorkflowResult


@dataclass(frozen=True, slots=True)
class DatasetBuildRequest:
    """Request for building an ``AnalysisReadyPhosphoDataset``."""

    phospho: DatasetInput
    site_metadata: DatasetInput | None = None
    sample_metadata: DatasetInput | None = None
    total: DatasetInput | None = None
    organism: Organism | None = None
    transformation_state: TransformationState | None = None

    def __post_init__(self) -> None:
        self._validate_input(self.phospho, field_name="phospho")
        self._validate_optional_input(self.site_metadata, field_name="site_metadata")
        self._validate_optional_input(
            self.sample_metadata,
            field_name="sample_metadata",
        )
        self._validate_optional_input(self.total, field_name="total")
        if self.organism is not None and not isinstance(self.organism, Organism):
            raise PhosPyInputError("dataset build request organism must be an Organism")
        if self.transformation_state is not None and not isinstance(
            self.transformation_state, TransformationState
        ):
            raise PhosPyInputError(
                "dataset build request transformation_state must be a TransformationState"
            )

    @staticmethod
    def _validate_input(value: object, *, field_name: str) -> None:
        if isinstance(value, (pd.DataFrame, str, Path)):
            return
        raise PhosPyInputError(
            f"dataset build request {field_name} must be a DataFrame, str, or Path"
        )

    @classmethod
    def _validate_optional_input(cls, value: object | None, *, field_name: str) -> None:
        if value is None:
            return
        cls._validate_input(value, field_name=field_name)


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

    def __post_init__(self) -> None:
        if not isinstance(self.dataset, AnalysisReadyPhosphoDataset):
            raise WorkflowValidationError(
                "kinase workflow request dataset must be AnalysisReadyPhosphoDataset"
            )
        if not isinstance(self.references, (ReferencePreset, ReferenceBundle)):
            raise WorkflowValidationError(
                "kinase workflow request references must be ReferencePreset or ReferenceBundle"
            )
        if not isinstance(self.scoring_config, KinaseScoringConfig):
            raise WorkflowValidationError(
                "kinase workflow request scoring_config must be KinaseScoringConfig"
            )
        if not isinstance(self.prediction_config, KinasePredictionConfig):
            raise WorkflowValidationError(
                "kinase workflow request prediction_config must be KinasePredictionConfig"
            )
        if self.activity_config is not None and not isinstance(
            self.activity_config, KinaseActivityConfig
        ):
            raise WorkflowValidationError(
                "kinase workflow request activity_config must be KinaseActivityConfig or None"
            )


@dataclass(frozen=True, slots=True)
class SignalomeWorkflowRequest:
    """Request for the public signalome workflow."""

    kinase_result: SimpleKinaseWorkflowResult
    config: SignalomeConfig = field(default_factory=SignalomeConfig)

    def __post_init__(self) -> None:
        from phospy.api.results import SimpleKinaseWorkflowResult

        if not isinstance(self.kinase_result, SimpleKinaseWorkflowResult):
            raise WorkflowValidationError(
                "signalome workflow request kinase_result must be SimpleKinaseWorkflowResult"
            )
        if not isinstance(self.config, SignalomeConfig):
            raise WorkflowValidationError(
                "signalome workflow request config must be SignalomeConfig"
            )


__all__ = [
    "DatasetBuildRequest",
    "SignalomeWorkflowRequest",
    "SimpleKinaseWorkflowRequest",
]
