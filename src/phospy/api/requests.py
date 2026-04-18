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
from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import WorkflowValidationError
from phospy.references.models import Organism, ReferenceBundle, ReferencePreset
from phospy.validation.datasets.inputs import DatasetInputSourceValidator

if TYPE_CHECKING:
    from phospy.api.results import KinaseWorkflowResult

_DATASET_INPUT_VALIDATOR = DatasetInputSourceValidator()


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

    def __post_init__(self) -> None:
        _DATASET_INPUT_VALIDATOR.run(self.phospho, field_name="phospho")
        _DATASET_INPUT_VALIDATOR.run(self.site_metadata, field_name="site_metadata")
        _DATASET_INPUT_VALIDATOR.run(
            self.sample_metadata,
            field_name="sample_metadata",
            allow_none=True,
        )
        _DATASET_INPUT_VALIDATOR.run(
            self.total,
            field_name="total",
            allow_none=True,
        )
        if self.organism is not None and not isinstance(self.organism, Organism):
            raise PhosPyInputError("dataset build request organism must be an Organism")


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
    """Request for the public signalome workflow.

    Signalome execution requires resolvable protein identity per interpreted site:
    ``dataset.site_metadata.protein_id`` when present, otherwise non-empty protein
    prefixes in the interpreted site identifiers.
    """

    kinase_result: KinaseWorkflowResult
    config: SignalomeConfig = field(default_factory=SignalomeConfig)

    def __post_init__(self) -> None:
        from phospy.api.results import KinaseWorkflowResult

        if not isinstance(self.kinase_result, KinaseWorkflowResult):
            raise WorkflowValidationError(
                "signalome workflow request kinase_result must be KinaseWorkflowResult"
            )
        if not isinstance(self.config, SignalomeConfig):
            raise WorkflowValidationError(
                "signalome workflow request config must be SignalomeConfig"
            )


__all__ = [
    "DatasetBuildRequest",
    "KinaseWorkflowRequest",
    "SignalomeWorkflowRequest",
]
