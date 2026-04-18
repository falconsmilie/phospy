"""Public result models."""

from __future__ import annotations

from dataclasses import InitVar, dataclass

import pandas as pd

from phospy._frame_ownership import own_optional_dataframe
from phospy.activities.models import KinaseActivityResult
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.validation import WorkflowValidationError
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

    def __post_init__(self) -> None:
        if not isinstance(self.dataset, AnalysisReadyPhosphoDataset):
            raise WorkflowValidationError(
                "kinase workflow result dataset must be AnalysisReadyPhosphoDataset"
            )
        if not isinstance(self.references, ReferenceBundle):
            raise WorkflowValidationError(
                "kinase workflow result references must be ReferenceBundle"
            )
        if not isinstance(self.scoring_result, KinaseScoringResult):
            raise WorkflowValidationError(
                "kinase workflow result scoring_result must be KinaseScoringResult"
            )
        if not isinstance(self.prediction_result, KinasePredictionResult):
            raise WorkflowValidationError(
                "kinase workflow result prediction_result must be KinasePredictionResult"
            )
        if self.activity_result is not None and not isinstance(
            self.activity_result, KinaseActivityResult
        ):
            raise WorkflowValidationError(
                "kinase workflow result activity_result must be KinaseActivityResult or None"
            )


@dataclass(frozen=True, slots=True)
class SignalomeWorkflowResult:
    """Top-level public signalome workflow result."""

    dataset: AnalysisReadyPhosphoDataset
    kinase_result: KinaseWorkflowResult
    module_assignments: SignalomeAssignments
    signalome_modules: SignalomeModules
    kinase_network: KinaseNetwork
    expanded_signalome: pd.DataFrame | None = None
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        if not isinstance(self.dataset, AnalysisReadyPhosphoDataset):
            raise WorkflowValidationError(
                "signalome workflow result dataset must be AnalysisReadyPhosphoDataset"
            )
        if not isinstance(self.kinase_result, KinaseWorkflowResult):
            raise WorkflowValidationError(
                "signalome workflow result kinase_result must be KinaseWorkflowResult"
            )
        if not isinstance(self.module_assignments, SignalomeAssignments):
            raise WorkflowValidationError(
                "signalome workflow result module_assignments must be SignalomeAssignments"
            )
        if not isinstance(self.signalome_modules, SignalomeModules):
            raise WorkflowValidationError(
                "signalome workflow result signalome_modules must be SignalomeModules"
            )
        if not isinstance(self.kinase_network, KinaseNetwork):
            raise WorkflowValidationError(
                "signalome workflow result kinase_network must be KinaseNetwork"
            )
        expanded_signalome = own_optional_dataframe(
            self.expanded_signalome,
            field_name="signalome_result.expanded_signalome",
            error_type=WorkflowValidationError,
            assume_owned=_assume_owned,
        )
        object.__setattr__(self, "expanded_signalome", expanded_signalome)

    @classmethod
    def _from_owned(
        cls,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        kinase_result: KinaseWorkflowResult,
        module_assignments: SignalomeAssignments,
        signalome_modules: SignalomeModules,
        kinase_network: KinaseNetwork,
        expanded_signalome: pd.DataFrame | None = None,
    ) -> SignalomeWorkflowResult:
        return cls(
            dataset=dataset,
            kinase_result=kinase_result,
            module_assignments=module_assignments,
            signalome_modules=signalome_modules,
            kinase_network=kinase_network,
            expanded_signalome=expanded_signalome,
            _assume_owned=True,
        )


__all__ = [
    "KinaseActivityResult",
    "KinasePredictionResult",
    "KinaseScoringResult",
    "KinaseWorkflowResult",
    "SignalomeWorkflowResult",
]
