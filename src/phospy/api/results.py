"""Public result models."""

from __future__ import annotations

from dataclasses import InitVar, dataclass, field

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
    SignalomeModuleSelectionDiagnostics,
    default_signalome_module_selection_diagnostics,
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
    module_selection_diagnostics: SignalomeModuleSelectionDiagnostics = field(
        default_factory=default_signalome_module_selection_diagnostics
    )
    expanded_signalome: pd.DataFrame | None = None
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
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
        module_selection_diagnostics: SignalomeModuleSelectionDiagnostics | None = None,
        expanded_signalome: pd.DataFrame | None = None,
    ) -> SignalomeWorkflowResult:
        return cls(
            dataset=dataset,
            kinase_result=kinase_result,
            module_assignments=module_assignments,
            signalome_modules=signalome_modules,
            kinase_network=kinase_network,
            module_selection_diagnostics=(
                default_signalome_module_selection_diagnostics()
                if module_selection_diagnostics is None
                else module_selection_diagnostics
            ),
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
