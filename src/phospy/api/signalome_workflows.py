from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

from ..prediction.results import KinasePredictionResult, PredMatResult
from ..prediction.scoring import KinaseScoringResult
from ..signalomes import SignalomeModuleSelectionPolicy
from ..signalomes.analysis import execute_signalome_inputs
from ..signalomes.results import SignalomeResult
from ..validation.requests.signalome import SignalomeInputs, validate_signalome_request

__all__ = ["SignalomeWorkflow"]


class SignalomeWorkflow:
    """Construct signalomes from validated scoring and prediction outputs."""

    def _validate_request(
        self,
        *,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult | PredMatResult,
        expression_matrix: pd.DataFrame,
        kinases_of_interest: Sequence[str],
        site_to_protein: Mapping[str, str] | None = None,
        kinase_network_threshold: float = 0.9,
        signalome_cutoff: float = 0.5,
        module_count: int | None = None,
        min_kinase_module_share_percent: float = 1.0,
        module_selection_policy: SignalomeModuleSelectionPolicy | None = None,
    ) -> SignalomeInputs:
        return validate_signalome_request(
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            expression_matrix=expression_matrix,
            kinases_of_interest=kinases_of_interest,
            site_to_protein=site_to_protein,
            kinase_network_threshold=kinase_network_threshold,
            signalome_cutoff=signalome_cutoff,
            module_count=module_count,
            min_kinase_module_share_percent=min_kinase_module_share_percent,
            module_selection_policy=module_selection_policy,
        )

    def run(
        self,
        *,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult | PredMatResult,
        expression_matrix: pd.DataFrame,
        kinases_of_interest: Sequence[str],
        site_to_protein: Mapping[str, str] | None = None,
        kinase_network_threshold: float = 0.9,
        signalome_cutoff: float = 0.5,
        module_count: int | None = None,
        min_kinase_module_share_percent: float = 1.0,
        module_selection_policy: SignalomeModuleSelectionPolicy | None = None,
    ) -> SignalomeResult:
        request = self._validate_request(
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            expression_matrix=expression_matrix,
            kinases_of_interest=kinases_of_interest,
            site_to_protein=site_to_protein,
            kinase_network_threshold=kinase_network_threshold,
            signalome_cutoff=signalome_cutoff,
            module_count=module_count,
            min_kinase_module_share_percent=min_kinase_module_share_percent,
            module_selection_policy=module_selection_policy,
        )
        return self.run_validated(request)

    def run_validated(
        self,
        request: SignalomeInputs,
    ) -> SignalomeResult:
        return execute_signalome_inputs(request)
