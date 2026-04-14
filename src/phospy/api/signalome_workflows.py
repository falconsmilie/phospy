from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

from ..prediction.results import KinasePredictionResult, PredMatResult
from ..prediction.scoring import KinaseScoringResult
from ..signalomes.analysis import execute_signalome_inputs
from ..signalomes.results import SignalomeResult
from ..validation.requests.signalome import SignalomeInputs, validate_signalome_request
from .contracts import SignalomeRunConfig

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
        config: SignalomeRunConfig | None = None,
    ) -> SignalomeInputs:
        resolved_config = SignalomeRunConfig.from_value(config)
        return validate_signalome_request(
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            expression_matrix=expression_matrix,
            kinases_of_interest=kinases_of_interest,
            site_to_protein=site_to_protein,
            kinase_network_threshold=resolved_config.kinase_network_threshold,
            kinase_network_policy=resolved_config.kinase_network_policy,
            signalome_cutoff=resolved_config.signalome_cutoff,
            module_count=resolved_config.module_count,
            min_kinase_module_share_percent=(
                resolved_config.min_kinase_module_share_percent
            ),
            module_selection_policy=resolved_config.module_selection_policy,
        )

    def run(
        self,
        *,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult | PredMatResult,
        expression_matrix: pd.DataFrame,
        kinases_of_interest: Sequence[str],
        site_to_protein: Mapping[str, str] | None = None,
        config: SignalomeRunConfig | None = None,
    ) -> SignalomeResult:
        request = self._validate_request(
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            expression_matrix=expression_matrix,
            kinases_of_interest=kinases_of_interest,
            site_to_protein=site_to_protein,
            config=config,
        )
        return self.run_validated(request)

    def run_validated(
        self,
        request: SignalomeInputs,
    ) -> SignalomeResult:
        return execute_signalome_inputs(request)
