from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

from ..api.contracts import PredictionRunConfig
from ..internal.defaults import DEFAULT_MOTIF_FLANK_SIZE
from ..internal.types import PREDICTION_SVM_MODE_DEFAULT, PredictionSvmMode
from ..prediction.engines import (
    KinaseWorkflowExecutionResult,
    KinaseWorkflowExecutor,
)
from ..references import ReferenceBundle
from ..validation.requests.workflow import WorkflowInputs

__all__ = ["KinaseWorkflow"]


def _validate_workflow_inputs(
    *,
    executor: KinaseWorkflowExecutor,
    phospho_matrix: pd.DataFrame,
    substrate_map: Mapping[str, Sequence[str]] | None = None,
    site_sequences: Mapping[str, str] | pd.Series | None = None,
    motif_sequences: Mapping[str, Sequence[str]] | None = None,
    reference_bundle: ReferenceBundle | None = None,
    prediction_config: PredictionRunConfig | None = None,
) -> WorkflowInputs:
    resolved_config = PredictionRunConfig.from_value(prediction_config)
    return executor.validate_request(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        site_sequences=site_sequences,
        motif_sequences=motif_sequences,
        reference_bundle=reference_bundle,
        min_substrates=resolved_config.min_substrates,
        min_motif_size=resolved_config.min_motif_size,
        allow_profile_only_fallback=resolved_config.allow_profile_only_fallback,
        ensemble_size=resolved_config.ensemble_size,
        top=resolved_config.top,
        score_threshold=resolved_config.score_threshold,
        inclusion=resolved_config.inclusion,
        n_iterations=resolved_config.n_iterations,
        random_state=resolved_config.random_state,
        svm_mode=resolved_config.svm_mode,
        profile_policy=resolved_config.profile_policy,
    )


class KinaseWorkflow:
    """Run the native kinase scoring and prediction workflow end to end."""

    def __init__(
        self,
        flank_size: int = DEFAULT_MOTIF_FLANK_SIZE,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = PREDICTION_SVM_MODE_DEFAULT,
    ) -> None:
        self._executor = KinaseWorkflowExecutor(
            flank_size=flank_size,
            kernel=kernel,
            svm_mode=svm_mode,
        )

    def _validate_request(
        self,
        *,
        phospho_matrix: pd.DataFrame,
        substrate_map: Mapping[str, Sequence[str]] | None = None,
        site_sequences: Mapping[str, str] | pd.Series | None = None,
        motif_sequences: Mapping[str, Sequence[str]] | None = None,
        reference_bundle: ReferenceBundle | None = None,
        prediction_config: PredictionRunConfig | None = None,
    ) -> WorkflowInputs:
        return _validate_workflow_inputs(
            executor=self._executor,
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            motif_sequences=motif_sequences,
            reference_bundle=reference_bundle,
            prediction_config=prediction_config,
        )

    def run(
        self,
        phospho_matrix: pd.DataFrame,
        substrate_map: Mapping[str, Sequence[str]] | None = None,
        site_sequences: Mapping[str, str] | pd.Series | None = None,
        motif_sequences: Mapping[str, Sequence[str]] | None = None,
        reference_bundle: ReferenceBundle | None = None,
        prediction_config: PredictionRunConfig | None = None,
    ) -> KinaseWorkflowExecutionResult:
        request = self._validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            motif_sequences=motif_sequences,
            reference_bundle=reference_bundle,
            prediction_config=prediction_config,
        )
        return self.run_validated(request)

    def run_validated(self, request: WorkflowInputs) -> KinaseWorkflowExecutionResult:
        return self._executor.execute_validated_request(request)
