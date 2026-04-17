"""Internal validator for signalome workflow requests."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.api.configs import SignalomeConfig
from phospy.api.requests import SignalomeWorkflowRequest
from phospy.api.results import SimpleKinaseWorkflowResult
from phospy.errors.validation import WorkflowValidationError


class SignalomeWorkflowValidator:
    """Validate `SignalomeWorkflowRequest` before interpretation."""

    _MIN_CUTOFF = 0.0
    _MAX_CUTOFF = 1.0

    def run(self, request: SignalomeWorkflowRequest) -> SignalomeWorkflowRequest:
        if not isinstance(request, SignalomeWorkflowRequest):
            raise WorkflowValidationError(
                "signalome workflow input must be a SignalomeWorkflowRequest"
            )
        if not isinstance(request.kinase_result, SimpleKinaseWorkflowResult):
            raise WorkflowValidationError(
                "signalome workflow request kinase_result must be SimpleKinaseWorkflowResult"
            )
        if not isinstance(request.config, SignalomeConfig):
            raise WorkflowValidationError(
                "signalome workflow request config must be SignalomeConfig"
            )
        cutoff = request.config.signalome_cutoff
        if isinstance(cutoff, bool) or not isinstance(cutoff, (int, float)):
            raise WorkflowValidationError(
                "signalome workflow request config.signalome_cutoff must be a float "
                "between 0.0 and 1.0"
            )
        if not self._MIN_CUTOFF <= float(cutoff) <= self._MAX_CUTOFF:
            raise WorkflowValidationError(
                "signalome workflow request config.signalome_cutoff must be between "
                "0.0 and 1.0"
            )

        prediction_matrix = self._validated_numeric_matrix(
            request.kinase_result.prediction_result.pred_mat,
            context="kinase_result.prediction_result.pred_mat",
        )
        if prediction_matrix.shape[1] == 0:
            raise WorkflowValidationError(
                "signalome workflow request kinase_result.prediction_result.pred_mat "
                "must contain at least one kinase column"
            )

        scoring_result = request.kinase_result.scoring_result
        score_matrix_source = scoring_result.combined_scores
        if score_matrix_source is None:
            score_matrix_source = scoring_result.profile_scores
        score_matrix = self._validated_numeric_matrix(
            score_matrix_source,
            context="kinase_result.scoring_result.combined_scores",
        )
        if score_matrix.shape[1] == 0:
            raise WorkflowValidationError(
                "signalome workflow request kinase_result.scoring_result.combined_scores "
                "must contain at least one kinase column"
            )
        return request

    @staticmethod
    def _validated_numeric_matrix(matrix: object, *, context: str) -> pd.DataFrame:
        if not isinstance(matrix, pd.DataFrame):
            raise WorkflowValidationError(
                f"signalome workflow request {context} must be a pandas DataFrame"
            )
        if matrix.empty:
            raise WorkflowValidationError(
                f"signalome workflow request {context} must not be empty"
            )
        try:
            numeric_matrix = matrix.astype(float)
        except (TypeError, ValueError) as exc:
            raise WorkflowValidationError(
                f"signalome workflow request {context} must contain numeric values"
            ) from exc
        if not np.isfinite(numeric_matrix.to_numpy(dtype=float, copy=False)).all():
            raise WorkflowValidationError(
                f"signalome workflow request {context} must contain finite numeric values"
            )
        return numeric_matrix
