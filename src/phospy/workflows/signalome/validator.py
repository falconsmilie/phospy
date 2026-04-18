"""Internal validator for signalome workflow requests."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.api.requests import SignalomeWorkflowRequest
from phospy.api.results import KinaseWorkflowResult
from phospy.errors.validation import WorkflowValidationError
from phospy.validation.common.dataframes import (
    require_dataframe,
    require_numeric_dataframe,
)
from phospy.validation.workflows.configs import WorkflowConfigValidator


class SignalomeWorkflowValidator:
    """Validate `SignalomeWorkflowRequest` before interpretation."""

    def __init__(
        self, *, config_validator: WorkflowConfigValidator | None = None
    ) -> None:
        self._config_validator = config_validator or WorkflowConfigValidator()

    def run(self, request: SignalomeWorkflowRequest) -> SignalomeWorkflowRequest:
        if not isinstance(request, SignalomeWorkflowRequest):
            raise WorkflowValidationError(
                "signalome workflow input must be a SignalomeWorkflowRequest"
            )
        if not isinstance(request.kinase_result, KinaseWorkflowResult):
            raise WorkflowValidationError(
                "signalome workflow request kinase_result must be KinaseWorkflowResult"
            )
        self._config_validator.run_signalome(request.config)

        prediction_matrix = require_dataframe(
            request.kinase_result.prediction_result.pred_mat,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            allow_empty=False,
            error_type=WorkflowValidationError,
        )
        require_numeric_dataframe(
            prediction_matrix,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            error_type=WorkflowValidationError,
        )
        self._require_no_missing_or_infinite(
            prediction_matrix,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
        )
        if prediction_matrix.shape[1] == 0:
            raise WorkflowValidationError(
                "signalome workflow request kinase_result.prediction_result.pred_mat "
                "must contain at least one kinase column"
            )

        scoring_result = request.kinase_result.scoring_result
        score_matrix = require_dataframe(
            scoring_result.profile_scores,
            field_name=(
                "signalome workflow request kinase_result.scoring_result.profile_scores"
            ),
            allow_empty=False,
            error_type=WorkflowValidationError,
        )
        require_numeric_dataframe(
            score_matrix,
            field_name=(
                "signalome workflow request kinase_result.scoring_result.profile_scores"
            ),
            error_type=WorkflowValidationError,
        )
        self._require_no_missing_or_infinite(
            score_matrix,
            field_name=(
                "signalome workflow request kinase_result.scoring_result.profile_scores"
            ),
        )
        if score_matrix.shape[1] == 0:
            raise WorkflowValidationError(
                "signalome workflow request kinase_result.scoring_result.profile_scores "
                "must contain at least one kinase column"
            )
        return request

    @staticmethod
    def _require_no_missing_or_infinite(
        frame: pd.DataFrame, *, field_name: str
    ) -> None:
        if frame.isna().to_numpy().any():
            raise WorkflowValidationError(
                f"{field_name} must not contain missing values"
            )
        if np.isinf(frame.to_numpy(copy=False)).any():
            raise WorkflowValidationError(
                f"{field_name} must contain finite numeric values"
            )
