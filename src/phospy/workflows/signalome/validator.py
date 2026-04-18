"""Internal validator for signalome workflow requests."""

from __future__ import annotations

from phospy.api.requests import SignalomeWorkflowRequest
from phospy.api.results import KinaseWorkflowResult
from phospy.errors.validation import WorkflowValidationError
from phospy.validation.common.missing_values import MissingValuePolicy
from phospy.validation.common.numeric_frames import require_numeric_matrix
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

        prediction_matrix = require_numeric_matrix(
            request.kinase_result.prediction_result.pred_mat,
            field_name=(
                "signalome workflow request kinase_result.prediction_result.pred_mat"
            ),
            allow_empty=False,
            missing_value_policy=MissingValuePolicy.FORBID,
            error_type=WorkflowValidationError,
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
        score_matrix = require_numeric_matrix(
            score_matrix_source,
            field_name=(
                "signalome workflow request "
                "kinase_result.scoring_result.combined_scores"
            ),
            allow_empty=False,
            missing_value_policy=MissingValuePolicy.FORBID,
            error_type=WorkflowValidationError,
        )
        if score_matrix.shape[1] == 0:
            raise WorkflowValidationError(
                "signalome workflow request kinase_result.scoring_result.combined_scores "
                "must contain at least one kinase column"
            )
        return request
