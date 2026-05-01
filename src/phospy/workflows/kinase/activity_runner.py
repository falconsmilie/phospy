"""Activity orchestration for kinase workflow execution."""

from __future__ import annotations

from collections.abc import Callable

from phospy.activities.models import KinaseActivityResult
from phospy.activities.scoring import compute_activity_from_inputs
from phospy.prediction.models import KinasePredictionResult
from phospy.validation.workflows.activity import KinaseActivityInputValidator
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)


class KinaseActivityRunner:
    """Run optional activity-stage execution from prediction outputs."""

    def __init__(
        self,
        *,
        activity_input_validator: KinaseActivityInputValidator | None = None,
        compute_activity: Callable[..., KinaseActivityResult] = (
            compute_activity_from_inputs
        ),
    ) -> None:
        self._activity_input_validator = (
            activity_input_validator or KinaseActivityInputValidator()
        )
        self._compute_activity = compute_activity

    def run(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        config: ResolvedKinaseExecutionConfig,
        prediction_result: KinasePredictionResult,
    ) -> KinaseActivityResult | None:
        activity_config = config.activity
        if activity_config is None:
            return None
        validated_inputs = self._activity_input_validator.run(
            pred_mat=prediction_result.pred_mat,
            phospho_matrix=request.activity_phospho_matrix,
            threshold=activity_config.threshold,
            min_substrates=activity_config.min_substrates,
            top_n_substrates=activity_config.top_n_substrates,
        )
        return self._compute_activity(validated_inputs)


__all__ = ["KinaseActivityRunner"]
