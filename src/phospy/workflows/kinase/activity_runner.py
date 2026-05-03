"""Activity orchestration for kinase workflow execution."""

from __future__ import annotations

from phospy.activities.methods import (
    KSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG,
    KseaZScoreActivityMethod,
    SimplifiedWeightedSubstrateActivityMethod,
)
from phospy.activities.models import KinaseActivityResult
from phospy.api.configs import (
    KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
    KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
)
from phospy.errors.workflows import WorkflowBoundaryError
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
    ) -> None:
        self._activity_input_validator = (
            activity_input_validator or KinaseActivityInputValidator()
        )

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
        if (
            activity_config.method
            == KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY
        ):
            validated_inputs = self._activity_input_validator.run(
                pred_mat=prediction_result.pred_mat,
                phospho_matrix=request.activity_phospho_matrix,
                threshold=activity_config.threshold,
                min_substrates=activity_config.min_substrates,
                top_n_substrates=activity_config.top_n_substrates,
            )
            return SimplifiedWeightedSubstrateActivityMethod(
                threshold=float(activity_config.threshold),
                min_substrates=int(activity_config.min_substrates),
                top_n_substrates=int(activity_config.top_n_substrates),
            ).run(validated_inputs)
        if activity_config.method == KINASE_ACTIVITY_METHOD_KSEA_ZSCORE:
            validated_inputs = self._activity_input_validator.run(
                pred_mat=prediction_result.pred_mat,
                phospho_matrix=request.activity_phospho_matrix,
                threshold=activity_config.ksea_evidence_threshold,
                min_substrates=activity_config.ksea_min_substrates,
                top_n_substrates=1,
            )
            return KseaZScoreActivityMethod(
                evidence_threshold=float(activity_config.ksea_evidence_threshold),
                min_substrates=int(activity_config.ksea_min_substrates),
                p_value_method=str(activity_config.ksea_p_value_method),
                adjust_p_values=bool(activity_config.ksea_adjust_p_values),
                q_value_method=KSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG,
            ).run(validated_inputs)
        raise WorkflowBoundaryError(
            seam="kinase.activity.method",
            next_action="select a supported activity method in activity_config.method",
            details={"method": str(activity_config.method)},
            message_prefix="kinase workflow boundary validation failed",
        )


__all__ = ["KinaseActivityRunner"]
