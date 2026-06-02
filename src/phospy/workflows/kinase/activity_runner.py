"""Activity orchestration for kinase workflow execution."""

from __future__ import annotations

import pandas as pd

from phospy.contracts.configs import (
    KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
    KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
)
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.activities.methods import (
    KSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG,
    KseaZScoreActivityMethod,
    SimplifiedWeightedSubstrateActivityMethod,
)
from phospy.science.activities.models import KinaseActivityResult
from phospy.science.prediction.models import KinasePredictionResult
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
        site_identity_map = _require_site_identity_map(request.site_identity_map)
        if (
            activity_config.method
            == KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY
        ):
            # Prediction matrix kinase columns are expected to already be
            # normalized upstream by the reference-ingestion boundary.
            validated_inputs = self._activity_input_validator.run(
                pred_mat=prediction_result._borrow_pred_mat_frame(),
                phospho_matrix=request.activity_phospho_matrix,
                threshold=activity_config.threshold,
                min_substrates=activity_config.min_substrates,
                top_n_substrates=activity_config.top_n_substrates,
            )
            result = SimplifiedWeightedSubstrateActivityMethod(
                threshold=float(activity_config.threshold),
                min_substrates=int(activity_config.min_substrates),
                top_n_substrates=int(activity_config.top_n_substrates),
            ).run(validated_inputs)
            return self._annotate_activity_result(
                activity_result=result,
                site_identity_map=site_identity_map,
            )
        if activity_config.method == KINASE_ACTIVITY_METHOD_KSEA_ZSCORE:
            validated_inputs = self._activity_input_validator.run(
                pred_mat=prediction_result._borrow_pred_mat_frame(),
                phospho_matrix=request.activity_phospho_matrix,
                threshold=activity_config.ksea_evidence_threshold,
                min_substrates=activity_config.ksea_min_substrates,
                top_n_substrates=1,
            )
            result = KseaZScoreActivityMethod(
                evidence_threshold=float(activity_config.ksea_evidence_threshold),
                min_substrates=int(activity_config.ksea_min_substrates),
                p_value_method=str(activity_config.ksea_p_value_method),
                adjust_p_values=bool(activity_config.ksea_adjust_p_values),
                q_value_method=KSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG,
            ).run(validated_inputs)
            return self._annotate_activity_result(
                activity_result=result,
                site_identity_map=site_identity_map,
            )
        raise WorkflowBoundaryError(
            seam="kinase.activity.method",
            next_action="select a supported activity method in activity_config.method",
            details={"method": str(activity_config.method)},
            message_prefix="kinase workflow boundary validation failed",
        )

    @staticmethod
    def _annotate_activity_result(
        *,
        activity_result: KinaseActivityResult,
        site_identity_map: pd.DataFrame,
    ) -> KinaseActivityResult:
        target_table = activity_result.target_table.copy(deep=True)
        if target_table.empty:
            target_table.loc[:, "site_key"] = pd.Series(dtype="object")
            target_table.loc[:, "display_id"] = pd.Series(dtype="object")
        else:
            display_lookup = {
                str(site_key): str(display_id)
                for site_key, display_id in site_identity_map.loc[
                    :, ["site_key", "display_id"]
                ].itertuples(index=False)
            }
            site_keys = target_table.loc[:, "site_id"].astype(str)
            display_ids = site_keys.map(lambda value: display_lookup.get(value, value))
            target_table.loc[:, "site_key"] = site_keys
            target_table.loc[:, "display_id"] = display_ids
            target_table.loc[:, "site_id"] = display_ids
        return KinaseActivityResult(
            weighted_activity=activity_result.weighted_activity,
            thresholded_substrate_mean_activity=(
                activity_result.thresholded_substrate_mean_activity
            ),
            thresholded_substrate_counts=activity_result.thresholded_substrate_counts,
            target_counts=activity_result.target_counts,
            target_table=target_table,
            threshold_membership_diagnostics=(
                activity_result.threshold_membership_diagnostics
            ),
            activity_substrate_counts=activity_result.activity_substrate_counts,
            statistics_table=activity_result.statistics_table,
            method_summary=activity_result.method_summary,
            activity_method=activity_result.activity_method,
        )


def _require_site_identity_map(site_identity_map: pd.DataFrame | None) -> pd.DataFrame:
    if site_identity_map is None:
        raise WorkflowBoundaryError(
            seam="kinase.activity.site_identity_map",
            next_action="ensure kinase workflow interpretation resolves site identity mapping",
            message_prefix="kinase workflow boundary validation failed",
        )
    return site_identity_map


__all__ = ["KinaseActivityRunner"]
