"""Activity orchestration for kinase workflow execution."""

from __future__ import annotations

import pandas as pd

from phospy.contracts.configs import (
    KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
    KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
    KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
)
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.activities.method_contracts import (
    kinase_activity_method_universe_contract,
)
from phospy.science.activities.methods import (
    KSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG,
    KseaZScoreActivityMethod,
    SimplifiedWeightedSubstrateActivityMethod,
    SsgseaSubstrateEnrichmentActivityMethod,
)
from phospy.science.activities.models import KinaseActivityResult
from phospy.science.activities.semantics import ActivityInputMatrix
from phospy.science.prediction.internal_view import KinasePredictionInternalView
from phospy.science.prediction.models import KinasePredictionResult
from phospy.science.transformations.models import QuantitativeMeaning
from phospy.validation.workflows.activity import KinaseActivityInputValidator
from phospy.workflows.kinase.component_models import KinaseScoringRunResult
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.membership import build_ksea_membership_selection


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
        scoring_execution: KinaseScoringRunResult | None = None,
    ) -> KinaseActivityResult | None:
        activity_config = config.activity
        if activity_config is None:
            return None
        site_identity_map = _require_site_identity_map(request.site_identity_map)
        prediction_view = KinasePredictionInternalView(prediction_result)
        universe_contract = kinase_activity_method_universe_contract(
            activity_config.method
        )
        _validate_prediction_membership_universe(
            pred_mat=prediction_view.pred_mat,
            request=request,
            method_id=universe_contract.method_id,
        )
        if (
            activity_config.method
            == KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY
        ):
            activity_matrix = request.activity_phospho_matrix
            membership_pred_mat = _prediction_membership_matrix_for_site_universe(
                pred_mat=prediction_view.pred_mat,
                site_index=activity_matrix.index,
            )
            # Prediction matrix kinase columns are expected to already be
            # normalized upstream by the reference-ingestion boundary.
            validated_inputs = self._activity_input_validator.run(
                pred_mat=membership_pred_mat,
                phospho_matrix=activity_matrix,
                threshold=activity_config.threshold,
                min_substrates=activity_config.min_substrates,
                top_n_substrates=activity_config.top_n_substrates,
                activity_input=_activity_input_from_dataset(
                    request=request,
                    method=activity_config.method,
                    matrix=activity_matrix,
                ),
                min_fraction=0.0,
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
            activity_matrix = request.ksea_background_phospho_matrix
            membership_pred_mat = _prediction_membership_matrix_for_site_universe(
                pred_mat=prediction_view.pred_mat,
                site_index=activity_matrix.index,
            )
            membership_selection = (
                None
                if scoring_execution is None
                else build_ksea_membership_selection(
                    request=request,
                    config=config,
                    scoring_execution=scoring_execution,
                    prediction_result=prediction_result,
                    evidence_threshold=float(activity_config.ksea_evidence_threshold),
                    membership_matrix=membership_pred_mat,
                )
            )
            validated_inputs = self._activity_input_validator.run(
                pred_mat=membership_pred_mat,
                phospho_matrix=activity_matrix,
                threshold=activity_config.ksea_evidence_threshold,
                min_substrates=activity_config.ksea_min_substrates,
                top_n_substrates=1,
                activity_input=_activity_input_from_dataset(
                    request=request,
                    method=activity_config.method,
                    matrix=activity_matrix,
                ),
                membership_selection=membership_selection,
                min_fraction=0.0,
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
        if activity_config.method == KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT:
            activity_matrix = request.ssgsea_effect_matrix
            membership = _kinase_substrate_membership_for_site_universe(
                membership=request.reference_membership_map,
                site_index=activity_matrix.index,
            )
            result = SsgseaSubstrateEnrichmentActivityMethod(
                min_substrates=int(activity_config.ssgsea_min_substrates),
                ranking_direction=str(activity_config.ssgsea_ranking_direction),
                permutation_count=int(activity_config.ssgsea_permutations),
                random_seed=activity_config.ssgsea_random_seed,
                adjust_p_values=bool(activity_config.ssgsea_adjust_p_values),
            ).run(
                activity_input=_activity_input_from_dataset(
                    request=request,
                    method=activity_config.method,
                    matrix=activity_matrix,
                ),
                kinase_substrate_membership=membership,
            )
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
        return KinaseActivityResult._from_owned(
            activity_matrix=activity_result.activity_matrix,
            p_value_matrix=activity_result.p_value_matrix,
            q_value_matrix=activity_result.q_value_matrix,
            confidence_interval_low=activity_result.confidence_interval_low,
            confidence_interval_high=activity_result.confidence_interval_high,
            substrate_count_matrix=activity_result.substrate_count_matrix,
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
            method_diagnostics=activity_result.method_diagnostics,
            policy_provenance=activity_result.policy_provenance,
            activity_method=activity_result.activity_method,
            input_semantics=activity_result.input_semantics,
            profile_metadata=activity_result.profile_metadata,
            membership_selection=activity_result.membership_selection,
        )


def _require_site_identity_map(site_identity_map: pd.DataFrame | None) -> pd.DataFrame:
    if site_identity_map is None:
        raise WorkflowBoundaryError(
            seam="kinase.activity.site_identity_map",
            next_action="ensure kinase workflow interpretation resolves site identity mapping",
            message_prefix="kinase workflow boundary validation failed",
        )
    return site_identity_map


def _activity_input_from_dataset(
    *,
    request: ResolvedKinaseWorkflowRequest,
    method: str,
    matrix: pd.DataFrame,
) -> ActivityInputMatrix:
    quantity = request.dataset.intensity_scale_state.quantity
    if method == KINASE_ACTIVITY_METHOD_KSEA_ZSCORE:
        if quantity in {
            QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
            QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
        }:
            return ActivityInputMatrix.sample_level_abundance(
                matrix,
                field_name="kinase_request.activity_phospho_matrix",
                _assume_owned=True,
            )
        if quantity is QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE:
            return ActivityInputMatrix.contrast_log_fold_change(
                matrix,
                field_name="kinase_request.activity_phospho_matrix",
                _assume_owned=True,
            )
        if quantity is QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE:
            return ActivityInputMatrix.standardised_effect(
                matrix,
                field_name="kinase_request.activity_phospho_matrix",
                _assume_owned=True,
            )
        observed = None if quantity is None else quantity.value
        raise WorkflowBoundaryError(
            seam="kinase.activity.method_quantitative_input_contract",
            next_action=(
                "provide phosphosite_log_abundance, phospho_total_log_ratio, "
                "contrast_log2_fold_change, or differential_effect_size values "
                "for KSEA-style activity; the workflow does not transform "
                "linear abundance data"
            ),
            details={"method": str(method), "quantitative_meaning": observed},
            message_prefix="kinase workflow boundary validation failed",
        )
    if method == KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT:
        if quantity is QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE:
            return ActivityInputMatrix.contrast_log_fold_change(
                matrix,
                field_name="kinase_request.activity_phospho_matrix",
                _assume_owned=True,
            )
        if quantity is QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE:
            return ActivityInputMatrix.standardised_effect(
                matrix,
                field_name="kinase_request.activity_phospho_matrix",
                _assume_owned=True,
            )
        observed = None if quantity is None else quantity.value
        raise WorkflowBoundaryError(
            seam="kinase.activity.method_quantitative_input_contract",
            next_action=(
                "provide contrast_log2_fold_change or differential_effect_size "
                "values for ssGSEA-style activity; the workflow does not "
                "transform abundance data into effect values"
            ),
            details={"method": str(method), "quantitative_meaning": observed},
            message_prefix="kinase workflow boundary validation failed",
        )
    return ActivityInputMatrix.sample_level_abundance(
        matrix,
        field_name="kinase_request.activity_phospho_matrix",
        _assume_owned=True,
    )


def _validate_prediction_membership_universe(
    *,
    pred_mat: pd.DataFrame,
    request: ResolvedKinaseWorkflowRequest,
    method_id: str,
) -> None:
    if request.site_universes is None:
        raise WorkflowBoundaryError(
            seam="kinase.activity.site_universes",
            next_action="ensure interpreter resolves typed kinase site universes",
            details={"method_id": method_id},
            message_prefix="kinase workflow boundary validation failed",
        )
    allowed = set(
        request.site_universes.predicted_membership_sites.astype(str).tolist()
    )
    unexpected = [
        str(site_id)
        for site_id in pred_mat.index.astype(str).tolist()
        if str(site_id) not in allowed
    ]
    if unexpected:
        raise WorkflowBoundaryError(
            seam="kinase.activity.predicted_membership_universe",
            next_action=(
                "ensure prediction_result.pred_mat rows come from the resolved "
                "predicted_membership_sites universe"
            ),
            details={
                "method_id": method_id,
                "unexpected_site_examples": unexpected[:5],
            },
            message_prefix="kinase workflow boundary validation failed",
        )


def _prediction_membership_matrix_for_site_universe(
    *,
    pred_mat: pd.DataFrame,
    site_index: pd.Index,
) -> pd.DataFrame:
    return pred_mat.reindex(index=site_index.copy())


def _kinase_substrate_membership_for_site_universe(
    *,
    membership: pd.DataFrame,
    site_index: pd.Index,
) -> pd.DataFrame:
    if "substrate_site" not in membership.columns:
        return membership.iloc[0:0].copy(deep=True)
    site_values = set(str(value) for value in site_index.astype(str).tolist())
    substrate_sites = membership.loc[:, "substrate_site"].astype(str)
    return membership.loc[substrate_sites.isin(site_values), :].copy(deep=True)


__all__ = ["KinaseActivityRunner"]
