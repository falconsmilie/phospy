"""Internal executor for the kinase workflow."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.activities.models import KinaseActivityResult
from phospy.activities.scoring import compute_activity_from_inputs
from phospy.api.results import KinaseWorkflowResult
from phospy.errors.workflows import WorkflowBoundaryError, WorkflowStageError
from phospy.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.validation.workflows.activity import KinaseActivityInputValidator
from phospy.workflows.kinase.contracts import ResolvedKinaseWorkflowRequest
from phospy.workflows.kinase.science import (
    build_kinase_profiles,
    build_prediction_outputs,
    rank_kinases_for_prediction,
    score_profile_correlations,
)


@dataclass(frozen=True, slots=True)
class _ScoringExecution:
    scoring_result: KinaseScoringResult
    score_matrix: pd.DataFrame
    quantified_substrates: dict[str, list[str]]


class KinaseWorkflowExecutor:
    """Run stage logic and assemble `KinaseWorkflowResult`."""

    def __init__(
        self,
        *,
        activity_input_validator: KinaseActivityInputValidator | None = None,
    ) -> None:
        self._activity_input_validator = (
            activity_input_validator or KinaseActivityInputValidator()
        )

    def run(self, request: ResolvedKinaseWorkflowRequest) -> KinaseWorkflowResult:
        scoring_execution = self._run_scoring_stage(request)
        prediction_result = self._run_prediction_stage(
            request=request,
            scoring_execution=scoring_execution,
        )
        activity_result = self._run_activity_stage(
            request=request,
            prediction_result=prediction_result,
        )
        return KinaseWorkflowResult(
            dataset=request.dataset,
            references=request.references,
            scoring_result=scoring_execution.scoring_result,
            prediction_result=prediction_result,
            activity_result=activity_result,
        )

    def _run_scoring_stage(
        self, request: ResolvedKinaseWorkflowRequest
    ) -> _ScoringExecution:
        # Current supported route is profile-based scoring over sequence-backed
        # phosphosite rows. Motif/combined scoring lanes are intentionally not
        # materialized here.
        scoring_phospho = request.dataset.phospho.loc[request.scoring_site_index, :]
        profile_build = build_kinase_profiles(
            phospho=scoring_phospho,
            kinase_substrate_map=request.kinase_substrate_map,
            min_substrates=request.scoring_config.min_substrates,
        )
        if profile_build.profile_matrix.empty:
            raise WorkflowStageError(
                "kinase workflow internal invariant failed at seam="
                "kinase.executor.scoring_profiles; interpreter should reject "
                "requests with zero eligible kinases before scoring"
            )
        profile_scores = score_profile_correlations(
            phospho=scoring_phospho,
            profile_matrix=profile_build.profile_matrix,
        )
        scoring_result = KinaseScoringResult(profile_scores=profile_scores)
        return _ScoringExecution(
            scoring_result=scoring_result,
            score_matrix=profile_scores,
            quantified_substrates=profile_build.quantified_substrates,
        )

    def _run_prediction_stage(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        scoring_execution: _ScoringExecution,
    ) -> KinasePredictionResult:
        kinase_ranking = rank_kinases_for_prediction(
            score_matrix=scoring_execution.score_matrix,
            quantified_substrates=scoring_execution.quantified_substrates,
        )
        selected_kinases = kinase_ranking.head(
            request.prediction_config.ensemble_size
        ).index
        if selected_kinases.empty:
            self._raise_boundary_error(
                seam="kinase.executor.prediction_ensemble",
                next_action=(
                    "provide dataset.phospho with at least two non-constant "
                    "sample columns or lower scoring_config.min_substrates "
                    "(scientific floor: min_substrates >= 2)"
                ),
                eligible_kinases=len(scoring_execution.quantified_substrates),
                ranked_kinases=int(kinase_ranking.size),
                prediction_config_ensemble_size=request.prediction_config.ensemble_size,
                prediction_config_top_k=request.prediction_config.top_k,
                dataset_samples=request.dataset.phospho.shape[1],
            )
        pred_mat, substrate_list = build_prediction_outputs(
            score_matrix=scoring_execution.score_matrix,
            selected_kinases=selected_kinases,
            quantified_substrates=scoring_execution.quantified_substrates,
            top_k=request.prediction_config.top_k,
        )
        return KinasePredictionResult(
            pred_mat=pred_mat,
            substrate_list=substrate_list,
        )

    def _run_activity_stage(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        prediction_result: KinasePredictionResult,
    ) -> KinaseActivityResult | None:
        activity_config = request.activity_config
        if activity_config is None or not activity_config.enabled:
            return None
        validated_inputs = self._activity_input_validator.run(
            pred_mat=prediction_result.pred_mat,
            phospho_matrix=request.activity_phospho_matrix,
            threshold=activity_config.threshold,
            min_substrates=activity_config.min_substrates,
            top_n_substrates=activity_config.top_n_substrates,
        )
        return compute_activity_from_inputs(validated_inputs)

    @staticmethod
    def _raise_boundary_error(
        *,
        seam: str,
        next_action: str,
        **details: int | float,
    ) -> None:
        details_text = ", ".join(f"{key}={value}" for key, value in details.items())
        raise WorkflowBoundaryError(
            "kinase workflow boundary validation failed at "
            f"seam={seam}; {details_text}; next_action={next_action}"
        )
