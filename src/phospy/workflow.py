from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd

from .motifs import MotifScoringResult
from .prediction import KinasePredictionResult, KinasePredictor
from .profiles import KinaseProfileResult, build_kinase_substrate_profiles
from .scoring import KinaseScorer, KinaseScoringResult
from .types import PredictionSvmMode
from .validation.compatibility import (
    ValidatedKinaseWorkflowInputs,
    build_workflow_request_inputs,
)
from .validation.requests import KinaseWorkflowRequest


@dataclass(slots=True)
class KinaseWorkflowResult:
    profile_result: KinaseProfileResult
    motif_result: MotifScoringResult | None
    scoring_result: KinaseScoringResult
    prediction_result: KinasePredictionResult


@dataclass(frozen=True, slots=True)
class WorkflowExecutionPlan:
    validated_inputs: ValidatedKinaseWorkflowInputs
    predictor_svm_mode: PredictionSvmMode
    kernel: str


class WorkflowExecutionPlanner:
    def __init__(
        self,
        *,
        flank_size: int,
        kernel: str,
        default_svm_mode: PredictionSvmMode,
    ) -> None:
        self.flank_size = flank_size
        self.kernel = kernel
        self.default_svm_mode = default_svm_mode

    def plan(self, request: KinaseWorkflowRequest) -> WorkflowExecutionPlan:
        validated_inputs = build_workflow_request_inputs(
            request,
            flank_size=self.flank_size,
        )
        return WorkflowExecutionPlan(
            validated_inputs=validated_inputs,
            predictor_svm_mode=(
                self.default_svm_mode if request.svm_mode is None else request.svm_mode
            ),
            kernel=self.kernel,
        )


class WorkflowExecutionRunner:
    def execute(self, plan: WorkflowExecutionPlan) -> KinaseWorkflowResult:
        request = plan.validated_inputs.request
        phospho_matrix = plan.validated_inputs.phospho_matrix

        profile_result = build_kinase_substrate_profiles(
            substrate_map=request.substrate_map,
            phospho_matrix=phospho_matrix,
            min_substrates=request.min_substrates,
        )

        scorer = KinaseScorer(profile_result.profile_matrix)
        motif_result: MotifScoringResult | None = None

        if plan.validated_inputs.motif_scorer is not None:
            motif_result = plan.validated_inputs.motif_scorer.score_sequences(
                seqs=request.site_sequences,
                site_index=phospho_matrix.index,
                min_motif_size=request.min_motif_size,
            )
            scoring_result = scorer.score(
                phospho_matrix=phospho_matrix,
                motif_scores=motif_result.motif_scores,
                motif_sizes=motif_result.motif_sizes,
                profile_sizes=profile_result.substrate_counts.astype(float),
                allow_profile_only_fallback=request.allow_profile_only_fallback,
            )
        else:
            scoring_result = scorer.score(phospho_matrix=phospho_matrix)

        predictor = KinasePredictor(
            kernel=plan.kernel,
            svm_mode=plan.predictor_svm_mode,
        )
        prediction_result = predictor.predict_from_scoring_result(
            scoring_result=scoring_result,
            ensemble_size=request.ensemble_size,
            top=request.top,
            score_threshold=request.score_threshold,
            inclusion=request.inclusion,
            n_iterations=request.n_iterations,
            random_state=request.random_state,
            allow_profile_only_fallback=request.allow_profile_only_fallback,
            svm_mode=request.svm_mode,
        )

        return KinaseWorkflowResult(
            profile_result=profile_result,
            motif_result=motif_result,
            scoring_result=scoring_result,
            prediction_result=prediction_result,
        )


class KinaseWorkflow:
    """Run the native kinase scoring and prediction workflow end to end."""

    def __init__(
        self,
        flank_size: int = 7,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = "default",
    ) -> None:
        self.flank_size = flank_size
        self.kernel = kernel
        self.svm_mode = svm_mode
        self.execution_planner = WorkflowExecutionPlanner(
            flank_size=flank_size,
            kernel=kernel,
            default_svm_mode=svm_mode,
        )
        self.execution_runner = WorkflowExecutionRunner()

    def run(
        self,
        phospho_matrix: pd.DataFrame,
        substrate_map: Mapping[str, Sequence[str]],
        site_sequences: Mapping[str, str] | pd.Series | None = None,
        motif_sequences: Mapping[str, Sequence[str]] | None = None,
        min_substrates: int = 1,
        min_motif_size: int = 1,
        allow_profile_only_fallback: bool = False,
        ensemble_size: int = 10,
        top: int = 50,
        score_threshold: float = 0.8,
        inclusion: int = 20,
        n_iterations: int = 5,
        random_state: int | None = None,
        svm_mode: PredictionSvmMode | None = None,
    ) -> KinaseWorkflowResult:
        request = KinaseWorkflowRequest.validate_request(
            phospho_matrix=phospho_matrix,
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            motif_sequences=motif_sequences,
            min_substrates=min_substrates,
            min_motif_size=min_motif_size,
            allow_profile_only_fallback=allow_profile_only_fallback,
            ensemble_size=ensemble_size,
            top=top,
            score_threshold=score_threshold,
            inclusion=inclusion,
            n_iterations=n_iterations,
            random_state=random_state,
            svm_mode=svm_mode,
        )
        return self.run_request(request)

    def run_request(self, request: KinaseWorkflowRequest) -> KinaseWorkflowResult:
        return self.execution_runner.execute(self.execution_planner.plan(request))
