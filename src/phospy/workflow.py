from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd

from .motifs import MotifScoringResult
from .prediction import KinasePredictionResult, KinasePredictor
from .profiles import KinaseProfileResult, build_kinase_substrate_profiles
from .scoring import KinaseScorer, KinaseScoringResult
from .types import PredictionSvmMode
from .validation.workflow import (
    ValidatedWorkflowRequest,
    validate_workflow_request,
)


@dataclass(slots=True)
class KinaseWorkflowResult:
    """Detached snapshot bundle for native workflow outputs.

    The nested result objects are produced outputs of a workflow run. They are
    not live views into workflow request inputs or mutable scorer/predictor
    internals.
    """

    profile_result: KinaseProfileResult
    motif_result: MotifScoringResult | None
    scoring_result: KinaseScoringResult
    prediction_result: KinasePredictionResult


@dataclass(frozen=True, slots=True)
class WorkflowExecutionPlan:
    request: ValidatedWorkflowRequest
    kernel: str
    predictor_svm_mode: PredictionSvmMode

    @property
    def validated_inputs(self) -> ValidatedWorkflowRequest:
        """Backward-compatible alias for older planner consumers."""
        return self.request


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

    def plan(
        self,
        request: ValidatedWorkflowRequest,
    ) -> WorkflowExecutionPlan:
        if not isinstance(request, ValidatedWorkflowRequest):
            msg = (
                "WorkflowExecutionPlanner.plan requires a ValidatedWorkflowRequest. "
                "Call KinaseWorkflow.validate_request(...) first."
            )
            raise TypeError(msg)
        return WorkflowExecutionPlan(
            request=request,
            kernel=self.kernel,
            predictor_svm_mode=request.predictor_svm_mode,
        )


class WorkflowExecutionRunner:
    def execute(self, plan: WorkflowExecutionPlan) -> KinaseWorkflowResult:
        request = plan.request.request
        phospho_matrix = plan.request.phospho_matrix

        profile_result = build_kinase_substrate_profiles(
            substrate_map=request.substrate_map,
            phospho_matrix=phospho_matrix,
            min_substrates=request.min_substrates,
        )

        scorer = KinaseScorer(profile_result.profile_matrix)
        motif_result: MotifScoringResult | None = None

        if plan.request.motif_scorer is not None:
            motif_result = plan.request.motif_scorer.score_sequences(
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

    def validate_request(
        self,
        *,
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
    ) -> ValidatedWorkflowRequest:
        return validate_workflow_request(
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
            flank_size=self.flank_size,
            default_svm_mode=self.svm_mode,
        )

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
        request = self.validate_request(
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

    def run_request(
        self,
        request: ValidatedWorkflowRequest,
    ) -> KinaseWorkflowResult:
        if not isinstance(request, ValidatedWorkflowRequest):
            msg = (
                "run_request requires a ValidatedWorkflowRequest. "
                "Call validate_request(...) first."
            )
            raise TypeError(msg)
        return self.execution_runner.execute(self.execution_planner.plan(request))
