from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd

from .motifs import MotifScoringResult
from .prediction import KinasePredictionResult, KinasePredictor, PredMatResult
from .profiles import KinaseProfileResult, build_kinase_substrate_profiles
from .scoring import KinaseScorer, KinaseScoringResult
from .signalomes import SignalomeResult, build_signalome_result
from .types import PredictionSvmMode
from .validation.signalomes import (
    ValidatedSignalomeRequest,
    validate_signalome_request,
)
from .validation.workflow import (
    ValidatedWorkflowRequest,
    validate_workflow_request,
)

__all__ = [
    "KinaseWorkflow",
    "KinaseWorkflowResult",
    "PredMatWorkflow",
    "PredMatWorkflowResult",
    "SignalomeWorkflow",
]


@dataclass(slots=True)
class KinaseWorkflowResult:
    """Workflow outputs for a single native scoring and prediction run."""

    profile_result: KinaseProfileResult
    motif_result: MotifScoringResult | None
    scoring_result: KinaseScoringResult
    prediction_result: KinasePredictionResult


@dataclass(slots=True)
class PredMatWorkflowResult:
    """Stable result bundle for one public predMat generation run.

    The canonical predMat contract is exposed through ``pred_mat_result``.
    """

    scoring_result: KinaseScoringResult
    prediction_result: KinasePredictionResult
    pred_mat_result: PredMatResult

    def close(self) -> None:
        """Release owned trace resources, if any are attached downstream."""

        self.prediction_result.close()

    def __enter__(self) -> PredMatWorkflowResult:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: object | None,
    ) -> None:
        self.close()


@dataclass(frozen=True, slots=True)
class _WorkflowPlan:
    request: ValidatedWorkflowRequest
    kernel: str
    predictor_svm_mode: PredictionSvmMode


class _WorkflowPlanner:
    def __init__(self, *, kernel: str) -> None:
        self.kernel = kernel

    def plan(
        self,
        request: ValidatedWorkflowRequest,
    ) -> _WorkflowPlan:
        return _WorkflowPlan(
            request=request,
            kernel=self.kernel,
            predictor_svm_mode=request.predictor_svm_mode,
        )


class _WorkflowRunner:
    def execute(self, plan: _WorkflowPlan) -> KinaseWorkflowResult:
        request = plan.request.request
        phospho_matrix = plan.request.phospho_matrix

        profile_result = build_kinase_substrate_profiles(
            substrate_map=request.substrate_map,
            phospho_matrix=phospho_matrix,
            min_substrates=request.min_substrates,
        )

        scorer = KinaseScorer(profile_result.profile_matrix)
        motif_result: MotifScoringResult | None = None
        scoring_matrix = phospho_matrix

        if plan.request.motif_scorer is not None:
            scoring_matrix = phospho_matrix.loc[list(plan.request.scoring_site_index)]
            motif_result = plan.request.motif_scorer.score_sequences(
                seqs=request.site_sequences,
                site_index=plan.request.scoring_site_index,
                min_motif_size=request.min_motif_size,
            )
            scoring_result = scorer.score(
                phospho_matrix=scoring_matrix,
                motif_scores=motif_result.motif_scores,
                motif_sizes=motif_result.motif_sizes,
                profile_sizes=profile_result.substrate_counts.astype(float),
                allow_profile_only_fallback=request.allow_profile_only_fallback,
            )
        else:
            scoring_result = scorer.score(phospho_matrix=scoring_matrix)

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


class _WorkflowFacadeBase:
    def __init__(
        self,
        flank_size: int = 7,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = "default",
    ) -> None:
        self.flank_size = flank_size
        self.kernel = kernel
        self.svm_mode = svm_mode
        self.planner = _WorkflowPlanner(kernel=kernel)
        self.runner = _WorkflowRunner()

    def _validate_request(
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

    def _execute_request(
        self,
        request: ValidatedWorkflowRequest,
    ) -> KinaseWorkflowResult:
        return self.runner.execute(self.planner.plan(request))


class KinaseWorkflow(_WorkflowFacadeBase):
    """Run the native kinase scoring and prediction workflow end to end."""

    def __init__(
        self,
        flank_size: int = 7,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = "default",
    ) -> None:
        super().__init__(
            flank_size=flank_size,
            kernel=kernel,
            svm_mode=svm_mode,
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
        request = self._validate_request(
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
        return self._run_request(request)

    def _run_request(
        self,
        request: ValidatedWorkflowRequest,
    ) -> KinaseWorkflowResult:
        return self._execute_request(request)


class PredMatWorkflow(_WorkflowFacadeBase):
    """Generate a predMat from phosphosite and sequence inputs."""

    def __init__(
        self,
        flank_size: int = 7,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = "default",
    ) -> None:
        super().__init__(
            flank_size=flank_size,
            kernel=kernel,
            svm_mode=svm_mode,
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
    ) -> PredMatWorkflowResult:
        request = self._validate_request(
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
        return self._run_request(request)

    def _run_request(
        self,
        request: ValidatedWorkflowRequest,
    ) -> PredMatWorkflowResult:
        result = self._execute_request(request)
        return PredMatWorkflowResult(
            scoring_result=result.scoring_result,
            prediction_result=result.prediction_result,
            pred_mat_result=result.prediction_result.pred_mat_result,
        )


class SignalomeWorkflow:
    """Construct signalomes from validated scoring and prediction outputs."""

    def _validate_request(
        self,
        *,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult | PredMatResult,
        expression_matrix: pd.DataFrame,
        kinases_of_interest: Sequence[str],
        site_to_protein: Mapping[str, str] | None = None,
        kinase_network_threshold: float = 0.9,
        signalome_cutoff: float = 0.5,
        module_count: int | None = None,
        min_kinase_module_share_percent: float = 1.0,
    ) -> ValidatedSignalomeRequest:
        return validate_signalome_request(
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            expression_matrix=expression_matrix,
            kinases_of_interest=kinases_of_interest,
            site_to_protein=site_to_protein,
            kinase_network_threshold=kinase_network_threshold,
            signalome_cutoff=signalome_cutoff,
            module_count=module_count,
            min_kinase_module_share_percent=min_kinase_module_share_percent,
        )

    def run(
        self,
        *,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult | PredMatResult,
        expression_matrix: pd.DataFrame,
        kinases_of_interest: Sequence[str],
        site_to_protein: Mapping[str, str] | None = None,
        kinase_network_threshold: float = 0.9,
        signalome_cutoff: float = 0.5,
        module_count: int | None = None,
        min_kinase_module_share_percent: float = 1.0,
    ) -> SignalomeResult:
        request = self._validate_request(
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            expression_matrix=expression_matrix,
            kinases_of_interest=kinases_of_interest,
            site_to_protein=site_to_protein,
            kinase_network_threshold=kinase_network_threshold,
            signalome_cutoff=signalome_cutoff,
            module_count=module_count,
            min_kinase_module_share_percent=min_kinase_module_share_percent,
        )
        return self._run_request(request)

    def _run_request(
        self,
        request: ValidatedSignalomeRequest,
    ) -> SignalomeResult:
        return build_signalome_result(
            scoring_matrix=request.scoring_matrix,
            pred_mat=request.pred_mat,
            expression_matrix=request.expression_matrix,
            kinases_of_interest=request.request.kinases_of_interest,
            site_to_protein=request.site_to_protein,
            kinase_network_threshold=request.request.kinase_network_threshold,
            signalome_cutoff=request.request.signalome_cutoff,
            module_count=request.request.module_count,
            min_kinase_module_share_percent=(
                request.request.min_kinase_module_share_percent
            ),
        )
