from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd

from .motifs import MotifScoringResult
from .prediction import KinasePredictionResult, KinasePredictor, PredMatResult
from .profiles import KinaseProfileResult, build_kinase_substrate_profiles
from .scoring import KinaseScorer, KinaseScoringResult
from .signalome_construction import execute_validated_signalome_request
from .signalomes import SignalomeResult
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

    The recommended predMat contract is exposed through ``pred_mat_result``.
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


class _WorkflowBase:
    def __init__(
        self,
        flank_size: int = 7,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = "default",
    ) -> None:
        self.flank_size = flank_size
        self.kernel = kernel
        self.svm_mode = svm_mode

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

    def _execute_validated_request(
        self,
        request: ValidatedWorkflowRequest,
    ) -> KinaseWorkflowResult:
        raw_request = request.request
        phospho_matrix = request.phospho_matrix

        profile_result = build_kinase_substrate_profiles(
            substrate_map=raw_request.substrate_map,
            phospho_matrix=phospho_matrix,
            min_substrates=raw_request.min_substrates,
        )

        scorer = KinaseScorer(profile_result.profile_matrix)
        motif_result: MotifScoringResult | None = None
        scoring_matrix = phospho_matrix

        if request.motif_scorer is not None:
            scoring_matrix = phospho_matrix.loc[list(request.scoring_site_index)]
            motif_result = request.motif_scorer.score_sequences(
                seqs=raw_request.site_sequences,
                site_index=request.scoring_site_index,
                min_motif_size=raw_request.min_motif_size,
            )
            scoring_result = scorer.score(
                phospho_matrix=scoring_matrix,
                motif_scores=motif_result.motif_scores,
                motif_sizes=motif_result.motif_sizes,
                profile_sizes=profile_result.substrate_counts.astype(float),
                allow_profile_only_fallback=raw_request.allow_profile_only_fallback,
            )
        else:
            scoring_result = scorer.score(phospho_matrix=scoring_matrix)

        predictor = KinasePredictor(
            kernel=self.kernel,
            svm_mode=request.predictor_svm_mode,
        )
        prediction_result = predictor.predict_from_scoring_result(
            scoring_result=scoring_result,
            ensemble_size=raw_request.ensemble_size,
            top=raw_request.top,
            score_threshold=raw_request.score_threshold,
            inclusion=raw_request.inclusion,
            n_iterations=raw_request.n_iterations,
            random_state=raw_request.random_state,
            allow_profile_only_fallback=raw_request.allow_profile_only_fallback,
            svm_mode=raw_request.svm_mode,
        )

        return KinaseWorkflowResult(
            profile_result=profile_result,
            motif_result=motif_result,
            scoring_result=scoring_result,
            prediction_result=prediction_result,
        )


class KinaseWorkflow(_WorkflowBase):
    """Run the native kinase scoring and prediction workflow end to end."""

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
        return self.run_validated(request)

    def run_validated(
        self,
        request: ValidatedWorkflowRequest,
    ) -> KinaseWorkflowResult:
        return self._execute_validated_request(request)


class PredMatWorkflow(_WorkflowBase):
    """Generate a predMat from phosphosite and sequence inputs."""

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
        return self.run_validated(request)

    def run_validated(
        self,
        request: ValidatedWorkflowRequest,
    ) -> PredMatWorkflowResult:
        result = self._execute_validated_request(request)
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
        return self.run_validated(request)

    def run_validated(
        self,
        request: ValidatedSignalomeRequest,
    ) -> SignalomeResult:
        return execute_validated_signalome_request(request)
