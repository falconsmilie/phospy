from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd

from ..motifs import MotifScoringResult, ReferenceBundle
from ..profiles import KinaseProfileResult, build_kinase_substrate_profiles
from ..scoring import KinaseScorer, KinaseScoringResult
from ..types import PredictionSvmMode
from ..validation.requests import ValidatedWorkflowRequest, validate_workflow_request
from .models import KinasePredictionResult
from .service import KinasePredictor


@dataclass(frozen=True, slots=True)
class KinaseWorkflowExecutionResult:
    """Trusted execution outputs for one validated kinase workflow request."""

    profile_result: KinaseProfileResult
    motif_result: MotifScoringResult | None
    scoring_result: KinaseScoringResult
    prediction_result: KinasePredictionResult


class KinaseWorkflowExecutor:
    """Prediction-domain executor for the validated kinase workflow path."""

    def __init__(
        self,
        *,
        flank_size: int = 7,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = "default",
    ) -> None:
        self.flank_size = flank_size
        self.kernel = kernel
        self.svm_mode = svm_mode

    def validate_request(
        self,
        *,
        phospho_matrix: pd.DataFrame,
        substrate_map: Mapping[str, Sequence[str]] | None = None,
        site_sequences: Mapping[str, str] | pd.Series | None = None,
        motif_sequences: Mapping[str, Sequence[str]] | None = None,
        reference_bundle: ReferenceBundle | None = None,
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
            reference_bundle=reference_bundle,
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

    def execute_validated_request(
        self,
        request: ValidatedWorkflowRequest,
    ) -> KinaseWorkflowExecutionResult:
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

        return KinaseWorkflowExecutionResult(
            profile_result=profile_result,
            motif_result=motif_result,
            scoring_result=scoring_result,
            prediction_result=prediction_result,
        )
