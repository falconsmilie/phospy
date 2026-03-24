from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd

from .motifs import KinaseMotifScorer, MotifScoringResult
from .prediction import (
    KinasePredictionResult,
    KinasePredictor,
    PredictionSvmMode,
)
from .profiles import (
    AggregationMethod,
    KinaseProfileBuilder,
    KinaseProfileResult,
)
from .scoring import KinaseScorer, KinaseScoringResult


@dataclass(slots=True)
class KinaseWorkflowResult:
    profile_result: KinaseProfileResult
    motif_result: MotifScoringResult | None
    scoring_result: KinaseScoringResult
    prediction_result: KinasePredictionResult


class KinaseWorkflow:
    """Run the native kinase scoring and prediction workflow end to end."""

    def __init__(
        self,
        aggregation: AggregationMethod = "median",
        flank_size: int = 7,
        kernel: str = "rbf",
        svm_mode: PredictionSvmMode = "default",
    ) -> None:
        self.aggregation = aggregation
        self.flank_size = flank_size
        self.kernel = kernel
        self.svm_mode = svm_mode

    def run(
        self,
        phospho_matrix: pd.DataFrame,
        substrate_map: Mapping[str, Sequence[str]],
        site_sequences: Mapping[str, str] | Sequence[str] | pd.Series | None = None,
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
        if not substrate_map:
            msg = "substrate_map must not be empty"
            raise ValueError(msg)

        if motif_sequences is not None and not motif_sequences:
            msg = (
                "motif_sequences must not be empty; pass None and set "
                "allow_profile_only_fallback=True for profile-only prediction"
            )
            raise ValueError(msg)

        if motif_sequences is None and not allow_profile_only_fallback:
            msg = (
                "motif_sequences are required for end-to-end prediction unless "
                "allow_profile_only_fallback=True"
            )
            raise ValueError(msg)

        if motif_sequences is not None and site_sequences is None:
            msg = "site_sequences are required when motif_sequences are provided"
            raise ValueError(msg)

        profile_builder = KinaseProfileBuilder(aggregation=self.aggregation)
        profile_result = profile_builder.build(
            substrate_map=substrate_map,
            phospho_matrix=phospho_matrix,
            min_substrates=min_substrates,
        )

        scorer = KinaseScorer(profile_result.profile_matrix)
        motif_result: MotifScoringResult | None = None

        if motif_sequences is not None:
            motif_scorer = KinaseMotifScorer.from_substrate_sequences(
                motif_sequences=motif_sequences,
                flank_size=self.flank_size,
            )
            motif_result = motif_scorer.score_sequences(
                seqs=site_sequences,
                site_index=phospho_matrix.index,
                min_motif_size=min_motif_size,
            )
            scoring_result = scorer.score(
                phospho_matrix=phospho_matrix,
                motif_scores=motif_result.motif_scores,
                motif_sizes=motif_result.motif_sizes,
                profile_sizes=profile_result.substrate_counts.astype(float),
                allow_profile_only_fallback=allow_profile_only_fallback,
            )
        else:
            scoring_result = scorer.score(phospho_matrix=phospho_matrix)

        predictor = KinasePredictor(
            kernel=self.kernel,
            svm_mode=self.svm_mode if svm_mode is None else svm_mode,
        )
        prediction_result = predictor.predict_from_scoring_result(
            scoring_result=scoring_result,
            ensemble_size=ensemble_size,
            top=top,
            score_threshold=score_threshold,
            inclusion=inclusion,
            n_iterations=n_iterations,
            random_state=random_state,
            allow_profile_only_fallback=allow_profile_only_fallback,
            svm_mode=svm_mode,
        )

        return KinaseWorkflowResult(
            profile_result=profile_result,
            motif_result=motif_result,
            scoring_result=scoring_result,
            prediction_result=prediction_result,
        )
