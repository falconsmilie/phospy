"""Prediction and scoring stage result models."""

from __future__ import annotations

from dataclasses import InitVar, dataclass

import pandas as pd

from phospy._frame_ownership import own_optional_dataframe
from phospy.errors.validation import PhosPyValidationError
from phospy.prediction.sequence_validation import SequenceValidationResult
from phospy.tables.kinase import KinasePredictionMatrix, KinaseScoreMatrix


@dataclass(frozen=True, slots=True)
class KinaseScoringResult:
    """Scoring-stage outputs.

    `profile_scores` and `rank_weighted_fusion_scores` define the supported
    downstream lane. `motif_scores` and `score_fusion_weights` are optional
    diagnostic tables controlled by
    `scoring_config.include_diagnostic_scoring_tables`.
    """

    profile_scores: pd.DataFrame
    motif_scores: pd.DataFrame | None = None
    rank_weighted_fusion_scores: pd.DataFrame | None = None
    score_fusion_weights: pd.DataFrame | None = None
    motif_sequence_validation: SequenceValidationResult | None = None
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        profile_scores = KinaseScoreMatrix(
            frame=self.profile_scores,
            field_name="scoring_result.profile_scores",
            _assume_owned=_assume_owned,
        ).frame
        motif_scores = (
            None
            if self.motif_scores is None
            else KinaseScoreMatrix(
                frame=self.motif_scores,
                field_name="scoring_result.motif_scores",
                _assume_owned=_assume_owned,
            ).frame
        )
        rank_weighted_fusion_scores = (
            None
            if self.rank_weighted_fusion_scores is None
            else KinaseScoreMatrix(
                frame=self.rank_weighted_fusion_scores,
                field_name="scoring_result.rank_weighted_fusion_scores",
                _assume_owned=_assume_owned,
            ).frame
        )
        score_fusion_weights = (
            None
            if self.score_fusion_weights is None
            else KinaseScoreMatrix(
                frame=self.score_fusion_weights,
                field_name="scoring_result.score_fusion_weights",
                _assume_owned=_assume_owned,
            ).frame
        )
        object.__setattr__(self, "profile_scores", profile_scores)
        object.__setattr__(self, "motif_scores", motif_scores)
        object.__setattr__(
            self, "rank_weighted_fusion_scores", rank_weighted_fusion_scores
        )
        object.__setattr__(self, "score_fusion_weights", score_fusion_weights)
        if self.motif_sequence_validation is not None and not isinstance(
            self.motif_sequence_validation,
            SequenceValidationResult,
        ):
            raise PhosPyValidationError(
                "scoring_result.motif_sequence_validation must be "
                "SequenceValidationResult or None"
            )

    @classmethod
    def _from_owned(
        cls,
        *,
        profile_scores: pd.DataFrame,
        motif_scores: pd.DataFrame | None = None,
        rank_weighted_fusion_scores: pd.DataFrame | None = None,
        score_fusion_weights: pd.DataFrame | None = None,
        motif_sequence_validation: SequenceValidationResult | None = None,
    ) -> KinaseScoringResult:
        return cls(
            profile_scores=profile_scores,
            motif_scores=motif_scores,
            rank_weighted_fusion_scores=rank_weighted_fusion_scores,
            score_fusion_weights=score_fusion_weights,
            motif_sequence_validation=motif_sequence_validation,
            _assume_owned=True,
        )


@dataclass(frozen=True, slots=True)
class KinasePredictionResult:
    """Prediction-stage outputs."""

    pred_mat: pd.DataFrame
    substrate_list: pd.DataFrame | None = None
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        pred_mat = KinasePredictionMatrix(
            frame=self.pred_mat,
            field_name="prediction_result.pred_mat",
            _assume_owned=_assume_owned,
        ).frame
        substrate_list = own_optional_dataframe(
            self.substrate_list,
            field_name="prediction_result.substrate_list",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        object.__setattr__(self, "pred_mat", pred_mat)
        object.__setattr__(self, "substrate_list", substrate_list)

    @classmethod
    def _from_owned(
        cls,
        *,
        pred_mat: pd.DataFrame,
        substrate_list: pd.DataFrame | None = None,
    ) -> KinasePredictionResult:
        return cls(
            pred_mat=pred_mat,
            substrate_list=substrate_list,
            _assume_owned=True,
        )
