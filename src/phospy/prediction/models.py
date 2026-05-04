"""Prediction and scoring stage result models."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from phospy._frame_ownership import (
    _borrow_dataframe,
    _borrow_optional_dataframe,
    export_dataframe,
    export_optional_dataframe,
    own_optional_dataframe,
)
from phospy.errors.validation import PhosPyValidationError
from phospy.prediction.motif_scoring import MotifLibraryValidationResult
from phospy.prediction.sequence_validation import SequenceValidationResult
from phospy.tables.kinase import KinasePredictionMatrix, KinaseScoreMatrix


@dataclass(frozen=True, slots=True, init=False)
class KinaseScoringResult:
    """Scoring-stage outputs.

    `profile_scores` and `rank_weighted_fusion_scores` define the supported
    downstream lane. `motif_scores` and `score_fusion_weights` are optional
    diagnostic tables controlled by
    `scoring_config.include_diagnostic_scoring_tables`.
    """

    motif_sequence_validation: SequenceValidationResult | None = None
    motif_library_validation: MotifLibraryValidationResult | None = None
    _profile_scores: pd.DataFrame = field(init=False, repr=False)
    _motif_scores: pd.DataFrame | None = field(init=False, repr=False)
    _rank_weighted_fusion_scores: pd.DataFrame | None = field(init=False, repr=False)
    _score_fusion_weights: pd.DataFrame | None = field(init=False, repr=False)

    def __init__(
        self,
        profile_scores: pd.DataFrame,
        motif_scores: pd.DataFrame | None = None,
        rank_weighted_fusion_scores: pd.DataFrame | None = None,
        score_fusion_weights: pd.DataFrame | None = None,
        motif_sequence_validation: SequenceValidationResult | None = None,
        motif_library_validation: MotifLibraryValidationResult | None = None,
        _assume_owned: bool = False,
    ) -> None:
        object.__setattr__(self, "motif_sequence_validation", motif_sequence_validation)
        object.__setattr__(self, "motif_library_validation", motif_library_validation)
        profile_scores = KinaseScoreMatrix(
            frame=profile_scores,
            field_name="scoring_result.profile_scores",
            _assume_owned=_assume_owned,
        ).frame
        motif_scores = (
            None
            if motif_scores is None
            else KinaseScoreMatrix(
                frame=motif_scores,
                field_name="scoring_result.motif_scores",
                _assume_owned=_assume_owned,
            ).frame
        )
        rank_weighted_fusion_scores = (
            None
            if rank_weighted_fusion_scores is None
            else KinaseScoreMatrix(
                frame=rank_weighted_fusion_scores,
                field_name="scoring_result.rank_weighted_fusion_scores",
                _assume_owned=_assume_owned,
            ).frame
        )
        score_fusion_weights = (
            None
            if score_fusion_weights is None
            else KinaseScoreMatrix(
                frame=score_fusion_weights,
                field_name="scoring_result.score_fusion_weights",
                require_site_index=False,
                _assume_owned=_assume_owned,
            ).frame
        )
        object.__setattr__(self, "_profile_scores", profile_scores)
        object.__setattr__(self, "_motif_scores", motif_scores)
        object.__setattr__(
            self, "_rank_weighted_fusion_scores", rank_weighted_fusion_scores
        )
        object.__setattr__(self, "_score_fusion_weights", score_fusion_weights)
        if motif_sequence_validation is not None and not isinstance(
            motif_sequence_validation,
            SequenceValidationResult,
        ):
            raise PhosPyValidationError(
                "scoring_result.motif_sequence_validation must be "
                "SequenceValidationResult or None"
            )
        if motif_library_validation is not None and not isinstance(
            motif_library_validation,
            MotifLibraryValidationResult,
        ):
            raise PhosPyValidationError(
                "scoring_result.motif_library_validation must be "
                "MotifLibraryValidationResult or None"
            )

    @property
    def profile_scores(self) -> pd.DataFrame:
        return export_dataframe(self._profile_scores)

    @property
    def motif_scores(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._motif_scores)

    @property
    def rank_weighted_fusion_scores(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._rank_weighted_fusion_scores)

    @property
    def score_fusion_weights(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._score_fusion_weights)

    def _borrow_profile_scores_frame(self) -> pd.DataFrame:
        """Package-private borrowed profile scores for internal workflows."""

        return _borrow_dataframe(self._profile_scores)

    def _borrow_motif_scores_frame(self) -> pd.DataFrame | None:
        """Package-private borrowed motif scores for internal workflows."""

        return _borrow_optional_dataframe(self._motif_scores)

    def _borrow_rank_weighted_fusion_scores_frame(self) -> pd.DataFrame | None:
        """Package-private borrowed fusion scores for internal workflows."""

        return _borrow_optional_dataframe(self._rank_weighted_fusion_scores)

    def _borrow_score_fusion_weights_frame(self) -> pd.DataFrame | None:
        """Package-private borrowed fusion weights for internal workflows."""

        return _borrow_optional_dataframe(self._score_fusion_weights)

    @classmethod
    def _from_owned(
        cls,
        *,
        profile_scores: pd.DataFrame,
        motif_scores: pd.DataFrame | None = None,
        rank_weighted_fusion_scores: pd.DataFrame | None = None,
        score_fusion_weights: pd.DataFrame | None = None,
        motif_sequence_validation: SequenceValidationResult | None = None,
        motif_library_validation: MotifLibraryValidationResult | None = None,
    ) -> KinaseScoringResult:
        return cls(
            profile_scores=profile_scores,
            motif_scores=motif_scores,
            rank_weighted_fusion_scores=rank_weighted_fusion_scores,
            score_fusion_weights=score_fusion_weights,
            motif_sequence_validation=motif_sequence_validation,
            motif_library_validation=motif_library_validation,
            _assume_owned=True,
        )

    def to_dataframe(self) -> pd.DataFrame:
        """Return a `profile_scores` snapshot isolated from this result."""

        return export_dataframe(self._profile_scores)

    def motif_scores_dataframe(self) -> pd.DataFrame | None:
        """Return an optional motif-score snapshot isolated from this result."""

        return export_optional_dataframe(self._motif_scores)

    def rank_weighted_fusion_scores_dataframe(self) -> pd.DataFrame | None:
        """Return an optional fusion-score snapshot isolated from this result."""

        return export_optional_dataframe(self._rank_weighted_fusion_scores)

    def score_fusion_weights_dataframe(self) -> pd.DataFrame | None:
        """Return an optional fusion-weight snapshot isolated from this result."""

        return export_optional_dataframe(self._score_fusion_weights)


@dataclass(frozen=True, slots=True, init=False)
class KinasePredictionResult:
    """Prediction-stage outputs."""

    _pred_mat: pd.DataFrame = field(init=False, repr=False)
    _substrate_list: pd.DataFrame | None = field(init=False, repr=False)

    def __init__(
        self,
        pred_mat: pd.DataFrame,
        substrate_list: pd.DataFrame | None = None,
        _assume_owned: bool = False,
    ) -> None:
        pred_mat = KinasePredictionMatrix(
            frame=pred_mat,
            field_name="prediction_result.pred_mat",
            _assume_owned=_assume_owned,
        ).frame
        substrate_list = own_optional_dataframe(
            substrate_list,
            field_name="prediction_result.substrate_list",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        object.__setattr__(self, "_pred_mat", pred_mat)
        object.__setattr__(self, "_substrate_list", substrate_list)

    @property
    def pred_mat(self) -> pd.DataFrame:
        return export_dataframe(self._pred_mat)

    @property
    def substrate_list(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._substrate_list)

    def _borrow_pred_mat_frame(self) -> pd.DataFrame:
        """Package-private borrowed prediction matrix for internal workflows."""

        return _borrow_dataframe(self._pred_mat)

    def _borrow_substrate_list_frame(self) -> pd.DataFrame | None:
        """Package-private borrowed substrate list for internal workflows."""

        return _borrow_optional_dataframe(self._substrate_list)

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

    def to_dataframe(self) -> pd.DataFrame:
        """Return a `pred_mat` snapshot isolated from this result."""

        return export_dataframe(self._pred_mat)

    def substrate_list_dataframe(self) -> pd.DataFrame | None:
        """Return an optional substrate-list snapshot isolated from this result."""

        return export_optional_dataframe(self._substrate_list)
