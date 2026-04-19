"""Prediction and scoring stage result models."""

from __future__ import annotations

from dataclasses import InitVar, dataclass

import pandas as pd

from phospy._frame_ownership import own_dataframe, own_optional_dataframe
from phospy.errors.validation import PhosPyValidationError


@dataclass(frozen=True, slots=True)
class KinaseScoringResult:
    """Scoring-stage outputs.

    `profile_scores` and `combined_scores` define the supported downstream lane.
    `motif_scores` and `weights` are optional diagnostic tables controlled by
    `scoring_config.include_diagnostic_scoring_tables`.
    """

    profile_scores: pd.DataFrame
    motif_scores: pd.DataFrame | None = None
    combined_scores: pd.DataFrame | None = None
    weights: pd.DataFrame | None = None
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        profile_scores = own_dataframe(
            self.profile_scores,
            field_name="scoring_result.profile_scores",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        motif_scores = own_optional_dataframe(
            self.motif_scores,
            field_name="scoring_result.motif_scores",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        combined_scores = own_optional_dataframe(
            self.combined_scores,
            field_name="scoring_result.combined_scores",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        weights = own_optional_dataframe(
            self.weights,
            field_name="scoring_result.weights",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        object.__setattr__(self, "profile_scores", profile_scores)
        object.__setattr__(self, "motif_scores", motif_scores)
        object.__setattr__(self, "combined_scores", combined_scores)
        object.__setattr__(self, "weights", weights)

    @classmethod
    def _from_owned(
        cls,
        *,
        profile_scores: pd.DataFrame,
        motif_scores: pd.DataFrame | None = None,
        combined_scores: pd.DataFrame | None = None,
        weights: pd.DataFrame | None = None,
    ) -> KinaseScoringResult:
        return cls(
            profile_scores=profile_scores,
            motif_scores=motif_scores,
            combined_scores=combined_scores,
            weights=weights,
            _assume_owned=True,
        )


@dataclass(frozen=True, slots=True)
class KinasePredictionResult:
    """Prediction-stage outputs."""

    pred_mat: pd.DataFrame
    substrate_list: pd.DataFrame | None = None
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        pred_mat = own_dataframe(
            self.pred_mat,
            field_name="prediction_result.pred_mat",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
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
