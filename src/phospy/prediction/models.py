"""Prediction and scoring stage result models."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
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
from phospy.prediction.scoring import (
    KINASE_SCORE_SOURCE_SUMMARY_COLUMNS,
    KINASE_SCORE_SOURCE_VALUES,
)
from phospy.prediction.sequence_validation import SequenceValidationResult
from phospy.tables.base import require_canonical_label_index
from phospy.tables.kinase import KinasePredictionMatrix, KinaseScoreMatrix
from phospy.validation.common.dataframes import (
    require_canonical_site_index,
    require_dataframe,
    require_exact_index_match,
    require_finite_numeric_dataframe,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
)


@dataclass(frozen=True, slots=True, init=False)
class KinaseScoringResult:
    """Scoring-stage outputs.

    `profile_scores` and `rank_weighted_fusion_scores` define the supported
    downstream lane. `motif_scores` and `score_fusion_weights` are optional
    diagnostic tables controlled by
    `scoring_config.include_diagnostic_scoring_tables`.
    `score_source_summary` is a compact per-kinase evidence-source diagnostic.
    `score_source_matrix` is an optional per-site/per-kinase evidence-source
    diagnostic table.
    """

    motif_sequence_validation: SequenceValidationResult | None = None
    motif_library_validation: MotifLibraryValidationResult | None = None
    _profile_scores: pd.DataFrame = field(init=False, repr=False)
    _motif_scores: pd.DataFrame | None = field(init=False, repr=False)
    _rank_weighted_fusion_scores: pd.DataFrame | None = field(init=False, repr=False)
    _score_fusion_weights: pd.DataFrame | None = field(init=False, repr=False)
    _score_source_matrix: pd.DataFrame | None = field(init=False, repr=False)
    _score_source_summary: pd.DataFrame | None = field(init=False, repr=False)

    def __init__(
        self,
        profile_scores: pd.DataFrame,
        motif_scores: pd.DataFrame | None = None,
        rank_weighted_fusion_scores: pd.DataFrame | None = None,
        score_fusion_weights: pd.DataFrame | None = None,
        score_source_matrix: pd.DataFrame | None = None,
        score_source_summary: pd.DataFrame | None = None,
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
        score_source_matrix = own_optional_dataframe(
            score_source_matrix,
            field_name="scoring_result.score_source_matrix",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        score_source_summary = own_optional_dataframe(
            score_source_summary,
            field_name="scoring_result.score_source_summary",
            error_type=PhosPyValidationError,
            assume_owned=_assume_owned,
        )
        score_source_matrix = _validate_score_source_matrix(
            score_source_matrix=score_source_matrix,
            profile_scores=profile_scores,
            rank_weighted_fusion_scores=rank_weighted_fusion_scores,
        )
        score_source_summary = _validate_score_source_summary(
            score_source_summary=score_source_summary,
            profile_scores=profile_scores,
            rank_weighted_fusion_scores=rank_weighted_fusion_scores,
        )
        object.__setattr__(self, "_profile_scores", profile_scores)
        object.__setattr__(self, "_motif_scores", motif_scores)
        object.__setattr__(
            self, "_rank_weighted_fusion_scores", rank_weighted_fusion_scores
        )
        object.__setattr__(self, "_score_fusion_weights", score_fusion_weights)
        object.__setattr__(self, "_score_source_matrix", score_source_matrix)
        object.__setattr__(self, "_score_source_summary", score_source_summary)
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

    @property
    def score_source_matrix(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._score_source_matrix)

    @property
    def score_source_summary(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._score_source_summary)

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

    def _borrow_score_source_matrix_frame(self) -> pd.DataFrame | None:
        """Package-private borrowed score-source matrix for internal workflows."""

        return _borrow_optional_dataframe(self._score_source_matrix)

    def _borrow_score_source_summary_frame(self) -> pd.DataFrame | None:
        """Package-private borrowed score-source summary for internal workflows."""

        return _borrow_optional_dataframe(self._score_source_summary)

    @classmethod
    def _from_owned(
        cls,
        *,
        profile_scores: pd.DataFrame,
        motif_scores: pd.DataFrame | None = None,
        rank_weighted_fusion_scores: pd.DataFrame | None = None,
        score_fusion_weights: pd.DataFrame | None = None,
        score_source_matrix: pd.DataFrame | None = None,
        score_source_summary: pd.DataFrame | None = None,
        motif_sequence_validation: SequenceValidationResult | None = None,
        motif_library_validation: MotifLibraryValidationResult | None = None,
    ) -> KinaseScoringResult:
        return cls(
            profile_scores=profile_scores,
            motif_scores=motif_scores,
            rank_weighted_fusion_scores=rank_weighted_fusion_scores,
            score_fusion_weights=score_fusion_weights,
            score_source_matrix=score_source_matrix,
            score_source_summary=score_source_summary,
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

    def score_source_matrix_dataframe(self) -> pd.DataFrame | None:
        """Return an optional score-source matrix snapshot isolated from this result."""

        return export_optional_dataframe(self._score_source_matrix)

    def score_source_summary_dataframe(self) -> pd.DataFrame | None:
        """Return an optional score-source summary snapshot isolated from this result."""

        return export_optional_dataframe(self._score_source_summary)


def _validate_score_source_matrix(
    *,
    score_source_matrix: pd.DataFrame | None,
    profile_scores: pd.DataFrame,
    rank_weighted_fusion_scores: pd.DataFrame | None,
) -> pd.DataFrame | None:
    if score_source_matrix is None:
        return None
    require_dataframe(
        score_source_matrix,
        field_name="scoring_result.score_source_matrix",
        allow_empty=False,
        error_type=PhosPyValidationError,
    )
    require_unique_columns(
        score_source_matrix,
        field_name="scoring_result.score_source_matrix",
        error_type=PhosPyValidationError,
    )
    require_unique_index(
        score_source_matrix,
        field_name="scoring_result.score_source_matrix",
        error_type=PhosPyValidationError,
    )
    require_canonical_label_index(
        score_source_matrix.columns,
        field_name="scoring_result.score_source_matrix.columns",
        error_type=PhosPyValidationError,
    )
    require_canonical_site_index(
        score_source_matrix.index,
        field_name="scoring_result.score_source_matrix.index",
        error_type=PhosPyValidationError,
    )
    expected = (
        rank_weighted_fusion_scores
        if rank_weighted_fusion_scores is not None
        else profile_scores
    )
    require_exact_index_match(
        left=score_source_matrix.index,
        right=expected.index,
        left_name="scoring_result.score_source_matrix.index",
        right_name=(
            "scoring_result.rank_weighted_fusion_scores.index"
            if rank_weighted_fusion_scores is not None
            else "scoring_result.profile_scores.index"
        ),
        error_type=PhosPyValidationError,
    )
    require_exact_index_match(
        left=score_source_matrix.columns,
        right=expected.columns,
        left_name="scoring_result.score_source_matrix.columns",
        right_name=(
            "scoring_result.rank_weighted_fusion_scores.columns"
            if rank_weighted_fusion_scores is not None
            else "scoring_result.profile_scores.columns"
        ),
        error_type=PhosPyValidationError,
    )
    raw_values = score_source_matrix.to_numpy(dtype=object, copy=False).ravel()
    invalid_values = sorted(
        {
            str(value)
            for value in raw_values
            if not isinstance(value, str) or value not in KINASE_SCORE_SOURCE_VALUES
        }
    )
    if invalid_values:
        preview = ", ".join(invalid_values[:5])
        suffix = "" if len(invalid_values) <= 5 else " ..."
        raise PhosPyValidationError(
            "scoring_result.score_source_matrix contains unsupported source labels: "
            f"{preview}{suffix}"
        )
    return score_source_matrix


def _validate_score_source_summary(
    *,
    score_source_summary: pd.DataFrame | None,
    profile_scores: pd.DataFrame,
    rank_weighted_fusion_scores: pd.DataFrame | None,
) -> pd.DataFrame | None:
    if score_source_summary is None:
        return None
    require_dataframe(
        score_source_summary,
        field_name="scoring_result.score_source_summary",
        allow_empty=False,
        error_type=PhosPyValidationError,
    )
    require_unique_columns(
        score_source_summary,
        field_name="scoring_result.score_source_summary",
        error_type=PhosPyValidationError,
    )
    require_unique_index(
        score_source_summary,
        field_name="scoring_result.score_source_summary",
        error_type=PhosPyValidationError,
    )
    require_numeric_dataframe(
        score_source_summary,
        field_name="scoring_result.score_source_summary",
        error_type=PhosPyValidationError,
    )
    require_finite_numeric_dataframe(
        score_source_summary,
        field_name="scoring_result.score_source_summary",
        error_type=PhosPyValidationError,
        allow_missing=False,
    )
    require_canonical_label_index(
        score_source_summary.index,
        field_name="scoring_result.score_source_summary.index",
        error_type=PhosPyValidationError,
    )
    require_canonical_label_index(
        score_source_summary.columns,
        field_name="scoring_result.score_source_summary.columns",
        error_type=PhosPyValidationError,
    )
    observed_columns = tuple(score_source_summary.columns.astype(str))
    if observed_columns != KINASE_SCORE_SOURCE_SUMMARY_COLUMNS:
        raise PhosPyValidationError(
            "scoring_result.score_source_summary columns must match "
            f"{list(KINASE_SCORE_SOURCE_SUMMARY_COLUMNS)}"
        )
    expected_index = (
        rank_weighted_fusion_scores.columns
        if rank_weighted_fusion_scores is not None
        else profile_scores.columns
    )
    require_exact_index_match(
        left=score_source_summary.index,
        right=expected_index,
        left_name="scoring_result.score_source_summary.index",
        right_name=(
            "scoring_result.rank_weighted_fusion_scores.columns"
            if rank_weighted_fusion_scores is not None
            else "scoring_result.profile_scores.columns"
        ),
        error_type=PhosPyValidationError,
    )
    if (score_source_summary.to_numpy(dtype=float) < 0.0).any():
        raise PhosPyValidationError(
            "scoring_result.score_source_summary must contain non-negative counts"
        )
    counts = score_source_summary.loc[
        :,
        [
            "fused_motif_profile_evidence_count",
            "profile_only_motif_missing_or_constant_count",
            "profile_only_no_motif_overlap_count",
            "unavailable_no_score_count",
        ],
    ].to_numpy(dtype=float)
    if not np.allclose(counts, np.floor(counts)):
        raise PhosPyValidationError(
            "scoring_result.score_source_summary must contain integer count values"
        )
    component_counts = score_source_summary.loc[
        :,
        [
            "fused_motif_profile_evidence_count",
            "profile_only_motif_missing_or_constant_count",
            "profile_only_no_motif_overlap_count",
        ],
    ].sum(axis=1)
    total_counts = score_source_summary.loc[:, "total_sites_count"]
    unavailable_counts = score_source_summary.loc[:, "unavailable_no_score_count"]
    if not np.allclose(
        component_counts.to_numpy(dtype=float)
        + unavailable_counts.to_numpy(dtype=float),
        total_counts.to_numpy(dtype=float),
    ):
        raise PhosPyValidationError(
            "scoring_result.score_source_summary component counts must sum to total_sites_count"
        )
    with_score_counts = score_source_summary.loc[:, "sites_with_score_count"]
    if not np.allclose(
        component_counts.to_numpy(dtype=float),
        with_score_counts.to_numpy(dtype=float),
    ):
        raise PhosPyValidationError(
            "scoring_result.score_source_summary sites_with_score_count must equal "
            "fused + profile-only counts"
        )
    return score_source_summary


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
