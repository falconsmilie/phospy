"""Prediction and scoring stage result models."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from phospy.errors.validation import PhosPyValidationError
from phospy.frames.ownership import (
    _borrow_dataframe,
    _borrow_optional_dataframe,
    export_dataframe,
    export_optional_dataframe,
    own_optional_dataframe,
)
from phospy.frames.validation import (
    require_dataframe,
    require_exact_index_match,
    require_finite_numeric_dataframe,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
)
from phospy.policies import coerce_policy_enum
from phospy.provenance.immutability import freeze_json_mapping_with_error_type
from phospy.science.configs import (
    KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
    KINASE_SCORING_MODES,
    normalize_kinase_scoring_mode,
)
from phospy.science.prediction.motif_scoring import (
    KinaseLibraryMotifScoringResult,
    MotifLibraryValidationResult,
)
from phospy.science.prediction.scoring import (
    KINASE_SCORE_SOURCE_SUMMARY_COLUMNS,
    KINASE_SCORE_SOURCE_VALUES,
)
from phospy.science.prediction.sequence_validation import SequenceValidationResult
from phospy.science.scoring.policy_models import (
    DownstreamScoreSource,
    ProfileSelfInclusionPolicy,
)
from phospy.science.sites.validation import require_site_key_index
from phospy.science.tables.base import require_canonical_label_index
from phospy.science.tables.kinase import (
    KinasePredictionMatrix,
    KinaseProfileScoreDiagnosticTable,
    KinaseScoreMatrix,
)


@dataclass(frozen=True, slots=True, init=False)
class KinaseScoringResult:
    """Scoring-stage outputs.

    `score_source` identifies the authoritative downstream score matrix used by
    prediction. The default remains the historical PhosR-inspired rank-weighted
    lane when it is present, otherwise profile scores. This lane is PhosPy
    scoring, not an exact PhosR implementation or numerical compatibility mode.

    `profile_scores` and `rank_weighted_fusion_scores` define the historical
    PhosPy downstream lane. `motif_scores` and `score_fusion_weights` are
    optional diagnostic tables controlled by
    `scoring_config.include_diagnostic_scoring_tables`.
    `kinase_library_motif_scores` is populated when Kinase Library motif scoring
    is explicitly selected.
    `score_source_summary` is a compact per-kinase evidence-source diagnostic.
    `score_source_matrix` is an optional per-site/per-kinase evidence-source
    diagnostic table.
    `profile_score_diagnostics` records sparse profile-scoring diagnostics such
    as leave-one-out cells that could not be scored after self-exclusion.
    `profile_self_inclusion_policy` records whether known substrate sites were
    allowed to contribute to their own kinase profile scores.
    """

    motif_sequence_validation: SequenceValidationResult | None = None
    motif_library_validation: MotifLibraryValidationResult | None = None
    scoring_mode: str = KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED
    score_scale: str = "relative_support_score_unit_interval"
    score_scale_metadata: Mapping[str, object] | None = None
    profile_self_inclusion_policy: ProfileSelfInclusionPolicy = (
        ProfileSelfInclusionPolicy.ALLOW
    )
    _profile_scores: pd.DataFrame = field(init=False, repr=False)
    _motif_scores: pd.DataFrame | None = field(init=False, repr=False)
    _rank_weighted_fusion_scores: pd.DataFrame | None = field(init=False, repr=False)
    _kinase_library_motif_scores: pd.DataFrame | None = field(
        init=False,
        repr=False,
    )
    _combined_profile_motif_scores: pd.DataFrame | None = field(
        init=False,
        repr=False,
    )
    _score_fusion_weights: pd.DataFrame | None = field(init=False, repr=False)
    _score_source_matrix: pd.DataFrame | None = field(init=False, repr=False)
    _score_source_summary: pd.DataFrame | None = field(init=False, repr=False)
    _profile_score_diagnostics: pd.DataFrame | None = field(
        init=False,
        repr=False,
    )
    _kinase_library_site_diagnostics: pd.DataFrame | None = field(
        init=False,
        repr=False,
    )
    _kinase_library_kinase_diagnostics: pd.DataFrame | None = field(
        init=False,
        repr=False,
    )
    _score_source: DownstreamScoreSource = field(init=False, repr=False)

    def __init__(
        self,
        profile_scores: pd.DataFrame,
        motif_scores: pd.DataFrame | None = None,
        rank_weighted_fusion_scores: pd.DataFrame | None = None,
        kinase_library_motif_scores: pd.DataFrame | None = None,
        combined_profile_motif_scores: pd.DataFrame | None = None,
        score_fusion_weights: pd.DataFrame | None = None,
        score_source_matrix: pd.DataFrame | None = None,
        score_source_summary: pd.DataFrame | None = None,
        profile_score_diagnostics: pd.DataFrame | None = None,
        kinase_library_site_diagnostics: pd.DataFrame | None = None,
        kinase_library_kinase_diagnostics: pd.DataFrame | None = None,
        scoring_mode: str = KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
        score_source: DownstreamScoreSource | str | None = None,
        score_scale: str | None = None,
        score_scale_metadata: Mapping[str, object] | None = None,
        profile_self_inclusion_policy: ProfileSelfInclusionPolicy | str = (
            ProfileSelfInclusionPolicy.ALLOW
        ),
        motif_sequence_validation: SequenceValidationResult | None = None,
        motif_library_validation: MotifLibraryValidationResult | None = None,
    ) -> None:
        object.__setattr__(self, "motif_sequence_validation", motif_sequence_validation)
        object.__setattr__(self, "motif_library_validation", motif_library_validation)
        resolved_scoring_mode = _validate_scoring_mode(scoring_mode)
        resolved_score_source = _resolve_score_source(
            score_source=score_source,
            rank_weighted_fusion_scores=rank_weighted_fusion_scores,
            kinase_library_motif_scores=kinase_library_motif_scores,
            combined_profile_motif_scores=combined_profile_motif_scores,
        )
        resolved_score_scale = _resolve_score_scale(
            score_scale=score_scale,
            score_source=resolved_score_source,
        )
        object.__setattr__(self, "scoring_mode", resolved_scoring_mode)
        object.__setattr__(self, "score_scale", resolved_score_scale)
        object.__setattr__(
            self,
            "score_scale_metadata",
            _own_score_scale_metadata(score_scale_metadata),
        )
        object.__setattr__(
            self,
            "profile_self_inclusion_policy",
            _validate_profile_self_inclusion_policy(profile_self_inclusion_policy),
        )
        object.__setattr__(self, "_score_source", resolved_score_source)
        profile_scores = KinaseScoreMatrix(
            frame=profile_scores,
            field_name="scoring_result.profile_scores",
            _assume_owned=False,
        ).frame
        motif_scores = (
            None
            if motif_scores is None
            else KinaseScoreMatrix(
                frame=motif_scores,
                field_name="scoring_result.motif_scores",
                _assume_owned=False,
            ).frame
        )
        rank_weighted_fusion_scores = (
            None
            if rank_weighted_fusion_scores is None
            else KinaseScoreMatrix(
                frame=rank_weighted_fusion_scores,
                field_name="scoring_result.rank_weighted_fusion_scores",
                _assume_owned=False,
            ).frame
        )
        kinase_library_motif_scores = (
            None
            if kinase_library_motif_scores is None
            else KinaseScoreMatrix(
                frame=kinase_library_motif_scores,
                field_name="scoring_result.kinase_library_motif_scores",
                enforce_unit_interval=True,
                _assume_owned=False,
            ).frame
        )
        combined_profile_motif_scores = (
            None
            if combined_profile_motif_scores is None
            else KinaseScoreMatrix(
                frame=combined_profile_motif_scores,
                field_name="scoring_result.combined_profile_motif_scores",
                enforce_unit_interval=True,
                _assume_owned=False,
            ).frame
        )
        score_fusion_weights = (
            None
            if score_fusion_weights is None
            else KinaseScoreMatrix(
                frame=score_fusion_weights,
                field_name="scoring_result.score_fusion_weights",
                require_site_index=False,
                _assume_owned=False,
            ).frame
        )
        score_source_matrix = own_optional_dataframe(
            score_source_matrix,
            field_name="scoring_result.score_source_matrix",
            error_type=PhosPyValidationError,
            assume_owned=False,
        )
        score_source_summary = own_optional_dataframe(
            score_source_summary,
            field_name="scoring_result.score_source_summary",
            error_type=PhosPyValidationError,
            assume_owned=False,
        )
        profile_score_diagnostics = _own_optional_profile_score_diagnostics(
            profile_score_diagnostics,
            assume_owned=False,
        )
        kinase_library_site_diagnostics = own_optional_dataframe(
            kinase_library_site_diagnostics,
            field_name="scoring_result.kinase_library_site_diagnostics",
            error_type=PhosPyValidationError,
            assume_owned=False,
        )
        kinase_library_kinase_diagnostics = own_optional_dataframe(
            kinase_library_kinase_diagnostics,
            field_name="scoring_result.kinase_library_kinase_diagnostics",
            error_type=PhosPyValidationError,
            assume_owned=False,
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
        _validate_authoritative_score_source(
            score_source=resolved_score_source,
            profile_scores=profile_scores,
            rank_weighted_fusion_scores=rank_weighted_fusion_scores,
            kinase_library_motif_scores=kinase_library_motif_scores,
            combined_profile_motif_scores=combined_profile_motif_scores,
        )
        object.__setattr__(self, "_profile_scores", profile_scores)
        object.__setattr__(self, "_motif_scores", motif_scores)
        object.__setattr__(
            self, "_rank_weighted_fusion_scores", rank_weighted_fusion_scores
        )
        object.__setattr__(
            self, "_kinase_library_motif_scores", kinase_library_motif_scores
        )
        object.__setattr__(
            self, "_combined_profile_motif_scores", combined_profile_motif_scores
        )
        object.__setattr__(self, "_score_fusion_weights", score_fusion_weights)
        object.__setattr__(self, "_score_source_matrix", score_source_matrix)
        object.__setattr__(self, "_score_source_summary", score_source_summary)
        object.__setattr__(
            self,
            "_profile_score_diagnostics",
            profile_score_diagnostics,
        )
        object.__setattr__(
            self,
            "_kinase_library_site_diagnostics",
            kinase_library_site_diagnostics,
        )
        object.__setattr__(
            self,
            "_kinase_library_kinase_diagnostics",
            kinase_library_kinase_diagnostics,
        )
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
    def kinase_library_motif_scores(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._kinase_library_motif_scores)

    @property
    def combined_profile_motif_scores(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._combined_profile_motif_scores)

    @property
    def score_fusion_weights(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._score_fusion_weights)

    @property
    def score_source_matrix(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._score_source_matrix)

    @property
    def score_source_summary(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._score_source_summary)

    @property
    def profile_score_diagnostics(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._profile_score_diagnostics)

    @property
    def kinase_library_site_diagnostics(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._kinase_library_site_diagnostics)

    @property
    def kinase_library_kinase_diagnostics(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._kinase_library_kinase_diagnostics)

    @property
    def score_source(self) -> str:
        return self._score_source.value

    @property
    def authoritative_scores(self) -> pd.DataFrame:
        return export_dataframe(self._borrow_authoritative_scores_frame())

    def _borrow_profile_scores_frame(self) -> pd.DataFrame:
        """Package-private profile-score snapshot for internal views."""

        return _borrow_dataframe(self._profile_scores)

    def _borrow_motif_scores_frame(self) -> pd.DataFrame | None:
        """Package-private motif-score snapshot for internal views."""

        return _borrow_optional_dataframe(self._motif_scores)

    def _borrow_rank_weighted_fusion_scores_frame(self) -> pd.DataFrame | None:
        """Package-private fusion-score snapshot for internal views."""

        return _borrow_optional_dataframe(self._rank_weighted_fusion_scores)

    def _borrow_kinase_library_motif_scores_frame(self) -> pd.DataFrame | None:
        """Package-private Kinase Library score snapshot for internal views."""

        return _borrow_optional_dataframe(self._kinase_library_motif_scores)

    def _borrow_combined_profile_motif_scores_frame(self) -> pd.DataFrame | None:
        """Package-private combined profile/motif score snapshot."""

        return _borrow_optional_dataframe(self._combined_profile_motif_scores)

    def _borrow_score_fusion_weights_frame(self) -> pd.DataFrame | None:
        """Package-private fusion-weight snapshot for internal views."""

        return _borrow_optional_dataframe(self._score_fusion_weights)

    def _borrow_score_source_matrix_frame(self) -> pd.DataFrame | None:
        """Package-private score-source matrix snapshot for internal views."""

        return _borrow_optional_dataframe(self._score_source_matrix)

    def _borrow_score_source_summary_frame(self) -> pd.DataFrame | None:
        """Package-private score-source summary snapshot for internal views."""

        return _borrow_optional_dataframe(self._score_source_summary)

    def _borrow_profile_score_diagnostics_frame(self) -> pd.DataFrame | None:
        """Package-private profile-score diagnostic snapshot."""

        return _borrow_optional_dataframe(self._profile_score_diagnostics)

    def _borrow_kinase_library_site_diagnostics_frame(self) -> pd.DataFrame | None:
        """Package-private Kinase Library site-diagnostics snapshot."""

        return _borrow_optional_dataframe(self._kinase_library_site_diagnostics)

    def _borrow_kinase_library_kinase_diagnostics_frame(self) -> pd.DataFrame | None:
        """Package-private Kinase Library kinase-diagnostics snapshot."""

        return _borrow_optional_dataframe(self._kinase_library_kinase_diagnostics)

    def _borrow_authoritative_scores_frame(self) -> pd.DataFrame:
        """Package-private authoritative downstream score snapshot."""

        if self._score_source is DownstreamScoreSource.RANK_WEIGHTED_FUSION_SCORES:
            if self._rank_weighted_fusion_scores is not None:
                return _borrow_dataframe(self._rank_weighted_fusion_scores)
        if self._score_source is DownstreamScoreSource.KINASE_LIBRARY_MOTIF_SCORES:
            if self._kinase_library_motif_scores is not None:
                return _borrow_dataframe(self._kinase_library_motif_scores)
        if self._score_source is DownstreamScoreSource.COMBINED_PROFILE_MOTIF_SCORES:
            if self._combined_profile_motif_scores is not None:
                return _borrow_dataframe(self._combined_profile_motif_scores)
        return _borrow_dataframe(self._profile_scores)

    @classmethod
    def _from_owned(
        cls,
        *,
        profile_scores: pd.DataFrame,
        motif_scores: pd.DataFrame | None = None,
        rank_weighted_fusion_scores: pd.DataFrame | None = None,
        kinase_library_motif_scores: pd.DataFrame | None = None,
        combined_profile_motif_scores: pd.DataFrame | None = None,
        score_fusion_weights: pd.DataFrame | None = None,
        score_source_matrix: pd.DataFrame | None = None,
        score_source_summary: pd.DataFrame | None = None,
        profile_score_diagnostics: pd.DataFrame | None = None,
        kinase_library_site_diagnostics: pd.DataFrame | None = None,
        kinase_library_kinase_diagnostics: pd.DataFrame | None = None,
        scoring_mode: str = KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
        score_source: DownstreamScoreSource | str | None = None,
        score_scale: str | None = None,
        score_scale_metadata: Mapping[str, object] | None = None,
        profile_self_inclusion_policy: ProfileSelfInclusionPolicy | str = (
            ProfileSelfInclusionPolicy.ALLOW
        ),
        motif_sequence_validation: SequenceValidationResult | None = None,
        motif_library_validation: MotifLibraryValidationResult | None = None,
    ) -> KinaseScoringResult:
        _require_frame_type(
            profile_scores,
            field_name="scoring_result.profile_scores",
        )
        _require_optional_frame_type(
            motif_scores,
            field_name="scoring_result.motif_scores",
        )
        _require_optional_frame_type(
            rank_weighted_fusion_scores,
            field_name="scoring_result.rank_weighted_fusion_scores",
        )
        _require_optional_frame_type(
            kinase_library_motif_scores,
            field_name="scoring_result.kinase_library_motif_scores",
        )
        _require_optional_frame_type(
            combined_profile_motif_scores,
            field_name="scoring_result.combined_profile_motif_scores",
        )
        _require_optional_frame_type(
            score_fusion_weights,
            field_name="scoring_result.score_fusion_weights",
        )
        _require_optional_frame_type(
            score_source_matrix,
            field_name="scoring_result.score_source_matrix",
        )
        _require_optional_frame_type(
            score_source_summary,
            field_name="scoring_result.score_source_summary",
        )
        _require_optional_frame_type(
            profile_score_diagnostics,
            field_name="scoring_result.profile_score_diagnostics",
        )
        _require_optional_frame_type(
            kinase_library_site_diagnostics,
            field_name="scoring_result.kinase_library_site_diagnostics",
        )
        _require_optional_frame_type(
            kinase_library_kinase_diagnostics,
            field_name="scoring_result.kinase_library_kinase_diagnostics",
        )
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

        result = object.__new__(cls)
        resolved_score_source = _resolve_score_source(
            score_source=score_source,
            rank_weighted_fusion_scores=rank_weighted_fusion_scores,
            kinase_library_motif_scores=kinase_library_motif_scores,
            combined_profile_motif_scores=combined_profile_motif_scores,
        )
        _validate_authoritative_score_source(
            score_source=resolved_score_source,
            profile_scores=profile_scores,
            rank_weighted_fusion_scores=rank_weighted_fusion_scores,
            kinase_library_motif_scores=kinase_library_motif_scores,
            combined_profile_motif_scores=combined_profile_motif_scores,
        )
        object.__setattr__(
            result, "motif_sequence_validation", motif_sequence_validation
        )
        object.__setattr__(result, "motif_library_validation", motif_library_validation)
        object.__setattr__(result, "scoring_mode", _validate_scoring_mode(scoring_mode))
        object.__setattr__(
            result,
            "score_scale",
            _resolve_score_scale(
                score_scale=score_scale,
                score_source=resolved_score_source,
            ),
        )
        object.__setattr__(
            result,
            "score_scale_metadata",
            _own_score_scale_metadata(score_scale_metadata),
        )
        object.__setattr__(
            result,
            "profile_self_inclusion_policy",
            _validate_profile_self_inclusion_policy(profile_self_inclusion_policy),
        )
        object.__setattr__(result, "_score_source", resolved_score_source)
        object.__setattr__(result, "_profile_scores", profile_scores)
        object.__setattr__(result, "_motif_scores", motif_scores)
        object.__setattr__(
            result,
            "_rank_weighted_fusion_scores",
            rank_weighted_fusion_scores,
        )
        object.__setattr__(
            result,
            "_kinase_library_motif_scores",
            kinase_library_motif_scores,
        )
        object.__setattr__(
            result,
            "_combined_profile_motif_scores",
            combined_profile_motif_scores,
        )
        object.__setattr__(result, "_score_fusion_weights", score_fusion_weights)
        object.__setattr__(result, "_score_source_matrix", score_source_matrix)
        object.__setattr__(result, "_score_source_summary", score_source_summary)
        object.__setattr__(
            result,
            "_profile_score_diagnostics",
            profile_score_diagnostics,
        )
        object.__setattr__(
            result,
            "_kinase_library_site_diagnostics",
            kinase_library_site_diagnostics,
        )
        object.__setattr__(
            result,
            "_kinase_library_kinase_diagnostics",
            kinase_library_kinase_diagnostics,
        )
        return result

    def to_dataframe(self) -> pd.DataFrame:
        """Return a `profile_scores` snapshot isolated from this result."""

        return export_dataframe(self._profile_scores)

    def motif_scores_dataframe(self) -> pd.DataFrame | None:
        """Return an optional motif-score snapshot isolated from this result."""

        return export_optional_dataframe(self._motif_scores)

    def rank_weighted_fusion_scores_dataframe(self) -> pd.DataFrame | None:
        """Return an optional fusion-score snapshot isolated from this result."""

        return export_optional_dataframe(self._rank_weighted_fusion_scores)

    def kinase_library_motif_scores_dataframe(self) -> pd.DataFrame | None:
        """Return optional Kinase Library motif scores isolated from this result."""

        return export_optional_dataframe(self._kinase_library_motif_scores)

    def combined_profile_motif_scores_dataframe(self) -> pd.DataFrame | None:
        """Return optional combined profile/motif scores isolated from this result."""

        return export_optional_dataframe(self._combined_profile_motif_scores)

    def score_fusion_weights_dataframe(self) -> pd.DataFrame | None:
        """Return an optional fusion-weight snapshot isolated from this result."""

        return export_optional_dataframe(self._score_fusion_weights)

    def score_source_matrix_dataframe(self) -> pd.DataFrame | None:
        """Return an optional score-source matrix snapshot isolated from this result."""

        return export_optional_dataframe(self._score_source_matrix)

    def score_source_summary_dataframe(self) -> pd.DataFrame | None:
        """Return an optional score-source summary snapshot isolated from this result."""

        return export_optional_dataframe(self._score_source_summary)

    def profile_score_diagnostics_dataframe(self) -> pd.DataFrame | None:
        """Return optional sparse profile-score diagnostics."""

        return export_optional_dataframe(self._profile_score_diagnostics)

    def kinase_library_site_diagnostics_dataframe(self) -> pd.DataFrame | None:
        """Return optional Kinase Library site diagnostics."""

        return export_optional_dataframe(self._kinase_library_site_diagnostics)

    def kinase_library_kinase_diagnostics_dataframe(self) -> pd.DataFrame | None:
        """Return optional Kinase Library matrix diagnostics."""

        return export_optional_dataframe(self._kinase_library_kinase_diagnostics)

    def authoritative_scores_dataframe(self) -> pd.DataFrame:
        """Return the selected downstream score matrix."""

        return export_dataframe(self._borrow_authoritative_scores_frame())


def _validate_scoring_mode(value: object) -> str:
    text = normalize_kinase_scoring_mode(value, warn_on_deprecated_alias=True)
    if text not in KINASE_SCORING_MODES:
        allowed = ", ".join(sorted(KINASE_SCORING_MODES))
        raise PhosPyValidationError(
            f"scoring_result.scoring_mode must be one of: {allowed}"
        )
    return text


def _validate_profile_self_inclusion_policy(
    value: object,
) -> ProfileSelfInclusionPolicy:
    return coerce_policy_enum(
        ProfileSelfInclusionPolicy,
        value,
        field_name="scoring_result.profile_self_inclusion_policy",
        error_type=PhosPyValidationError,
    )


def _resolve_score_source(
    *,
    score_source: DownstreamScoreSource | str | None,
    rank_weighted_fusion_scores: pd.DataFrame | None,
    kinase_library_motif_scores: pd.DataFrame | None,
    combined_profile_motif_scores: pd.DataFrame | None,
) -> DownstreamScoreSource:
    if score_source is not None:
        return DownstreamScoreSource.parse(
            score_source,
            field_name="scoring_result.score_source",
        )
    if combined_profile_motif_scores is not None:
        return DownstreamScoreSource.COMBINED_PROFILE_MOTIF_SCORES
    if kinase_library_motif_scores is not None and rank_weighted_fusion_scores is None:
        return DownstreamScoreSource.KINASE_LIBRARY_MOTIF_SCORES
    if rank_weighted_fusion_scores is not None:
        return DownstreamScoreSource.RANK_WEIGHTED_FUSION_SCORES
    return DownstreamScoreSource.PROFILE_SCORES


def _resolve_score_scale(
    *,
    score_scale: str | None,
    score_source: DownstreamScoreSource,
) -> str:
    if score_scale is not None:
        resolved = str(score_scale).strip()
        if resolved:
            return resolved
        raise PhosPyValidationError("scoring_result.score_scale must be non-empty")
    if score_source is DownstreamScoreSource.KINASE_LIBRARY_MOTIF_SCORES:
        return "kinase_library_motif_minmax_unit_interval"
    if score_source is DownstreamScoreSource.COMBINED_PROFILE_MOTIF_SCORES:
        return "combined_profile_kinase_library_motif_unit_interval"
    return "relative_support_score_unit_interval"


def _own_score_scale_metadata(
    value: Mapping[str, object] | None,
) -> Mapping[str, object] | None:
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise PhosPyValidationError(
            "scoring_result.score_scale_metadata must be a mapping or None"
        )
    return freeze_json_mapping_with_error_type(
        value,
        field_name="scoring_result.score_scale_metadata",
        error_type=PhosPyValidationError,
    )


def _own_optional_profile_score_diagnostics(
    table: pd.DataFrame | None,
    *,
    assume_owned: bool,
) -> pd.DataFrame | None:
    if table is None:
        return None
    return KinaseProfileScoreDiagnosticTable(
        frame=table,
        _assume_owned=assume_owned,
    ).frame


def _validate_authoritative_score_source(
    *,
    score_source: DownstreamScoreSource,
    profile_scores: pd.DataFrame,
    rank_weighted_fusion_scores: pd.DataFrame | None,
    kinase_library_motif_scores: pd.DataFrame | None,
    combined_profile_motif_scores: pd.DataFrame | None,
) -> None:
    if score_source is DownstreamScoreSource.PROFILE_SCORES:
        return
    if score_source is DownstreamScoreSource.RANK_WEIGHTED_FUSION_SCORES:
        _require_authoritative_matrix(
            rank_weighted_fusion_scores,
            profile_scores=profile_scores,
            field_name="scoring_result.rank_weighted_fusion_scores",
        )
        return
    if score_source is DownstreamScoreSource.KINASE_LIBRARY_MOTIF_SCORES:
        _require_authoritative_matrix(
            kinase_library_motif_scores,
            profile_scores=profile_scores,
            field_name="scoring_result.kinase_library_motif_scores",
        )
        return
    if score_source is DownstreamScoreSource.COMBINED_PROFILE_MOTIF_SCORES:
        _require_authoritative_matrix(
            combined_profile_motif_scores,
            profile_scores=profile_scores,
            field_name="scoring_result.combined_profile_motif_scores",
        )
        return
    raise PhosPyValidationError(
        f"scoring_result.score_source is unsupported: {score_source.value}"
    )


def _require_authoritative_matrix(
    matrix: pd.DataFrame | None,
    *,
    profile_scores: pd.DataFrame,
    field_name: str,
) -> None:
    if matrix is None:
        raise PhosPyValidationError(
            f"{field_name} is required by scoring_result.score_source"
        )
    require_exact_index_match(
        left=matrix.index,
        right=profile_scores.index,
        left_name=f"{field_name}.index",
        right_name="scoring_result.profile_scores.index",
        error_type=PhosPyValidationError,
    )


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
    require_site_key_index(
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
    ) -> None:
        pred_mat = KinasePredictionMatrix(
            frame=pred_mat,
            field_name="prediction_result.pred_mat",
            _assume_owned=False,
        ).frame
        substrate_list = own_optional_dataframe(
            substrate_list,
            field_name="prediction_result.substrate_list",
            error_type=PhosPyValidationError,
            assume_owned=False,
        )
        object.__setattr__(self, "_pred_mat", pred_mat)
        object.__setattr__(self, "_substrate_list", substrate_list)

    @property
    def pred_mat(self) -> pd.DataFrame:
        return export_dataframe(self._pred_mat)

    @property
    def substrate_list(self) -> pd.DataFrame | None:
        return _export_public_substrate_list(self._substrate_list)

    def _borrow_pred_mat_frame(self) -> pd.DataFrame:
        """Package-private prediction-matrix snapshot for internal views."""

        return _borrow_dataframe(self._pred_mat)

    def _borrow_substrate_list_frame(self) -> pd.DataFrame | None:
        """Package-private substrate-list snapshot for internal views."""

        return _borrow_optional_dataframe(self._substrate_list)

    @classmethod
    def _from_owned(
        cls,
        *,
        pred_mat: pd.DataFrame,
        substrate_list: pd.DataFrame | None = None,
    ) -> KinasePredictionResult:
        _require_frame_type(pred_mat, field_name="prediction_result.pred_mat")
        _require_optional_frame_type(
            substrate_list,
            field_name="prediction_result.substrate_list",
        )
        result = object.__new__(cls)
        object.__setattr__(result, "_pred_mat", pred_mat)
        object.__setattr__(result, "_substrate_list", substrate_list)
        return result

    def to_dataframe(self) -> pd.DataFrame:
        """Return a `pred_mat` snapshot isolated from this result."""

        return export_dataframe(self._pred_mat)

    def substrate_list_dataframe(self) -> pd.DataFrame | None:
        """Return an optional substrate-list snapshot isolated from this result."""

        return _export_public_substrate_list(self._substrate_list)


def _require_frame_type(value: object, *, field_name: str) -> None:
    if not isinstance(value, pd.DataFrame):
        raise PhosPyValidationError(f"{field_name} must be a pandas DataFrame")


def _require_optional_frame_type(value: object | None, *, field_name: str) -> None:
    if value is not None and not isinstance(value, pd.DataFrame):
        raise PhosPyValidationError(
            f"{field_name} must be a pandas DataFrame when provided"
        )


def _export_public_substrate_list(table: pd.DataFrame | None) -> pd.DataFrame | None:
    exported = export_optional_dataframe(table)
    if exported is None:
        return None
    legacy_columns = ["kinase", "substrate_site", "score", "rank"]
    if not all(column in exported.columns for column in legacy_columns):
        return exported
    if {"site_key", "display_id"}.issubset(set(exported.columns)):
        return exported
    return exported.loc[:, legacy_columns]


def _substrate_list_uses_encoded_site_keys(table: pd.DataFrame) -> bool:
    if "site_key" not in table.columns:
        return False
    values = table.loc[:, "site_key"].astype(str).tolist()
    if not values:
        return False
    return all(value.startswith("phospy:v1|") for value in values)


__all__ = [
    "KinaseLibraryMotifScoringResult",
    "KinasePredictionResult",
    "KinaseScoringResult",
]
