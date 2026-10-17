"""Passive Kinase result-boundary contracts used by downstream capabilities."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

import pandas as pd

from phospy.errors.validation import PhosPyValidationError
from phospy.frames.ownership import own_dataframe
from phospy.frames.table_schema import TableSchema
from phospy.kinase._matrix_validation import (
    _validate_kinase_score_like_matrix,  # pyright: ignore[reportPrivateUsage] - KinaseScoreMatrix and KinasePredictionMatrix intentionally share this singular Kinase-owned matrix validator.
)
from phospy.policies import PolicyEnum
from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)


class DownstreamScoreSource(PolicyEnum):
    PROFILE_SCORES = "profile_scores"
    RANK_WEIGHTED_FUSION_SCORES = "rank_weighted_fusion_scores"
    KINASE_LIBRARY_MOTIF_SCORES = "kinase_library_motif_scores"
    COMBINED_PROFILE_MOTIF_SCORES = "combined_profile_motif_scores"


@dataclass(frozen=True, slots=True)
class DownstreamScoreSelectionPolicy:
    """Versioned policy describing downstream score-lane selection semantics."""

    name: str
    version: str
    parameters: Mapping[str, object]
    description: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "parameters",
            MappingProxyType(
                {str(key): value for key, value in self.parameters.items()}
            ),
        )

    @property
    def record(self) -> ScientificPolicyRecord:
        return ScientificPolicyRecord(
            id=ScientificPolicyId.SIGNALOME_DOWNSTREAM_SCORE_SELECTION,
            name=self.name,
            version=self.version,
            description=self.description,
            parameters=self.parameters,
            assumptions=(
                "Downstream score-lane selection changes which upstream signal is "
                "interpreted by signalome stages.",
            ),
            output_scale="Selected downstream kinase-site score matrix.",
            quantitative_meaning="relative_downstream_support",
        )


SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY = (
    DownstreamScoreSelectionPolicy(
        name="signalome_downstream_score_rank_weighted_preferred_v1",
        version="1",
        parameters={
            "preferred_source": (
                DownstreamScoreSource.RANK_WEIGHTED_FUSION_SCORES.value
            ),
            "fallback_source": DownstreamScoreSource.PROFILE_SCORES.value,
        },
        description=(
            "Select PhosR-inspired rank-weighted fusion scores when available; "
            "otherwise use profile-only scores."
        ),
    )
)


@dataclass(frozen=True, slots=True, eq=False)
class KinaseScoreMatrix(TableSchema):
    """Schema wrapper for kinase scoring matrices."""

    field_name: str = field(
        default="scoring_result.profile_scores",
        repr=False,
        compare=False,
    )
    allow_missing: bool = field(default=True, repr=False, compare=False)
    enforce_unit_interval: bool = field(default=False, repr=False, compare=False)
    require_site_index: bool = field(default=True, repr=False, compare=False)

    _field_name = "scoring_result.profile_scores"
    _error_type = PhosPyValidationError

    def _own_frame(self, assume_owned: bool) -> pd.DataFrame:
        return own_dataframe(
            self.frame,
            field_name=self.field_name,
            error_type=self._error_type,
            assume_owned=assume_owned,
        )

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        _validate_kinase_score_like_matrix(
            frame=frame,
            field_name=self.field_name,
            error_type=self._error_type,
            allow_missing=self.allow_missing,
            enforce_unit_interval=self.enforce_unit_interval,
            require_site_index=self.require_site_index,
        )
        return frame


KINASE_NON_DEFAULT_REFERENCE_SOURCE_CAVEAT_CODE = "kinase_non_default_reference_source"
KINASE_REFERENCE_SCORE_FALLBACK_CAVEAT_CODE = "kinase_reference_score_fallback"
KINASE_SCORING_LIMITATION_CAVEAT_CODE = "kinase_non_phosr_equivalent_scoring"
KINASE_DOWNSTREAM_LIMITATION_CAVEAT_CODES = frozenset(
    {
        KINASE_NON_DEFAULT_REFERENCE_SOURCE_CAVEAT_CODE,
        KINASE_REFERENCE_SCORE_FALLBACK_CAVEAT_CODE,
        KINASE_SCORING_LIMITATION_CAVEAT_CODE,
    }
)


__all__ = [
    "DownstreamScoreSelectionPolicy",
    "DownstreamScoreSource",
    "KINASE_DOWNSTREAM_LIMITATION_CAVEAT_CODES",
    "KINASE_NON_DEFAULT_REFERENCE_SOURCE_CAVEAT_CODE",
    "KINASE_REFERENCE_SCORE_FALLBACK_CAVEAT_CODE",
    "KINASE_SCORING_LIMITATION_CAVEAT_CODE",
    "KinaseScoreMatrix",
    "SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY",
]
