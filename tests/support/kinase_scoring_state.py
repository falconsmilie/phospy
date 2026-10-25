"""Populated scoring-state fixtures for reloadable bundle tests."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.kinase.scoring.config import (
    KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
    KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
    KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED,
)
from phospy.kinase.scoring.results import KinaseScoringResult
from phospy.kinase.scoring.selection import (
    KINASE_SCORE_SOURCE_PROFILE_ONLY_NO_MOTIF_OVERLAP,
)
from phospy.kinase.scoring.tables import (
    KINASE_PROFILE_SCORE_DIAGNOSTIC_COLUMNS,
    KINASE_PROFILE_SCORE_DIAGNOSTIC_REASON_INSUFFICIENT_SUBSTRATES_AFTER_LEAVE_ONE_OUT,
    KINASE_PROFILE_SCORE_DIAGNOSTIC_STATUS_UNSCORED,
)

NONDEFAULT_SCORING_MODES = (
    KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
    KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
    KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF,
)


def complete_scoring_result_for_mode(
    base: KinaseScoringResult,
    *,
    scoring_mode: str,
) -> KinaseScoringResult:
    """Return a coherent, fully populated persistence fixture for one mode."""

    profile_scores = base.profile_scores
    motif_scores = profile_scores.mul(0.75)
    motif_scores.columns.name = profile_scores.columns.name
    kinase_library_scores = profile_scores.clip(lower=0.0, upper=1.0)
    kinase_library_scores.columns.name = profile_scores.columns.name
    combined_scores = profile_scores.mul(0.4).add(kinase_library_scores.mul(0.6))
    combined_scores.columns.name = profile_scores.columns.name

    rank_weighted = None
    kinase_library = None
    combined = None
    score_source = "profile_scores"
    score_scale = "fixture_profile_scale"
    if scoring_mode == KINASE_SCORING_MODE_PHOSR_RANK_WEIGHTED:
        rank_weighted = profile_scores.mul(0.9)
        rank_weighted.columns.name = profile_scores.columns.name
        score_source = "rank_weighted_fusion_scores"
        score_scale = "fixture_rank_weighted_scale"
    elif scoring_mode in {
        KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF,
        KINASE_SCORING_MODE_KINASE_LIBRARY_MOTIF_ONLY,
    }:
        kinase_library = kinase_library_scores
        score_source = "kinase_library_motif_scores"
        score_scale = "fixture_kinase_library_scale"
    elif scoring_mode == KINASE_SCORING_MODE_COMBINED_PROFILE_MOTIF:
        kinase_library = kinase_library_scores
        combined = combined_scores
        score_source = "combined_profile_motif_scores"
        score_scale = "fixture_combined_scale"
    else:
        raise AssertionError(f"unsupported scoring fixture mode: {scoring_mode}")

    expected = rank_weighted if rank_weighted is not None else profile_scores
    source_matrix = pd.DataFrame(
        KINASE_SCORE_SOURCE_PROFILE_ONLY_NO_MOTIF_OVERLAP,
        index=expected.index.copy(),
        columns=expected.columns.copy(),
        dtype=object,
    )
    summary = pd.DataFrame(
        {
            "fused_motif_profile_evidence_count": np.zeros(
                expected.shape[1], dtype=np.int64
            ),
            "profile_only_motif_missing_or_constant_count": np.zeros(
                expected.shape[1], dtype=np.int64
            ),
            "profile_only_no_motif_overlap_count": np.full(
                expected.shape[1], expected.shape[0], dtype=np.int64
            ),
            "unavailable_no_score_count": np.zeros(expected.shape[1], dtype=np.int64),
            "sites_with_score_count": np.full(
                expected.shape[1], expected.shape[0], dtype=np.int64
            ),
            "total_sites_count": np.full(
                expected.shape[1], expected.shape[0], dtype=np.int64
            ),
        },
        index=expected.columns.copy(),
    )
    summary.index.name = "kinase"
    profile_diagnostics = pd.DataFrame(
        [
            {
                "kinase": str(profile_scores.columns[0]),
                "substrate_site": str(profile_scores.index[0]),
                "status": KINASE_PROFILE_SCORE_DIAGNOSTIC_STATUS_UNSCORED,
                "reason": KINASE_PROFILE_SCORE_DIAGNOSTIC_REASON_INSUFFICIENT_SUBSTRATES_AFTER_LEAVE_ONE_OUT,
                "substrates_before_leave_one_out": 2,
                "substrates_after_leave_one_out": 1,
                "min_substrates": 2,
            }
        ],
        columns=pd.Index(KINASE_PROFILE_SCORE_DIAGNOSTIC_COLUMNS),
    )
    profile_diagnostics.index.name = "diagnostic_row"
    site_diagnostics = pd.DataFrame(
        {
            "status": ["valid_scored_site"],
            "sequence_length": np.array([15], dtype=np.int64),
        },
        index=pd.Index([profile_scores.index[0]], name=profile_scores.index.name),
    )
    kinase_diagnostics = pd.DataFrame(
        {
            "status": pd.Series(dtype=object),
            "matrix_count": pd.Series(dtype=np.int64),
        },
        index=pd.Index([], name="kinase"),
    )
    fusion_weights = pd.DataFrame(
        np.full(profile_scores.shape, 0.5, dtype=np.float64),
        index=profile_scores.index.copy(),
        columns=profile_scores.columns.copy(),
    )

    return KinaseScoringResult(
        profile_scores=profile_scores,
        motif_scores=motif_scores,
        rank_weighted_fusion_scores=rank_weighted,
        kinase_library_motif_scores=kinase_library,
        combined_profile_motif_scores=combined,
        score_fusion_weights=fusion_weights,
        score_source_matrix=source_matrix,
        score_source_summary=summary,
        profile_score_diagnostics=profile_diagnostics,
        kinase_library_site_diagnostics=site_diagnostics,
        kinase_library_kinase_diagnostics=kinase_diagnostics,
        scoring_mode=scoring_mode,
        score_source=score_source,
        score_scale=score_scale,
        score_scale_metadata=(
            {}
            if scoring_mode == KINASE_SCORING_MODE_KINASE_LIBRARY_CONTEXTUAL_MOTIF
            else {
                "fixture": "complete_scoring_state",
                "mode": scoring_mode,
                "nested": {"weights": [0.4, 0.6]},
            }
        ),
        profile_self_inclusion_policy="leave_one_out",
        motif_sequence_validation=base.motif_sequence_validation,
        motif_library_validation=base.motif_library_validation,
    )


__all__ = [
    "NONDEFAULT_SCORING_MODES",
    "complete_scoring_result_for_mode",
]
