from __future__ import annotations

import numpy as np

from phospy.science.signalomes.clustering.candidate_scoring import (
    _ProfileDegeneracySummary,
)
from phospy.science.signalomes.clustering.candidate_selection import (
    filter_cluster_candidates,
    select_best_candidate_count,
    select_threshold_candidate,
)
from phospy.science.signalomes.clustering.policies import (
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_TREE_ENGINE_EXACT,
)
from phospy.science.signalomes.models import SignalomeClusterCandidateScore


def test_select_best_candidate_count_prefers_lower_count_on_tie() -> None:
    assert select_best_candidate_count({2: 0.5, 3: 0.5, 4: 0.4}) == 2


def test_filter_cluster_candidates_applies_min_median_threshold() -> None:
    scores = {
        2: SignalomeClusterCandidateScore(
            min_median_correlation=0.60,
            mean_median_correlation=0.70,
        ),
        3: SignalomeClusterCandidateScore(
            min_median_correlation=0.40,
            mean_median_correlation=0.90,
        ),
    }

    filtered = filter_cluster_candidates(scores, threshold=0.5)

    assert filtered == {2: 0.70}


def test_select_threshold_candidate_returns_selection_payload() -> None:
    result = select_threshold_candidate(
        candidate_scores={
            2: SignalomeClusterCandidateScore(
                min_median_correlation=0.6,
                mean_median_correlation=0.7,
            ),
            3: SignalomeClusterCandidateScore(
                min_median_correlation=0.5,
                mean_median_correlation=0.7,
            ),
        },
        candidate_labels={
            2: np.asarray([0, 0, 1, 1], dtype=int),
            3: np.asarray([0, 1, 1, 2], dtype=int),
        },
        max_clusters=4,
        threshold=0.5,
        requested_module_count=None,
        reason="selected",
        profile_degeneracy=_ProfileDegeneracySummary(
            zero_variance_count=0,
            near_constant_count=0,
            excluded_count=0,
            excluded_mask=np.zeros(4, dtype=bool),
        ),
        correlation_exclusion_note="",
        approximation_note="",
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_mode=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        exact_cluster_tree_built=True,
        candidate_scoring_evaluated=True,
        candidate_scoring_skip_reason=None,
        candidate_scoring_sampling=None,
    )

    assert result is not None
    assert result.diagnostics.selected_module_count == 2
    assert result.diagnostics.threshold_used == 0.5
