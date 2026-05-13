from __future__ import annotations

import pytest

from phospy.errors.workflows import SignalomeScaleError
from phospy.science.signalomes.clustering.policies import (
    MAX_FULL_CORRELATION_SITE_COUNT,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
)
from phospy.science.signalomes.clustering.scale_guards import (
    raise_if_exact_tree_limit_exceeded,
    raise_if_full_candidate_scoring_limit_exceeded,
    resolve_max_exact_tree_sites,
)


def test_resolve_max_exact_tree_sites_uses_safe_default() -> None:
    assert resolve_max_exact_tree_sites(None) == MAX_FULL_CORRELATION_SITE_COUNT


def test_raise_if_exact_tree_limit_exceeded_raises_with_context() -> None:
    with pytest.raises(SignalomeScaleError, match="max_exact_tree_sites=10"):
        raise_if_exact_tree_limit_exceeded(
            n_sites=11,
            max_exact_tree_sites=10,
            candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        )


def test_full_candidate_scoring_limit_only_applies_to_full_policy() -> None:
    with pytest.raises(SignalomeScaleError, match="max_full_candidate_scoring_sites=5"):
        raise_if_full_candidate_scoring_limit_exceeded(
            n_sites=6,
            max_full_candidate_scoring_sites=5,
            max_exact_tree_sites=10,
            candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        )

    raise_if_full_candidate_scoring_limit_exceeded(
        n_sites=6,
        max_full_candidate_scoring_sites=5,
        max_exact_tree_sites=10,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    )
