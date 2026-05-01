"""Scale guards for exact tree construction and candidate scoring."""

from __future__ import annotations

from phospy.errors.workflows import SignalomeScaleError
from phospy.signalomes.clustering.policies import (
    MAX_FULL_CORRELATION_SITE_COUNT,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
)


def resolve_max_exact_tree_sites(max_exact_tree_sites: int | None) -> int:
    """Resolve exact-tree guard limit; `None` maps to the safe default limit."""

    resolved = (
        MAX_FULL_CORRELATION_SITE_COUNT
        if max_exact_tree_sites is None
        else int(max_exact_tree_sites)
    )
    if resolved < 1:
        raise ValueError("max_exact_tree_sites must be >= 1")
    return resolved


def raise_if_exact_tree_limit_exceeded(
    *,
    n_sites: int,
    max_exact_tree_sites: int | None,
    candidate_scoring_policy: str,
) -> None:
    """Raise when exact cluster-tree construction exceeds configured scale guard."""

    resolved_max_exact_tree_sites = resolve_max_exact_tree_sites(max_exact_tree_sites)
    if n_sites > int(resolved_max_exact_tree_sites):
        raise SignalomeScaleError(
            "Signalome exact cluster-tree construction received "
            f"{n_sites:,} sites, which exceeds max_exact_tree_sites="
            f"{int(resolved_max_exact_tree_sites):,} "
            "(tree_implementation='exact_cluster_tree'). "
            f"candidate_scoring_policy='{candidate_scoring_policy}' still "
            "requires exact cluster-tree construction in the current "
            "implementation."
        )


def raise_if_full_candidate_scoring_limit_exceeded(
    *,
    n_sites: int,
    max_full_candidate_scoring_sites: int,
    max_exact_tree_sites: int | None,
    candidate_scoring_policy: str,
) -> None:
    """Raise when full candidate-correlation scoring exceeds configured limit."""

    resolved_max_exact_tree_sites = resolve_max_exact_tree_sites(max_exact_tree_sites)
    if (
        candidate_scoring_policy == SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
        and n_sites > int(max_full_candidate_scoring_sites)
        and n_sites <= int(resolved_max_exact_tree_sites)
    ):
        raise SignalomeScaleError(
            "Signalome full candidate-correlation scoring would evaluate "
            f"{n_sites:,} sites, which exceeds configured "
            f"max_full_candidate_scoring_sites={int(max_full_candidate_scoring_sites):,}. "
            "Exact cluster-tree construction has not been attempted for this "
            "request. Use candidate_scoring_policy='sampled' for candidate "
            "module-count evaluation, reduce interpreted sites, or increase "
            "max_full_candidate_scoring_sites deliberately."
        )


__all__ = [
    "raise_if_exact_tree_limit_exceeded",
    "raise_if_full_candidate_scoring_limit_exceeded",
    "resolve_max_exact_tree_sites",
]
