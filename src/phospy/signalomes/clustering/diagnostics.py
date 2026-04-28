"""Helpers for clustering backend diagnostic metadata."""

from __future__ import annotations


def approximation_used_from_candidate_mode(
    *,
    candidate_scoring_mode: str,
    candidate_scoring_evaluated: bool,
) -> bool:
    """Return whether approximate candidate scoring was used."""

    return bool(candidate_scoring_evaluated and candidate_scoring_mode == "sampled")


__all__ = ["approximation_used_from_candidate_mode"]
