"""Deterministic numerical ranking kernels for SPS discovery."""

from __future__ import annotations

import math
from collections.abc import Sequence

from scipy.stats import chi2


def stable_ranks(
    scored_sites: Sequence[tuple[str, float]],
) -> dict[str, tuple[int, float]]:
    """Return deterministic competition ranks and PhosR-compatible midranks."""

    ordered = sorted(scored_sites, key=lambda item: (item[1], item[0]))
    ranks: dict[str, tuple[int, float]] = {}
    start = 0
    while start < len(ordered):
        stop = start + 1
        while stop < len(ordered) and ordered[stop][1] == ordered[start][1]:
            stop += 1
        competition_rank = start + 1
        average_rank = (float(start + 1) + float(stop)) / 2.0
        for site_key, _ in ordered[start:stop]:
            ranks[site_key] = (competition_rank, average_rank)
        start = stop
    return ranks


def fisher_style_consensus_score(quantiles: Sequence[float]) -> float:
    contribution_count = len(quantiles)
    statistic = -2.0 * math.fsum(math.log(value) for value in quantiles)
    degrees_of_freedom = 2 * (contribution_count - 1)
    return float(chi2.sf(statistic, degrees_of_freedom))


__all__: list[str] = []
