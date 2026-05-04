"""Candidate selection helpers for kinase prediction."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.policy_models import ThresholdMode


@dataclass(frozen=True, slots=True)
class CandidateShortfallDiagnostics:
    """Summary diagnostics when strict candidate rules leave no kinases."""

    kinase_count: int
    site_count: int
    effective_top: int
    qualifying_kinases: int
    max_qualifying_sites: int
    near_miss_kinases: tuple[tuple[str, int], ...]


def build_candidate_substrate_list(
    *,
    scores: pd.DataFrame,
    top: int,
    score_threshold: float,
    inclusion: int,
    allowed_sites_by_kinase: Mapping[str, list[str]] | None = None,
    threshold_mode: ThresholdMode | str = ThresholdMode.GREATER_THAN,
) -> dict[str, list[str]]:
    """Select per-kinase candidate phosphosites from a score matrix."""

    if top < 1:
        raise ValueError("top must be >= 1")
    if inclusion < 1:
        raise ValueError("inclusion must be >= 1")
    if not 0.0 <= float(score_threshold) <= 1.0:
        raise ValueError("score_threshold must be within [0.0, 1.0]")
    if scores.empty:
        return {}

    score_values = scores.to_numpy(dtype=float, copy=False)
    site_ids = scores.index.to_numpy(copy=False)
    kinase_labels = scores.columns.to_numpy(copy=False)
    top_k = min(int(top), score_values.shape[0])
    substrate_list: dict[str, list[str]] = {}
    resolved_threshold_mode = ThresholdMode.parse(
        threshold_mode,
        field_name="candidate_substrate_selection.threshold_mode",
    )

    for kinase_position, kinase in enumerate(kinase_labels):
        kinase_key = str(kinase)
        kinase_scores = score_values[:, kinase_position]
        finite_mask = np.isfinite(kinase_scores)
        if allowed_sites_by_kinase is not None:
            allowed_sites = set(allowed_sites_by_kinase.get(kinase_key, []))
            allowed_mask = np.fromiter(
                (site in allowed_sites for site in site_ids),
                dtype=bool,
                count=len(site_ids),
            )
            finite_mask &= allowed_mask
        valid_positions = np.flatnonzero(finite_mask)
        if valid_positions.size == 0:
            continue
        valid_top_k = min(top_k, int(valid_positions.size))
        if valid_top_k == 0:
            continue
        local_scores = kinase_scores[valid_positions]
        local_top_positions = np.argpartition(-local_scores, valid_top_k - 1)[
            :valid_top_k
        ]
        selected_positions = valid_positions[local_top_positions]
        selected_scores = kinase_scores[selected_positions]
        ranked_positions = selected_positions[
            np.lexsort((selected_positions, -selected_scores))
        ]
        qualifying_positions = ranked_positions[
            _threshold_pass_mask(
                values=kinase_scores[ranked_positions],
                score_threshold=float(score_threshold),
                threshold_mode=resolved_threshold_mode,
            )
        ]
        sites = site_ids[qualifying_positions].tolist()
        if len(sites) >= inclusion:
            substrate_list[kinase_key] = [str(site) for site in sites]
    return substrate_list


def summarize_candidate_shortfall(
    *,
    scores: pd.DataFrame,
    top: int,
    score_threshold: float,
    inclusion: int,
    max_examples: int = 3,
    threshold_mode: ThresholdMode | str = ThresholdMode.GREATER_THAN,
) -> CandidateShortfallDiagnostics:
    """Summarize candidate loss for actionable boundary diagnostics."""

    score_values = scores.to_numpy(dtype=float, copy=False)
    site_count, kinase_count = score_values.shape
    effective_top = min(int(top), site_count) if site_count > 0 else 0
    if site_count == 0 or kinase_count == 0 or effective_top == 0:
        return CandidateShortfallDiagnostics(
            kinase_count=kinase_count,
            site_count=site_count,
            effective_top=effective_top,
            qualifying_kinases=0,
            max_qualifying_sites=0,
            near_miss_kinases=(),
        )

    qualifying_counts = np.zeros(kinase_count, dtype=int)
    kinase_labels = scores.columns.to_numpy(dtype=object, copy=False)
    resolved_threshold_mode = ThresholdMode.parse(
        threshold_mode,
        field_name="candidate_substrate_selection.threshold_mode",
    )
    for kinase_position in range(kinase_count):
        kinase_scores = score_values[:, kinase_position]
        finite_positions = np.flatnonzero(np.isfinite(kinase_scores))
        if finite_positions.size == 0:
            continue
        valid_top = min(effective_top, int(finite_positions.size))
        top_positions = finite_positions[
            np.argpartition(-kinase_scores[finite_positions], valid_top - 1)[:valid_top]
        ]
        top_values = kinase_scores[top_positions]
        qualifying_counts[kinase_position] = int(
            np.count_nonzero(
                _threshold_pass_mask(
                    values=top_values,
                    score_threshold=float(score_threshold),
                    threshold_mode=resolved_threshold_mode,
                )
            )
        )

    near_miss_positions = np.flatnonzero(
        (qualifying_counts > 0) & (qualifying_counts < inclusion)
    )
    near_miss_kinases = tuple(
        sorted(
            (
                (str(kinase_labels[position]), int(qualifying_counts[position]))
                for position in near_miss_positions
            ),
            key=lambda item: (-item[1], item[0]),
        )[:max_examples]
    )
    return CandidateShortfallDiagnostics(
        kinase_count=kinase_count,
        site_count=site_count,
        effective_top=effective_top,
        qualifying_kinases=int(np.count_nonzero(qualifying_counts > 0)),
        max_qualifying_sites=int(qualifying_counts.max()) if kinase_count > 0 else 0,
        near_miss_kinases=near_miss_kinases,
    )


def _threshold_pass_mask(
    *,
    values: np.ndarray,
    score_threshold: float,
    threshold_mode: ThresholdMode,
) -> np.ndarray:
    if threshold_mode is ThresholdMode.GREATER_THAN:
        return values > score_threshold
    if threshold_mode is ThresholdMode.GREATER_THAN_OR_EQUAL:
        return values >= score_threshold
    return values > score_threshold


__all__ = [
    "CandidateShortfallDiagnostics",
    "build_candidate_substrate_list",
    "summarize_candidate_shortfall",
]
