"""Candidate-shortfall diagnostics for Kinase prediction."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.kinase.policies import ThresholdMode
from phospy.kinase.prediction.candidates import _threshold_pass_mask


@dataclass(frozen=True, slots=True)
class CandidateShortfallDiagnostics:
    """Summary diagnostics when strict candidate rules leave no kinases."""

    kinase_count: int
    site_count: int
    effective_top: int
    qualifying_kinases: int
    max_qualifying_sites: int
    near_miss_kinases: tuple[tuple[str, int], ...]


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


__all__ = ["CandidateShortfallDiagnostics", "summarize_candidate_shortfall"]
