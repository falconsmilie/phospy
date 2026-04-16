from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..errors import TableSchemaError
from ..internal.defaults import (
    DEFAULT_PREDICTION_INCLUSION,
    DEFAULT_PREDICTION_SCORE_THRESHOLD,
    DEFAULT_PREDICTION_TOP,
)
from ..validation.schema.tables import PredictionScoreMatrixSchema
from .validation import validate_positive_int


@dataclass(frozen=True, slots=True)
class CandidateShortfallDiagnostics:
    """Summary diagnostics for strict-threshold candidate shortfalls."""

    kinase_count: int
    site_count: int
    effective_top: int
    qualifying_kinases: int
    max_qualifying_sites: int
    near_miss_kinases: tuple[tuple[str, int], ...]


class CandidateSelector:
    def select(
        self,
        combined_scores: pd.DataFrame,
        *,
        top: int,
        score_threshold: float,
        inclusion: int,
    ) -> dict[str, list[str]]:
        return _build_candidate_substrate_list(
            combined_scores,
            top=top,
            score_threshold=score_threshold,
            inclusion=inclusion,
        )


def _build_candidate_substrate_list(
    combined_scores: pd.DataFrame,
    *,
    top: int,
    score_threshold: float,
    inclusion: int,
) -> dict[str, list[str]]:
    """Build per-kinase candidate lists from combined score matrix.

    Implementation note:
    We intentionally keep per-column ``np.argpartition`` instead of a whole-matrix
    vectorization pass. A reviewed full-matrix attempt regressed runtime across
    representative matrix sizes due to extra global sorting/materialization
    overhead. See ADR 0005 for benchmark evidence and decision rationale.
    """

    substrate_list: dict[str, list[str]] = {}
    score_values = combined_scores.to_numpy(dtype=float, copy=False)
    site_ids = combined_scores.index.to_numpy(copy=False)
    kinase_labels = combined_scores.columns.to_numpy(copy=False)
    top_k = min(int(top), score_values.shape[0])

    for kinase_position, kinase in enumerate(kinase_labels):
        kinase_scores = score_values[:, kinase_position]
        top_positions = np.argpartition(-kinase_scores, top_k - 1)[:top_k]
        top_values = kinase_scores[top_positions]
        ranked_positions = top_positions[np.lexsort((top_positions, -top_values))]
        qualifying_positions = ranked_positions[
            kinase_scores[ranked_positions] > score_threshold
        ]
        sites = site_ids[qualifying_positions].tolist()
        if len(sites) >= inclusion:
            substrate_list[kinase] = sites
    return substrate_list


def summarize_candidate_shortfall(
    combined_scores: pd.DataFrame,
    *,
    top: int,
    score_threshold: float,
    inclusion: int,
    max_examples: int = 3,
) -> CandidateShortfallDiagnostics:
    """Summarize why strict candidate filters produced no qualifying kinases."""

    score_values = combined_scores.to_numpy(dtype=float, copy=False)
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
    kinase_labels = combined_scores.columns.to_numpy(dtype=object, copy=False)
    for kinase_position in range(kinase_count):
        kinase_scores = score_values[:, kinase_position]
        top_positions = np.argpartition(-kinase_scores, effective_top - 1)[
            :effective_top
        ]
        top_values = kinase_scores[top_positions]
        qualifying_counts[kinase_position] = int(
            np.count_nonzero(top_values > score_threshold)
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
            key=lambda value: (-value[1], value[0]),
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


def build_candidate_substrate_list(
    combined_scores: pd.DataFrame,
    top: int = DEFAULT_PREDICTION_TOP,
    score_threshold: float = DEFAULT_PREDICTION_SCORE_THRESHOLD,
    inclusion: int = DEFAULT_PREDICTION_INCLUSION,
) -> dict[str, list[str]]:
    """Select candidate kinase substrates from the combined score matrix."""

    validate_positive_int(top, name="top")
    validate_positive_int(inclusion, name="inclusion")
    if not 0.0 <= float(score_threshold) <= 1.0:
        msg = "score_threshold must be between 0.0 and 1.0"
        raise TableSchemaError(msg)
    validated_scores = PredictionScoreMatrixSchema.validate(
        combined_scores,
        context="combined_scores",
    )
    return _build_candidate_substrate_list(
        validated_scores,
        top=top,
        score_threshold=score_threshold,
        inclusion=inclusion,
    )


__all__ = [
    "CandidateShortfallDiagnostics",
    "CandidateSelector",
    "build_candidate_substrate_list",
    "summarize_candidate_shortfall",
]
