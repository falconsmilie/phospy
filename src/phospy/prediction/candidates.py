from __future__ import annotations

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
    "CandidateSelector",
    "build_candidate_substrate_list",
]
