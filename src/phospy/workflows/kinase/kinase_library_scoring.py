"""Kinase Library motif scoring adapter for kinase workflow execution."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass

import pandas as pd

from phospy.science.prediction.motif_scoring import (
    SEQUENCE_SEMANTICS_CENTRED_SEQUENCE,
    KinaseLibraryMotifScoringResult,
    score_kinase_library_motifs,
)
from phospy.science.prediction.motif_scoring.scaling import minmax_scale_columns
from phospy.science.references.kinase_library_models import KinaseLibraryResource

KINASE_LIBRARY_WORKFLOW_SCORE_SCALE = "kinase_library_motif_minmax_unit_interval"


@dataclass(frozen=True, slots=True)
class KinaseLibraryWorkflowScoringResult:
    """Workflow-normalized Kinase Library motif scoring outputs."""

    scores: pd.DataFrame
    raw_scores: pd.DataFrame
    site_diagnostics: pd.DataFrame
    kinase_diagnostics: pd.DataFrame
    score_scale_metadata: Mapping[str, object]


class KinaseLibraryWorkflowScorer:
    """Run pure motif scoring and adapt it to workflow support scores.

    The pure scorer returns provider-scale raw values. This adapter preserves
    those raw scores in diagnostics and exposes a per-kinase min-max normalized
    matrix for workflow prediction support.
    """

    def __init__(
        self,
        *,
        score_motifs: Callable[..., KinaseLibraryMotifScoringResult] = (
            score_kinase_library_motifs
        ),
    ) -> None:
        self._score_motifs = score_motifs

    def run(
        self,
        *,
        resource: KinaseLibraryResource,
        site_sequences: pd.Series,
        site_identities: pd.Series,
        site_index: Sequence[object],
    ) -> KinaseLibraryWorkflowScoringResult:
        library_result = self._score_motifs(
            site_sequences=site_sequences,
            site_identities=site_identities,
            matrices=resource,
            site_index=tuple(str(site_id) for site_id in site_index),
            sequence_semantics=SEQUENCE_SEMANTICS_CENTRED_SEQUENCE,
            higher_is_better=True,
        )
        scores = minmax_scale_columns(library_result.raw_scores)
        metadata = {
            **library_result.score_scale_payload(),
            "resource_source_name": resource.source_name,
            "resource_source_version": resource.source_version,
            "resource_score_scale": resource.score_scale,
            "workflow_score_scale": KINASE_LIBRARY_WORKFLOW_SCORE_SCALE,
            "workflow_score_transform": "per-kinase min-max scaling of raw scores",
            "authoritative_matrix": "kinase_library_motif_scores",
        }
        return KinaseLibraryWorkflowScoringResult(
            scores=scores,
            raw_scores=library_result.raw_scores,
            site_diagnostics=library_result.site_diagnostics,
            kinase_diagnostics=library_result.kinase_diagnostics,
            score_scale_metadata=metadata,
        )


__all__ = [
    "KINASE_LIBRARY_WORKFLOW_SCORE_SCALE",
    "KinaseLibraryWorkflowScorer",
    "KinaseLibraryWorkflowScoringResult",
]
