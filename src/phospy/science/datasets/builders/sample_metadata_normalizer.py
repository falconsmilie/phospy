"""Sample-metadata focused collaborator for dataset convention normalisation."""

from __future__ import annotations

import pandas as pd

from phospy.science.datasets.builders.normalization_reporter import (
    SAMPLE_LABEL_INDEX_POLICY,
    DatasetConventionNormalisationReporter,
)


class SampleMetadataNormalizer:
    """Normalize sample-metadata index labels and align with phospho columns."""

    def __init__(
        self,
        *,
        reporter: DatasetConventionNormalisationReporter | None = None,
    ) -> None:
        self._reporter = reporter or DatasetConventionNormalisationReporter()

    def run(
        self,
        sample_metadata: pd.DataFrame | None,
        *,
        phospho_columns: pd.Index,
    ) -> pd.DataFrame | None:
        if sample_metadata is None:
            return None
        normalized = sample_metadata
        normalized.index = self._reporter.normalize_index_labels(
            normalized.index,
            field_name="dataset build request sample_metadata.index",
            policy=SAMPLE_LABEL_INDEX_POLICY,
        )
        if (
            not normalized.index.equals(phospho_columns)
            and normalized.index.isin(phospho_columns).all()
            and phospho_columns.isin(normalized.index).all()
        ):
            return normalized.reindex(phospho_columns)
        return normalized
