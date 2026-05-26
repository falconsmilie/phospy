"""Total-protein matrix normalizer for dataset convention handling."""

from __future__ import annotations

import pandas as pd

from phospy.science.datasets.builders.normalization_reporter import (
    SAMPLE_LABEL_INDEX_POLICY,
    DatasetConventionNormalisationReporter,
)


class TotalProteinMatrixNormalizer:
    """Normalize total-matrix labels and align sample columns with phospho."""

    def __init__(
        self,
        *,
        reporter: DatasetConventionNormalisationReporter | None = None,
    ) -> None:
        self._reporter = reporter or DatasetConventionNormalisationReporter()

    def run(
        self,
        total: pd.DataFrame | None,
        *,
        phospho_columns: pd.Index,
    ) -> pd.DataFrame | None:
        if total is None:
            return None
        normalized = total
        normalized.index = self._reporter.normalize_index_labels(
            normalized.index,
            field_name="dataset build request total.index",
            policy=SAMPLE_LABEL_INDEX_POLICY,
        )
        normalized.columns = self._reporter.normalize_index_labels(
            normalized.columns,
            field_name="dataset build request total.columns",
            policy=SAMPLE_LABEL_INDEX_POLICY,
        )
        if (
            not normalized.columns.equals(phospho_columns)
            and normalized.columns.isin(phospho_columns).all()
            and phospho_columns.isin(normalized.columns).all()
        ):
            normalized = normalized.reindex(columns=phospho_columns)
        return normalized
