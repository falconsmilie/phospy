"""Validator for the analysis-ready dataset boundary."""

from __future__ import annotations

import pandas as pd

from phospy.errors.validation import DatasetValidationError
from phospy.references.models import Organism
from phospy.tables.datasets import (
    PhosphoIntensityMatrix,
    SampleMetadataTable,
    SiteMetadataTable,
    TotalProteinMatrix,
)
from phospy.validation.common.dataframes import (
    require_dataframe,
    require_exact_index_match,
    require_finite_numeric_dataframe,
    require_non_empty_dataframe,
    require_numeric_dataframe,
    require_unique_columns,
)


class AnalysisReadyDatasetValidator:
    """Validate the public `AnalysisReadyPhosphoDataset` contract."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        total: pd.DataFrame | None,
        comparisons: pd.DataFrame | None,
        organism: Organism | None,
    ) -> None:
        phospho_frame = PhosphoIntensityMatrix(
            frame=phospho,
            _assume_owned=True,
        ).frame
        SiteMetadataTable(
            frame=site_metadata,
            expected_index=phospho_frame.index,
            _assume_owned=True,
        )

        if sample_metadata is not None:
            SampleMetadataTable(
                frame=sample_metadata,
                expected_index=phospho_frame.columns,
                _assume_owned=True,
            )

        if comparisons is not None:
            comparisons_frame = require_dataframe(
                comparisons,
                field_name="dataset.comparisons",
                allow_empty=True,
                error_type=DatasetValidationError,
            )
            require_non_empty_dataframe(
                comparisons_frame,
                field_name="dataset.comparisons",
                error_type=DatasetValidationError,
            )
            require_numeric_dataframe(
                comparisons_frame,
                field_name="dataset.comparisons",
                error_type=DatasetValidationError,
            )
            require_finite_numeric_dataframe(
                comparisons_frame,
                field_name="dataset.comparisons",
                error_type=DatasetValidationError,
                allow_missing=False,
            )
            require_unique_columns(
                comparisons_frame,
                field_name="dataset.comparisons",
                error_type=DatasetValidationError,
            )
            require_exact_index_match(
                left=comparisons_frame.index,
                right=phospho_frame.index,
                left_name="dataset.comparisons.index",
                right_name="dataset.phospho.index",
                error_type=DatasetValidationError,
            )

        if total is not None:
            TotalProteinMatrix(
                frame=total,
                expected_sample_index=phospho_frame.columns,
                _assume_owned=True,
            )

        if organism is not None and not isinstance(organism, Organism):
            raise DatasetValidationError(
                "dataset.organism must be an Organism enum value or None"
            )
