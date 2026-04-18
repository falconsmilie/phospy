"""Validator for the analysis-ready dataset boundary."""

from __future__ import annotations

import pandas as pd

from phospy.errors.validation import DatasetValidationError
from phospy.references.models import Organism
from phospy.transformations.models import TransformationState
from phospy.validation.common.dataframes import (
    require_columns,
    require_dataframe,
    require_exact_index_match,
    require_non_empty_string_column,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
)
from phospy.validation.common.missing_values import (
    MissingValuePolicy,
    require_missing_value_policy,
)
from phospy.validation.transformations.state import TransformationStateValidator


class AnalysisReadyDatasetValidator:
    """Validate the public `AnalysisReadyPhosphoDataset` contract."""

    _REQUIRED_SITE_COLUMNS = ("gene_symbol", "site", "site_sequence")

    def __init__(
        self,
        transformation_validator: TransformationStateValidator | None = None,
    ) -> None:
        self._transformation_validator = (
            transformation_validator or TransformationStateValidator()
        )

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        total: pd.DataFrame | None,
        organism: Organism | None,
        transformation_state: TransformationState,
    ) -> None:
        phospho_frame = require_dataframe(
            phospho,
            field_name="dataset.phospho",
            allow_empty=False,
            error_type=DatasetValidationError,
        )
        require_numeric_dataframe(
            phospho_frame,
            field_name="dataset.phospho",
            error_type=DatasetValidationError,
        )
        require_missing_value_policy(
            phospho_frame,
            field_name="dataset.phospho",
            policy=MissingValuePolicy.FORBID,
            error_type=DatasetValidationError,
        )
        require_unique_index(
            phospho_frame,
            field_name="dataset.phospho",
            error_type=DatasetValidationError,
        )
        require_unique_columns(
            phospho_frame,
            field_name="dataset.phospho",
            error_type=DatasetValidationError,
        )

        site_metadata_frame = require_dataframe(
            site_metadata,
            field_name="dataset.site_metadata",
            allow_empty=False,
            error_type=DatasetValidationError,
        )
        require_exact_index_match(
            left=site_metadata_frame.index,
            right=phospho_frame.index,
            left_name="dataset.site_metadata.index",
            right_name="dataset.phospho.index",
            error_type=DatasetValidationError,
        )
        require_columns(
            site_metadata_frame,
            field_name="dataset.site_metadata",
            required_columns=self._REQUIRED_SITE_COLUMNS,
            error_type=DatasetValidationError,
        )
        require_non_empty_string_column(
            site_metadata_frame,
            field_name="dataset.site_metadata",
            column_name="site_sequence",
            error_type=DatasetValidationError,
        )

        if sample_metadata is not None:
            sample_metadata_frame = require_dataframe(
                sample_metadata,
                field_name="dataset.sample_metadata",
                allow_empty=False,
                error_type=DatasetValidationError,
            )
            require_exact_index_match(
                left=sample_metadata_frame.index,
                right=phospho_frame.columns,
                left_name="dataset.sample_metadata.index",
                right_name="dataset.phospho.columns",
                error_type=DatasetValidationError,
            )

        if total is not None:
            total_frame = require_dataframe(
                total,
                field_name="dataset.total",
                allow_empty=False,
                error_type=DatasetValidationError,
            )
            require_numeric_dataframe(
                total_frame,
                field_name="dataset.total",
                error_type=DatasetValidationError,
            )
            require_missing_value_policy(
                total_frame,
                field_name="dataset.total",
                policy=MissingValuePolicy.FORBID,
                error_type=DatasetValidationError,
            )
            require_unique_index(
                total_frame,
                field_name="dataset.total",
                error_type=DatasetValidationError,
            )
            require_exact_index_match(
                left=total_frame.columns,
                right=phospho_frame.columns,
                left_name="dataset.total.columns",
                right_name="dataset.phospho.columns",
                error_type=DatasetValidationError,
            )

        if organism is not None and not isinstance(organism, Organism):
            raise DatasetValidationError(
                "dataset.organism must be an Organism enum value or None"
            )

        self._transformation_validator.run(
            transformation_state=transformation_state,
            has_total_matrix=total is not None,
        )
