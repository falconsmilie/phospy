"""Dataset scientific table wrappers."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from phospy.errors.validation import DatasetValidationError
from phospy.tables.base import TableSchema
from phospy.validation.common.dataframes import (
    require_canonical_site_index,
    require_columns,
    require_dataframe,
    require_exact_index_match,
    require_finite_numeric_dataframe,
    require_non_empty_string_column,
    require_numeric_dataframe,
    require_site_identity_coherence,
    require_unique_columns,
    require_unique_index,
)


@dataclass(frozen=True, slots=True)
class PhosphoIntensityMatrix(TableSchema):
    """Schema wrapper for ``dataset.phospho``."""

    allow_missing: bool = field(default=False, repr=False, compare=False)

    _field_name = "dataset.phospho"
    _error_type = DatasetValidationError

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        require_dataframe(
            frame,
            field_name=self._field_name,
            allow_empty=False,
            error_type=self._error_type,
        )
        require_unique_columns(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        require_numeric_dataframe(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        require_finite_numeric_dataframe(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
            allow_missing=self.allow_missing,
        )
        require_unique_index(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        require_canonical_site_index(
            frame.index,
            field_name=f"{self._field_name}.index",
            error_type=self._error_type,
            strict_supported_format=True,
        )
        return frame


@dataclass(frozen=True, slots=True)
class SiteMetadataTable(TableSchema):
    """Schema wrapper for ``dataset.site_metadata``."""

    expected_index: pd.Index | None = field(default=None, repr=False, compare=False)
    require_site_sequence: bool = field(default=False, repr=False, compare=False)

    _field_name = "dataset.site_metadata"
    _error_type = DatasetValidationError

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        require_dataframe(
            frame,
            field_name=self._field_name,
            allow_empty=False,
            error_type=self._error_type,
        )
        require_canonical_site_index(
            frame.index,
            field_name=f"{self._field_name}.index",
            error_type=self._error_type,
            strict_supported_format=True,
        )
        if self.expected_index is not None:
            require_exact_index_match(
                left=frame.index,
                right=self.expected_index,
                left_name=f"{self._field_name}.index",
                right_name="dataset.phospho.index",
                error_type=self._error_type,
            )
        require_columns(
            frame,
            field_name=self._field_name,
            required_columns=("gene_symbol", "site"),
            error_type=self._error_type,
        )
        if self.require_site_sequence:
            require_columns(
                frame,
                field_name=self._field_name,
                required_columns=("site_sequence",),
                error_type=self._error_type,
            )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name="gene_symbol",
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name="site",
            error_type=self._error_type,
        )
        if "site_sequence" in frame.columns:
            require_non_empty_string_column(
                frame,
                field_name=self._field_name,
                column_name="site_sequence",
                error_type=self._error_type,
            )
        if "protein_id" in frame.columns:
            require_non_empty_string_column(
                frame,
                field_name=self._field_name,
                column_name="protein_id",
                error_type=self._error_type,
            )
        require_site_identity_coherence(
            site_index=frame.index,
            site_metadata=frame,
            site_index_field_name=f"{self._field_name}.index",
            site_metadata_field_name=self._field_name,
            error_type=self._error_type,
        )
        return frame


@dataclass(frozen=True, slots=True)
class SampleMetadataTable(TableSchema):
    """Schema wrapper for ``dataset.sample_metadata``."""

    expected_index: pd.Index | None = field(default=None, repr=False, compare=False)

    _field_name = "dataset.sample_metadata"
    _error_type = DatasetValidationError

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        require_dataframe(
            frame,
            field_name=self._field_name,
            allow_empty=False,
            error_type=self._error_type,
        )
        require_unique_index(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        if self.expected_index is not None:
            require_exact_index_match(
                left=frame.index,
                right=self.expected_index,
                left_name=f"{self._field_name}.index",
                right_name="dataset.phospho.columns",
                error_type=self._error_type,
            )
        return frame


@dataclass(frozen=True, slots=True)
class TotalProteinMatrix(TableSchema):
    """Schema wrapper for ``dataset.total``."""

    expected_sample_index: pd.Index | None = field(
        default=None,
        repr=False,
        compare=False,
    )

    _field_name = "dataset.total"
    _error_type = DatasetValidationError

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        require_dataframe(
            frame,
            field_name=self._field_name,
            allow_empty=False,
            error_type=self._error_type,
        )
        require_numeric_dataframe(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        require_finite_numeric_dataframe(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
            allow_missing=False,
        )
        require_unique_index(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        require_unique_columns(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        if self.expected_sample_index is not None:
            require_exact_index_match(
                left=frame.columns,
                right=self.expected_sample_index,
                left_name=f"{self._field_name}.columns",
                right_name="dataset.phospho.columns",
                error_type=self._error_type,
            )
        return frame
