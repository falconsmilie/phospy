"""Reference scientific table wrappers."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.validation import ReferenceValidationError
from phospy.tables.base import TableSchema
from phospy.validation.common.dataframes import (
    require_canonical_site_index,
    require_canonical_site_series,
    require_canonical_string_column,
    require_columns,
    require_dataframe,
    require_non_empty_string_column,
    require_unique_index,
    require_unique_row_pairs,
)


@dataclass(frozen=True, slots=True)
class KinaseSubstrateReference(TableSchema):
    """Schema wrapper for ``references.kinase_substrate_map``."""

    _field_name = "references.kinase_substrate_map"
    _error_type = ReferenceValidationError

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        require_dataframe(
            frame,
            field_name=self._field_name,
            allow_empty=False,
            error_type=self._error_type,
        )
        require_columns(
            frame,
            field_name=self._field_name,
            required_columns=("kinase", "substrate_site"),
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name="kinase",
            error_type=self._error_type,
        )
        require_canonical_string_column(
            frame,
            field_name=self._field_name,
            column_name="kinase",
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name="substrate_site",
            error_type=self._error_type,
        )
        require_canonical_site_series(
            frame.loc[:, "substrate_site"],
            field_name=f"{self._field_name}.substrate_site",
            error_type=self._error_type,
        )
        require_unique_row_pairs(
            frame,
            field_name=self._field_name,
            column_names=("kinase", "substrate_site"),
            error_type=self._error_type,
        )
        return frame


@dataclass(frozen=True, slots=True)
class SiteSequenceReference(TableSchema):
    """Schema wrapper for ``references.site_sequences``."""

    _field_name = "references.site_sequences"
    _error_type = ReferenceValidationError

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
        require_columns(
            frame,
            field_name=self._field_name,
            required_columns=("site_sequence",),
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name="site_sequence",
            error_type=self._error_type,
        )
        require_canonical_string_column(
            frame,
            field_name=self._field_name,
            column_name="site_sequence",
            error_type=self._error_type,
        )
        require_canonical_site_index(
            frame.index,
            field_name=f"{self._field_name}.index",
            error_type=self._error_type,
        )
        return frame
