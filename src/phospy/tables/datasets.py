"""Dataset scientific table wrappers."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from phospy.errors.validation import DatasetValidationError
from phospy.frames.validation import (
    require_dataframe,
    require_exact_index_match,
    require_finite_numeric_dataframe,
    require_non_empty_string_column,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
)
from phospy.science.evidence.localisation import (
    validate_localisation_confidence_column,
    validate_localisation_probability_column,
)
from phospy.science.sites.identifiers import (
    canonicalize_site_components_series,
    canonicalize_site_identifier,
)
from phospy.science.sites.identity_contracts import (
    ANALYSIS_READY_DATASET_BASE_IDENTITY_CONTRACT,
    enforce_analysis_ready_site_key_index,
    enforce_phosphosite_identity_contract,
    enforce_site_key_matches_metadata,
)
from phospy.science.sites.metadata_validation import (
    enforce_site_identity_rows,
    validate_site_identity_metadata,
    validate_site_sequence_column,
)
from phospy.tables.base import (
    TableSchema,
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
        enforce_analysis_ready_site_key_index(
            frame.index,
            field_name=f"{self._field_name}.index",
            error_type=self._error_type,
        )
        return frame


@dataclass(frozen=True, slots=True)
class SiteMetadataTable(TableSchema):
    """Schema wrapper for ``dataset.site_metadata``."""

    expected_index: pd.Index | None = field(default=None, repr=False, compare=False)
    allow_opaque_site_values: bool = field(default=False, repr=False, compare=False)

    _field_name = "dataset.site_metadata"
    _error_type = DatasetValidationError

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        require_dataframe(
            frame,
            field_name=self._field_name,
            allow_empty=False,
            error_type=self._error_type,
        )
        enforce_analysis_ready_site_key_index(
            frame.index,
            field_name=f"{self._field_name}.index",
            error_type=self._error_type,
        )
        enforce_phosphosite_identity_contract(
            site_metadata=frame,
            field_name=self._field_name,
            contract=ANALYSIS_READY_DATASET_BASE_IDENTITY_CONTRACT,
            error_type=self._error_type,
            expected_index=self.expected_index,
            expected_index_field_name="dataset.phospho.index",
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name="site_key",
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name="display_id",
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name="organism",
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name="protein_namespace",
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name="protein_identifier",
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
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name="site_sequence",
            error_type=self._error_type,
        )
        validate_site_sequence_column(
            site_metadata=frame,
            field_name=self._field_name,
            error_type=self._error_type,
            column_name="site_sequence",
        )
        validate_localisation_probability_column(
            site_metadata=frame,
            field_name=self._field_name,
            error_type=self._error_type,
            column_name="localisation_probability",
        )
        validate_localisation_confidence_column(
            site_metadata=frame,
            field_name=self._field_name,
            error_type=self._error_type,
            column_name="localisation_confidence",
        )
        base_identity_frame = _drop_signalome_grouping_metadata(frame)
        validate_site_identity_metadata(
            site_metadata=base_identity_frame,
            field_name=self._field_name,
            error_type=self._error_type,
            allow_opaque_site_values=self.allow_opaque_site_values,
        )
        enforce_site_key_matches_metadata(
            site_metadata=frame,
            field_name=self._field_name,
            error_type=self._error_type,
            site_key_column="site_key",
        )
        identity_frame = _build_identity_coherence_frame(frame)
        try:
            enforce_site_identity_rows(
                site_metadata=identity_frame,
                field_name=self._field_name,
                error_type=self._error_type,
                allow_opaque_site_values=self.allow_opaque_site_values,
            )
        except self._error_type as exc:
            raise self._error_type(
                f"{self._field_name} site-identity coherence failed; {exc}"
            ) from exc
        return frame


@dataclass(frozen=True, slots=True)
class SampleMetadataTable(TableSchema):
    """Schema wrapper for ``dataset.sample_metadata``.

    This table validates metadata-table integrity and sample-index alignment only.
    It does not validate workflow experimental design semantics.
    """

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
        require_unique_columns(
            frame,
            field_name=self._field_name,
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


def _display_ids_are_canonical(frame: pd.DataFrame) -> bool:
    values = frame.loc[:, "display_id"].astype(str).tolist()
    for value in values:
        try:
            canonicalize_site_identifier(
                value,
                field_name="dataset.site_metadata.display_id",
                error_type=DatasetValidationError,
            )
        except DatasetValidationError:
            return False
    return True


def _build_identity_coherence_frame(frame: pd.DataFrame) -> pd.DataFrame:
    identity_frame = _drop_signalome_grouping_metadata(frame)
    if "display_id" not in identity_frame.columns:
        return identity_frame
    if _display_ids_are_canonical(identity_frame):
        reindexed = identity_frame.copy(deep=False)
        reindexed.index = pd.Index(
            reindexed.loc[:, "display_id"].astype(str).tolist(),
            name="display_id",
        )
        return reindexed
    try:
        canonical_display = canonicalize_site_components_series(
            gene_symbol=identity_frame.loc[:, "gene_symbol"],
            site=identity_frame.loc[:, "site"],
            field_name="dataset.site_metadata.gene_symbol/site",
            error_type=DatasetValidationError,
            output_name="display_id",
        )
        with_canonical_display = identity_frame.copy(deep=True)
        with_canonical_display.loc[:, "display_id"] = canonical_display.astype(
            str
        ).tolist()
        with_canonical_display.index = pd.Index(
            with_canonical_display.loc[:, "display_id"].astype(str).tolist(),
            name="display_id",
        )
        return with_canonical_display
    except DatasetValidationError:
        without_display = identity_frame.copy(deep=False)
        return without_display.drop(columns=["display_id"], errors="ignore")


def _drop_signalome_grouping_metadata(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.drop(columns=["protein_id"], errors="ignore")
