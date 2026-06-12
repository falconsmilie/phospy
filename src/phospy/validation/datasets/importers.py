"""Validation helpers for upstream phosphosite importers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import cast

import pandas as pd

from phospy.contracts.requests import PhosphositeImportRequest
from phospy.errors.input import PhosPyInputError, UnsupportedInputFormatError
from phospy.science.evidence.localisation import (
    SUPPORTED_LOCALISATION_CONFIDENCE_SCALES,
)
from phospy.validation.common.dataframes import require_unique_columns
from phospy.validation.datasets.inputs import DatasetInputSourceValidator


class PhosphositeImportRequestValidator:
    """Validate importer request structure before source-table translation."""

    def __init__(
        self,
        *,
        source_validator: DatasetInputSourceValidator | None = None,
    ) -> None:
        self._source_validator = source_validator or DatasetInputSourceValidator()

    def run(self, request: object) -> PhosphositeImportRequest:
        if not isinstance(request, PhosphositeImportRequest):
            raise PhosPyInputError(
                "phosphosite importer input must be a PhosphositeImportRequest"
            )
        self._source_validator.run(request.source, field_name="source")
        normalise_sample_column_mapping(request.sample_intensity_columns)
        for field_name, value in _request_column_fields(request).items():
            _validate_optional_column_name(value, field_name=field_name)
        _validate_required_column_name(
            request.gene_symbol_column,
            field_name="phosphosite import request gene_symbol_column",
        )
        _validate_required_column_name(
            request.site_column,
            field_name="phosphosite import request site_column",
        )
        _validate_localisation_scale(request.localisation_confidence_scale)
        _validate_source_name(request.source_name)
        _validate_peptide_evidence_column_group(request)
        return request


def normalise_sample_column_mapping(
    value: object,
) -> dict[str, str]:
    """Return explicit ``source_column -> sample_id`` mapping."""

    raw_mapping: dict[object, object]
    if isinstance(value, str):
        raise PhosPyInputError(
            "phosphosite import request sample_intensity_columns must be a mapping "
            "of source column to sample ID or a sequence of source column names"
        )
    if isinstance(value, Mapping):
        raw_mapping = dict(cast(Mapping[object, object], value))
    elif not isinstance(value, Sequence):
        raise PhosPyInputError(
            "phosphosite import request sample_intensity_columns must be a mapping "
            "of source column to sample ID or a sequence of source column names"
        )
    else:
        raw_mapping = {
            column: column for column in tuple(cast(Sequence[object], value))
        }
    if not raw_mapping:
        raise PhosPyInputError(
            "phosphosite import request sample_intensity_columns must not be empty"
        )

    mapping: dict[str, str] = {}
    for source_column, sample_id in raw_mapping.items():
        if not isinstance(source_column, str) or source_column.strip() == "":
            raise PhosPyInputError(
                "phosphosite import request sample intensity source columns must "
                "be non-empty strings"
            )
        if not isinstance(sample_id, str) or sample_id.strip() == "":
            raise PhosPyInputError(
                "phosphosite import request sample intensity sample IDs must be "
                "non-empty strings"
            )
        mapping[source_column.strip()] = sample_id.strip()
    if len(mapping) != len(raw_mapping):
        raise PhosPyInputError(
            "phosphosite import request sample_intensity_columns must not contain "
            "duplicate source column names after trimming"
        )
    if len(set(mapping.values())) != len(mapping):
        raise PhosPyInputError(
            "phosphosite import request sample_intensity_columns must map to "
            "unique sample IDs"
        )
    return mapping


def required_import_source_columns(
    request: PhosphositeImportRequest,
    sample_column_mapping: dict[str, str],
) -> tuple[str, ...]:
    """Return source columns that must be present for the mapped import."""

    columns = [
        request.gene_symbol_column,
        request.site_column,
        *sample_column_mapping.keys(),
    ]
    optional_columns = _request_column_fields(request)
    for value in optional_columns.values():
        if value is not None:
            columns.append(value)
    columns.extend(_required_peptide_evidence_columns(request))
    return tuple(dict.fromkeys(columns))


def require_import_source_columns(
    frame: pd.DataFrame,
    *,
    required_columns: tuple[str, ...],
    field_name: str = "phosphosite import source",
) -> None:
    require_unique_columns(frame, field_name=field_name, error_type=PhosPyInputError)
    missing = [column for column in required_columns if column not in frame.columns]
    if not missing:
        return
    joined = ", ".join(missing)
    raise PhosPyInputError(
        f"{field_name} is missing required columns: {joined}; "
        f"missing_column_count={int(len(missing))}, "
        f"column_count={int(frame.shape[1])}"
    )


def peptide_evidence_requested(request: PhosphositeImportRequest) -> bool:
    """Return whether request fields ask the importer to emit peptide evidence."""

    evidence_fields = (
        request.unique_feature_id_column,
        request.peptide_sequence_column,
        request.modified_peptide_sequence_column,
        request.peptide_site_string_column,
        request.peptide_site_id_column,
        request.peptide_row_id_column,
    )
    return any(value is not None for value in evidence_fields)


def _request_column_fields(request: PhosphositeImportRequest) -> dict[str, str | None]:
    return {
        "phosphosite import request row_id_column": request.row_id_column,
        "phosphosite import request protein_id_column": request.protein_id_column,
        "phosphosite import request protein_accession_column": (
            request.protein_accession_column
        ),
        "phosphosite import request protein_identifier_column": (
            request.protein_identifier_column
        ),
        "phosphosite import request protein_namespace_column": (
            request.protein_namespace_column
        ),
        "phosphosite import request organism_column": request.organism_column,
        "phosphosite import request isoform_id_column": request.isoform_id_column,
        "phosphosite import request site_sequence_column": (
            request.site_sequence_column
        ),
        "phosphosite import request display_id_column": request.display_id_column,
        "phosphosite import request site_key_column": request.site_key_column,
        "phosphosite import request localisation_confidence_column": (
            request.localisation_confidence_column
        ),
        "phosphosite import request peptide_row_id_column": (
            request.peptide_row_id_column
        ),
        "phosphosite import request unique_feature_id_column": (
            request.unique_feature_id_column
        ),
        "phosphosite import request peptide_sequence_column": (
            request.peptide_sequence_column
        ),
        "phosphosite import request modified_peptide_sequence_column": (
            request.modified_peptide_sequence_column
        ),
        "phosphosite import request peptide_site_string_column": (
            request.peptide_site_string_column
        ),
        "phosphosite import request peptide_site_id_column": (
            request.peptide_site_id_column
        ),
    }


def _validate_peptide_evidence_column_group(
    request: PhosphositeImportRequest,
) -> None:
    if not peptide_evidence_requested(request):
        return
    missing_fields = [
        field_name
        for field_name, value in {
            "unique_feature_id_column": request.unique_feature_id_column,
            "peptide_sequence_column": request.peptide_sequence_column,
            "modified_peptide_sequence_column": (
                request.modified_peptide_sequence_column
            ),
            "peptide_site_string_column": request.peptide_site_string_column,
        }.items()
        if value is None
    ]
    if missing_fields:
        joined = ", ".join(missing_fields)
        raise PhosPyInputError(
            "phosphosite import request peptide evidence mapping is incomplete; "
            f"missing fields: {joined}"
        )
    if (
        request.protein_accession_column is None
        and request.protein_id_column is None
        and request.protein_identifier_column is None
    ):
        raise PhosPyInputError(
            "phosphosite import request peptide evidence mapping requires protein "
            "context via protein_accession_column, protein_id_column, or "
            "protein_identifier_column"
        )


def _required_peptide_evidence_columns(
    request: PhosphositeImportRequest,
) -> tuple[str, ...]:
    if not peptide_evidence_requested(request):
        return ()
    values = (
        request.unique_feature_id_column,
        request.peptide_sequence_column,
        request.modified_peptide_sequence_column,
        request.peptide_site_string_column,
    )
    return tuple(value for value in values if value is not None)


def _validate_optional_column_name(value: str | None, *, field_name: str) -> None:
    if value is None:
        return
    _validate_required_column_name(value, field_name=field_name)


def _validate_required_column_name(value: object, *, field_name: str) -> None:
    if not isinstance(value, str) or value.strip() == "":
        raise PhosPyInputError(f"{field_name} must be a non-empty string")


def _validate_localisation_scale(value: object) -> None:
    if (
        isinstance(value, str)
        and value.strip() in SUPPORTED_LOCALISATION_CONFIDENCE_SCALES
    ):
        return
    supported = ", ".join(
        repr(scale) for scale in SUPPORTED_LOCALISATION_CONFIDENCE_SCALES
    )
    raise PhosPyInputError(
        "phosphosite import request localisation_confidence_scale must be one of: "
        f"{supported}"
    )


def _validate_source_name(value: object) -> None:
    if not isinstance(value, str) or value.strip() == "":
        raise UnsupportedInputFormatError(
            "phosphosite import request source_name must be a non-empty string"
        )


__all__ = [
    "PhosphositeImportRequestValidator",
    "normalise_sample_column_mapping",
    "peptide_evidence_requested",
    "required_import_source_columns",
    "require_import_source_columns",
]
