"""FragPipe source-column resolution."""
# pyright: reportMissingTypeStubs=false

from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.io.readers._table_parsing import (
    resolve_column,
    resolve_intensity_columns,
    resolve_required_column,
)
from phospy.io.readers.fragpipe.constants import (
    _CONTAMINANT_CANDIDATES,
    _DECOY_CANDIDATES,
    _GENE_SYMBOL_CANDIDATES,
    _MODIFIED_PEPTIDE_SEQUENCE_CANDIDATES,
    _PEPTIDE_SEQUENCE_CANDIDATES,
    _PROTEIN_ACCESSION_CANDIDATES,
    _PROTEIN_START_CANDIDATES,
    _PTMPROPHET_PROBABILITY_CANDIDATES,
    _ROW_ID_CANDIDATES,
    _SITE_CANDIDATES,
    _SITE_SEQUENCE_CANDIDATES,
    _UNIQUE_FEATURE_CANDIDATES,
)
from phospy.io.readers.fragpipe.models import (
    FragPipeColumnMapping,
    _ResolvedFragPipeColumns,
)
from phospy.validation.datasets.fragpipe import validate_optional_fragpipe_column_name


def _resolve_fragpipe_columns(
    source: pd.DataFrame,
    mapping: FragPipeColumnMapping,
    *,
    intensity_column_prefixes: Sequence[str],
) -> _ResolvedFragPipeColumns:
    if not isinstance(mapping, FragPipeColumnMapping):
        raise PhosPyInputError(
            "fragpipe import request column_mapping must be a FragPipeColumnMapping"
        )
    columns = pd.Index(source.columns.astype(str).tolist())
    protein_accession = _resolve_required_column(
        columns,
        explicit=mapping.protein_accession,
        candidates=_PROTEIN_ACCESSION_CANDIDATES,
        field_name="fragpipe column_mapping.protein_accession",
    )
    gene_symbol = _resolve_required_column(
        columns,
        explicit=mapping.gene_symbol,
        candidates=_GENE_SYMBOL_CANDIDATES,
        field_name="fragpipe column_mapping.gene_symbol",
    )
    peptide_sequence = _resolve_required_column(
        columns,
        explicit=mapping.peptide_sequence,
        candidates=_PEPTIDE_SEQUENCE_CANDIDATES,
        field_name="fragpipe column_mapping.peptide_sequence",
    )
    modified_peptide_sequence = _resolve_required_column(
        columns,
        explicit=mapping.modified_peptide_sequence,
        candidates=_MODIFIED_PEPTIDE_SEQUENCE_CANDIDATES,
        field_name="fragpipe column_mapping.modified_peptide_sequence",
    )
    ptmprophet_probabilities = _resolve_required_column(
        columns,
        explicit=mapping.ptmprophet_probabilities,
        candidates=_PTMPROPHET_PROBABILITY_CANDIDATES,
        field_name="fragpipe column_mapping.ptmprophet_probabilities",
    )
    intensity_columns = _resolve_intensity_columns(
        source,
        mapping.intensity_columns,
        intensity_column_prefixes=intensity_column_prefixes,
    )
    return _ResolvedFragPipeColumns(
        protein_accession=protein_accession,
        gene_symbol=gene_symbol,
        peptide_sequence=peptide_sequence,
        modified_peptide_sequence=modified_peptide_sequence,
        ptmprophet_probabilities=ptmprophet_probabilities,
        protein_start=_resolve_column(
            columns,
            explicit=mapping.protein_start,
            candidates=_PROTEIN_START_CANDIDATES,
            field_name="fragpipe column_mapping.protein_start",
            required=False,
        ),
        site=_resolve_column(
            columns,
            explicit=mapping.site,
            candidates=_SITE_CANDIDATES,
            field_name="fragpipe column_mapping.site",
            required=False,
        ),
        site_sequence=_resolve_column(
            columns,
            explicit=mapping.site_sequence,
            candidates=_SITE_SEQUENCE_CANDIDATES,
            field_name="fragpipe column_mapping.site_sequence",
            required=False,
        ),
        intensity_columns=intensity_columns,
        contaminant=_resolve_column(
            columns,
            explicit=mapping.contaminant,
            candidates=_CONTAMINANT_CANDIDATES,
            field_name="fragpipe column_mapping.contaminant",
            required=False,
        ),
        decoy=_resolve_column(
            columns,
            explicit=mapping.decoy,
            candidates=_DECOY_CANDIDATES,
            field_name="fragpipe column_mapping.decoy",
            required=False,
        ),
        row_id=_resolve_column(
            columns,
            explicit=mapping.row_id,
            candidates=_ROW_ID_CANDIDATES,
            field_name="fragpipe column_mapping.row_id",
            required=False,
        ),
        unique_feature_id=_resolve_column(
            columns,
            explicit=mapping.unique_feature_id,
            candidates=_UNIQUE_FEATURE_CANDIDATES,
            field_name="fragpipe column_mapping.unique_feature_id",
            required=False,
        ),
    )


def _resolve_required_column(
    columns: pd.Index,
    *,
    explicit: str | None,
    candidates: tuple[str, ...],
    field_name: str,
) -> str:
    return resolve_required_column(
        columns,
        explicit=explicit,
        candidates=candidates,
        field_name=field_name,
        importer_label="FragPipe",
        validate_column_name=validate_optional_fragpipe_column_name,
    )


def _resolve_column(
    columns: pd.Index,
    *,
    explicit: str | None,
    candidates: tuple[str, ...],
    field_name: str,
    required: bool,
) -> str | None:
    return resolve_column(
        columns,
        explicit=explicit,
        candidates=candidates,
        field_name=field_name,
        importer_label="FragPipe",
        required=required,
        validate_column_name=validate_optional_fragpipe_column_name,
    )


def _resolve_intensity_columns(
    source: pd.DataFrame,
    value: Mapping[str, str] | Sequence[str] | None,
    *,
    intensity_column_prefixes: Sequence[str],
) -> dict[str, str]:
    return resolve_intensity_columns(
        source,
        value,
        intensity_column_prefixes=intensity_column_prefixes,
        importer_label="FragPipe",
        request_label="fragpipe",
        mapping_class_name="FragPipeColumnMapping",
        reject_duplicate_inferred_sample_ids=False,
    )


def _resolved_columns_payload(resolved: _ResolvedFragPipeColumns) -> dict[str, object]:
    return {
        "protein_accession": resolved.protein_accession,
        "gene_symbol": resolved.gene_symbol,
        "peptide_sequence": resolved.peptide_sequence,
        "modified_peptide_sequence": resolved.modified_peptide_sequence,
        "ptmprophet_probabilities": resolved.ptmprophet_probabilities,
        "protein_start": resolved.protein_start,
        "site": resolved.site,
        "site_sequence": resolved.site_sequence,
        "intensity_columns": dict(resolved.intensity_columns),
        "contaminant": resolved.contaminant,
        "decoy": resolved.decoy,
        "row_id": resolved.row_id,
        "unique_feature_id": resolved.unique_feature_id,
    }


__all__ = ["_resolve_fragpipe_columns"]
