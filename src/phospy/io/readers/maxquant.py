"""MaxQuant phosphosite table importer."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace

import pandas as pd

from phospy.contracts.requests import PhosphositeImportRequest
from phospy.contracts.results import (
    IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
    IMPORTER_QUALITY_STATUS_REPORTED,
    ImporterFlaggedRowSummary,
    ImporterQualityCount,
    ImporterQualityReport,
    PhosphositeImportResult,
)
from phospy.errors.input import PhosPyInputError
from phospy.io.readers._table_parsing import (
    build_row_ids,
    build_unique_feature_ids,
    first_list_token,
    multi_value_count,
    optional_text,
    raise_for_forbidden_flags,
    require_non_empty_unique_columns,
    required_text,
    resolve_column,
    resolve_flag_series,
    resolve_intensity_columns,
    resolve_required_column,
    split_multi_value,
)
from phospy.io.readers.importers import (
    MappedPhosphositeTableImporter,
    _read_upstream_table,
)
from phospy.science.evidence.multi_site import parse_phospho_site_tokens
from phospy.validation.datasets.maxquant import (
    MAXQUANT_FLAG_POLICY_ERROR,
    MAXQUANT_FLAG_POLICY_FLAG,
    MAXQUANT_FLAG_POLICY_REMOVE,
    validate_maxquant_flag_policy,
    validate_optional_maxquant_column_name,
)

_ADAPTED_ROW_ID_COLUMN = "__phospy_maxquant_row_id"
_ADAPTED_PROTEIN_ACCESSION_COLUMN = "__phospy_maxquant_protein_accession"
_ADAPTED_PROTEIN_ID_COLUMN = "__phospy_maxquant_protein_id"
_ADAPTED_GENE_SYMBOL_COLUMN = "__phospy_maxquant_gene_symbol"
_ADAPTED_SITE_COLUMN = "__phospy_maxquant_site"
_ADAPTED_SITE_SEQUENCE_COLUMN = "__phospy_maxquant_site_sequence"
_ADAPTED_LOCALISATION_COLUMN = "__phospy_maxquant_localisation_confidence"
_ADAPTED_UNIQUE_FEATURE_ID_COLUMN = "__phospy_maxquant_feature_id"
_ADAPTED_PEPTIDE_SEQUENCE_COLUMN = "__phospy_maxquant_peptide_sequence"
_ADAPTED_MODIFIED_PEPTIDE_SEQUENCE_COLUMN = (
    "__phospy_maxquant_modified_peptide_sequence"
)
_ADAPTED_PEPTIDE_SITE_STRING_COLUMN = "__phospy_maxquant_peptide_site_string"
_MAXQUANT_CONTAMINANT_OUTPUT_COLUMN = "maxquant_potential_contaminant"
_MAXQUANT_REVERSE_OUTPUT_COLUMN = "maxquant_reverse"
_DEFAULT_INTENSITY_PREFIXES = (
    "Intensity ",
    "LFQ intensity ",
    "Reporter intensity corrected ",
)
_PROTEIN_ACCESSION_CANDIDATES = (
    "Leading proteins",
    "Leading protein",
    "Proteins",
    "Protein",
    "Protein IDs",
    "Majority protein IDs",
)
_GENE_SYMBOL_CANDIDATES = (
    "Gene names",
    "Gene name",
    "Genes",
    "Gene",
)
_SITE_CANDIDATES = (
    "Modified site",
    "Phosphosite",
    "Phospho site",
    "Site",
)
_AMINO_ACID_CANDIDATES = (
    "Amino acid",
    "Amino acids",
    "Modified amino acid",
)
_POSITION_CANDIDATES = (
    "Positions within proteins",
    "Position within protein",
    "Position",
    "Positions",
)
_LOCALISATION_CANDIDATES = (
    "Localization prob",
    "Localization probability",
    "Localisation prob",
    "Localisation probability",
    "Phospho (STY) Probabilities",
)
_PEPTIDE_SEQUENCE_CANDIDATES = (
    "Sequence",
    "Peptide sequence",
)
_MODIFIED_PEPTIDE_SEQUENCE_CANDIDATES = (
    "Modified sequence",
    "Modified peptide sequence",
)
_SITE_SEQUENCE_CANDIDATES = (
    "Sequence window",
    "Sequence Window",
    "Window sequence",
)
_UNIQUE_FEATURE_CANDIDATES = (
    "id",
    "ID",
    "Phospho (STY) site IDs",
    "Site IDs",
    "Evidence IDs",
)
_ROW_ID_CANDIDATES: tuple[str, ...] = ()
_POTENTIAL_CONTAMINANT_CANDIDATES = (
    "Potential contaminant",
    "Potential contaminant?",
)
_REVERSE_CANDIDATES = (
    "Reverse",
    "Reverse?",
)
_LOCALISATION_PROBABILITY_TOKEN_PATTERN = re.compile(
    r"[STYsty]\s*\(\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*\)"
)


@dataclass(frozen=True, slots=True)
class MaxQuantColumnMapping:
    """Optional source-column overrides for MaxQuant phosphosite imports."""

    protein_accession: str | None = None
    gene_symbol: str | None = None
    modified_site: str | None = None
    amino_acid: str | None = None
    site_position: str | None = None
    localisation_confidence: str | None = None
    peptide_sequence: str | None = None
    modified_peptide_sequence: str | None = None
    intensity_columns: Mapping[str, str] | Sequence[str] | None = None
    potential_contaminant: str | None = None
    reverse: str | None = None
    row_id: str | None = None
    unique_feature_id: str | None = None
    site_sequence: str | None = None


@dataclass(frozen=True, slots=True)
class MaxQuantPhosphositeImportRequest:
    """Request for importing MaxQuant-style phosphosite output.

    Localisation confidence is emitted as ``localisation_confidence`` on the
    common importer result. Values are normalised by the shared importer to
    probabilities in ``[0.0, 1.0]`` according to
    ``localisation_confidence_scale``.
    """

    source: object
    column_mapping: MaxQuantColumnMapping = field(default_factory=MaxQuantColumnMapping)
    contaminant_policy: str = MAXQUANT_FLAG_POLICY_REMOVE
    reverse_policy: str = MAXQUANT_FLAG_POLICY_REMOVE
    localisation_confidence_scale: str = "probability"
    intensity_column_prefixes: Sequence[str] = _DEFAULT_INTENSITY_PREFIXES
    source_name: str = "maxquant"


@dataclass(frozen=True, slots=True)
class _ResolvedMaxQuantColumns:
    protein_accession: str
    gene_symbol: str
    modified_site: str | None
    amino_acid: str | None
    site_position: str | None
    localisation_confidence: str
    peptide_sequence: str
    modified_peptide_sequence: str | None
    intensity_columns: dict[str, str]
    potential_contaminant: str | None
    reverse: str | None
    row_id: str | None
    unique_feature_id: str | None
    site_sequence: str | None


class MaxQuantPhosphositeImporter:
    """Import MaxQuant-style phosphosite output into common PhosPy candidates."""

    def __init__(
        self,
        *,
        mapped_importer: MappedPhosphositeTableImporter | None = None,
    ) -> None:
        self._mapped_importer = mapped_importer or MappedPhosphositeTableImporter()

    def run(
        self,
        request: MaxQuantPhosphositeImportRequest,
    ) -> PhosphositeImportResult:
        if not isinstance(request, MaxQuantPhosphositeImportRequest):
            raise PhosPyInputError(
                "MaxQuant importer input must be a MaxQuantPhosphositeImportRequest"
            )
        contaminant_policy = validate_maxquant_flag_policy(
            request.contaminant_policy,
            field_name="maxquant import request contaminant_policy",
        )
        reverse_policy = validate_maxquant_flag_policy(
            request.reverse_policy,
            field_name="maxquant import request reverse_policy",
        )
        source = _read_upstream_table(request.source)
        _require_non_empty_unique_columns(source)
        resolved = _resolve_maxquant_columns(
            source,
            request.column_mapping,
            intensity_column_prefixes=request.intensity_column_prefixes,
        )
        filtered, flags, filter_diagnostics, filter_warnings = _apply_flag_policies(
            source,
            resolved=resolved,
            contaminant_policy=contaminant_policy,
            reverse_policy=reverse_policy,
        )
        if filtered.empty:
            raise PhosPyInputError(
                "MaxQuant importer removed all rows after contaminant/reverse filtering"
            )

        adapted, adapter_diagnostics, adapter_warnings = _adapt_maxquant_source(
            filtered,
            resolved=resolved,
        )
        mapped_result = self._mapped_importer.run(
            PhosphositeImportRequest(
                source=adapted,
                sample_intensity_columns=resolved.intensity_columns,
                gene_symbol_column=_ADAPTED_GENE_SYMBOL_COLUMN,
                site_column=_ADAPTED_SITE_COLUMN,
                row_id_column=_ADAPTED_ROW_ID_COLUMN,
                protein_id_column=_ADAPTED_PROTEIN_ID_COLUMN,
                protein_accession_column=_ADAPTED_PROTEIN_ACCESSION_COLUMN,
                site_sequence_column=(
                    _ADAPTED_SITE_SEQUENCE_COLUMN
                    if _ADAPTED_SITE_SEQUENCE_COLUMN in adapted.columns
                    else None
                ),
                localisation_confidence_column=_ADAPTED_LOCALISATION_COLUMN,
                localisation_confidence_scale=request.localisation_confidence_scale,
                unique_feature_id_column=_ADAPTED_UNIQUE_FEATURE_ID_COLUMN,
                peptide_sequence_column=_ADAPTED_PEPTIDE_SEQUENCE_COLUMN,
                modified_peptide_sequence_column=(
                    _ADAPTED_MODIFIED_PEPTIDE_SEQUENCE_COLUMN
                ),
                peptide_site_string_column=_ADAPTED_PEPTIDE_SITE_STRING_COLUMN,
                source_name=request.source_name,
            )
        )
        return _augment_mapped_result(
            mapped_result,
            adapted=adapted,
            flags=flags,
            contaminant_policy=contaminant_policy,
            reverse_policy=reverse_policy,
            resolved=resolved,
            filter_diagnostics=filter_diagnostics,
            adapter_diagnostics=adapter_diagnostics,
            warnings=filter_warnings + adapter_warnings,
        )


def _resolve_maxquant_columns(
    source: pd.DataFrame,
    mapping: MaxQuantColumnMapping,
    *,
    intensity_column_prefixes: Sequence[str],
) -> _ResolvedMaxQuantColumns:
    if not isinstance(mapping, MaxQuantColumnMapping):
        raise PhosPyInputError(
            "maxquant import request column_mapping must be a MaxQuantColumnMapping"
        )
    columns = pd.Index(source.columns.astype(str).tolist())
    protein_accession = _resolve_required_column(
        columns,
        explicit=mapping.protein_accession,
        candidates=_PROTEIN_ACCESSION_CANDIDATES,
        field_name="maxquant column_mapping.protein_accession",
    )
    gene_symbol = _resolve_required_column(
        columns,
        explicit=mapping.gene_symbol,
        candidates=_GENE_SYMBOL_CANDIDATES,
        field_name="maxquant column_mapping.gene_symbol",
    )
    modified_site = _resolve_column(
        columns,
        explicit=mapping.modified_site,
        candidates=_SITE_CANDIDATES,
        field_name="maxquant column_mapping.modified_site",
        required=False,
    )
    amino_acid = _resolve_column(
        columns,
        explicit=mapping.amino_acid,
        candidates=_AMINO_ACID_CANDIDATES,
        field_name="maxquant column_mapping.amino_acid",
        required=False,
    )
    site_position = _resolve_column(
        columns,
        explicit=mapping.site_position,
        candidates=_POSITION_CANDIDATES,
        field_name="maxquant column_mapping.site_position",
        required=False,
    )
    if modified_site is None and (amino_acid is None or site_position is None):
        raise PhosPyInputError(
            "MaxQuant importer requires either a modified_site column containing "
            "tokens like 'S123' or both amino_acid and site_position columns"
        )
    localisation_confidence = _resolve_required_column(
        columns,
        explicit=mapping.localisation_confidence,
        candidates=_LOCALISATION_CANDIDATES,
        field_name="maxquant column_mapping.localisation_confidence",
    )
    peptide_sequence = _resolve_required_column(
        columns,
        explicit=mapping.peptide_sequence,
        candidates=_PEPTIDE_SEQUENCE_CANDIDATES,
        field_name="maxquant column_mapping.peptide_sequence",
    )
    modified_peptide_sequence = _resolve_column(
        columns,
        explicit=mapping.modified_peptide_sequence,
        candidates=_MODIFIED_PEPTIDE_SEQUENCE_CANDIDATES,
        field_name="maxquant column_mapping.modified_peptide_sequence",
        required=False,
    )
    intensity_columns = _resolve_intensity_columns(
        source,
        mapping.intensity_columns,
        intensity_column_prefixes=intensity_column_prefixes,
    )
    return _ResolvedMaxQuantColumns(
        protein_accession=protein_accession,
        gene_symbol=gene_symbol,
        modified_site=modified_site,
        amino_acid=amino_acid,
        site_position=site_position,
        localisation_confidence=localisation_confidence,
        peptide_sequence=peptide_sequence,
        modified_peptide_sequence=modified_peptide_sequence,
        intensity_columns=intensity_columns,
        potential_contaminant=_resolve_column(
            columns,
            explicit=mapping.potential_contaminant,
            candidates=_POTENTIAL_CONTAMINANT_CANDIDATES,
            field_name="maxquant column_mapping.potential_contaminant",
            required=False,
        ),
        reverse=_resolve_column(
            columns,
            explicit=mapping.reverse,
            candidates=_REVERSE_CANDIDATES,
            field_name="maxquant column_mapping.reverse",
            required=False,
        ),
        row_id=_resolve_column(
            columns,
            explicit=mapping.row_id,
            candidates=_ROW_ID_CANDIDATES,
            field_name="maxquant column_mapping.row_id",
            required=False,
        ),
        unique_feature_id=_resolve_column(
            columns,
            explicit=mapping.unique_feature_id,
            candidates=_UNIQUE_FEATURE_CANDIDATES,
            field_name="maxquant column_mapping.unique_feature_id",
            required=False,
        ),
        site_sequence=_resolve_column(
            columns,
            explicit=mapping.site_sequence,
            candidates=_SITE_SEQUENCE_CANDIDATES,
            field_name="maxquant column_mapping.site_sequence",
            required=False,
        ),
    )


def _apply_flag_policies(
    source: pd.DataFrame,
    *,
    resolved: _ResolvedMaxQuantColumns,
    contaminant_policy: str,
    reverse_policy: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object], tuple[str, ...]]:
    flags = pd.DataFrame(index=source.index.copy())
    contaminant_flags = resolve_flag_series(
        source,
        column=resolved.potential_contaminant,
        field_name="MaxQuant potential contaminant flag",
    )
    reverse_flags = resolve_flag_series(
        source,
        column=resolved.reverse,
        field_name="MaxQuant reverse flag",
    )
    flags[_MAXQUANT_CONTAMINANT_OUTPUT_COLUMN] = (
        [False] * int(source.shape[0])
        if contaminant_flags is None
        else contaminant_flags.tolist()
    )
    flags[_MAXQUANT_REVERSE_OUTPUT_COLUMN] = (
        [False] * int(source.shape[0])
        if reverse_flags is None
        else reverse_flags.tolist()
    )

    _raise_for_forbidden_flags(
        flags[_MAXQUANT_CONTAMINANT_OUTPUT_COLUMN],
        policy=contaminant_policy,
        label="potential contaminant",
    )
    _raise_for_forbidden_flags(
        flags[_MAXQUANT_REVERSE_OUTPUT_COLUMN],
        policy=reverse_policy,
        label="reverse",
    )
    keep_mask = pd.Series(True, index=source.index.copy(), dtype=bool)
    if contaminant_policy == MAXQUANT_FLAG_POLICY_REMOVE:
        keep_mask &= ~flags[_MAXQUANT_CONTAMINANT_OUTPUT_COLUMN].astype(bool)
    if reverse_policy == MAXQUANT_FLAG_POLICY_REMOVE:
        keep_mask &= ~flags[_MAXQUANT_REVERSE_OUTPUT_COLUMN].astype(bool)

    diagnostics = {
        "input_row_count": int(source.shape[0]),
        "potential_contaminant_column": resolved.potential_contaminant,
        "reverse_column": resolved.reverse,
        "potential_contaminant_policy": contaminant_policy,
        "reverse_policy": reverse_policy,
        "potential_contaminant_rows": int(
            flags[_MAXQUANT_CONTAMINANT_OUTPUT_COLUMN].astype(bool).sum()
        ),
        "reverse_rows": int(flags[_MAXQUANT_REVERSE_OUTPUT_COLUMN].astype(bool).sum()),
        "removed_rows": int((~keep_mask).sum()),
        "retained_row_count": int(keep_mask.sum()),
    }
    warnings: list[str] = []
    if resolved.potential_contaminant is None:
        warnings.append(
            "MaxQuant potential contaminant column was not found; contaminant "
            "policy could not filter or flag rows"
        )
    if resolved.reverse is None:
        warnings.append(
            "MaxQuant reverse column was not found; reverse policy could not "
            "filter or flag rows"
        )
    filtered = source.loc[keep_mask, :].copy(deep=True)
    filtered_flags = flags.loc[keep_mask, :].copy(deep=True)
    return filtered, filtered_flags, diagnostics, tuple(warnings)


def _adapt_maxquant_source(
    source: pd.DataFrame,
    *,
    resolved: _ResolvedMaxQuantColumns,
) -> tuple[pd.DataFrame, dict[str, object], tuple[str, ...]]:
    adapted = source.copy(deep=True)
    source_row_numbers = [position + 1 for position in range(int(source.shape[0]))]
    protein_values = [
        first_list_token(
            value,
            field_name=f"MaxQuant {resolved.protein_accession}",
            row_position=position,
        )
        for position, value in enumerate(source.loc[:, resolved.protein_accession])
    ]
    gene_values = [
        first_list_token(
            value,
            field_name=f"MaxQuant {resolved.gene_symbol}",
            row_position=position,
        )
        for position, value in enumerate(source.loc[:, resolved.gene_symbol])
    ]
    gene_group_rows = int(
        sum(
            1
            for value in source.loc[:, resolved.gene_symbol]
            if multi_value_count(value) > 1
        )
    )
    site_values: list[str] = []
    peptide_site_values: list[str] = []
    protein_group_rows = 0
    for position, row in enumerate(source.itertuples(index=False, name=None)):
        row_lookup: dict[str, object] = dict(
            zip(source.columns.astype(str).tolist(), row, strict=True)
        )
        if multi_value_count(row_lookup[resolved.protein_accession]) > 1:
            protein_group_rows += 1
        site_tokens = _resolve_row_site_tokens(
            row_lookup,
            resolved=resolved,
            row_position=position,
        )
        site_values.append(",".join(site_tokens))
        peptide_site_values.append(";".join(site_tokens))

    adapted[_ADAPTED_ROW_ID_COLUMN] = _build_row_ids(
        source=source,
        resolved=resolved,
        protein_values=protein_values,
        site_values=site_values,
        source_row_numbers=source_row_numbers,
    )
    adapted[_ADAPTED_PROTEIN_ACCESSION_COLUMN] = protein_values
    adapted[_ADAPTED_PROTEIN_ID_COLUMN] = protein_values
    adapted[_ADAPTED_GENE_SYMBOL_COLUMN] = gene_values
    adapted[_ADAPTED_SITE_COLUMN] = site_values
    adapted[_ADAPTED_LOCALISATION_COLUMN] = _normalise_localisation_source_values(
        source.loc[:, resolved.localisation_confidence]
    )
    adapted[_ADAPTED_UNIQUE_FEATURE_ID_COLUMN] = _build_unique_feature_ids(
        source=source,
        resolved=resolved,
        source_row_numbers=source_row_numbers,
    )
    adapted[_ADAPTED_PEPTIDE_SEQUENCE_COLUMN] = [
        required_text(
            value,
            field_name=f"MaxQuant {resolved.peptide_sequence}",
            row_position=position,
        )
        for position, value in enumerate(source.loc[:, resolved.peptide_sequence])
    ]
    adapted[_ADAPTED_MODIFIED_PEPTIDE_SEQUENCE_COLUMN] = (
        _resolve_modified_peptide_sequences(source, resolved=resolved)
    )
    adapted[_ADAPTED_PEPTIDE_SITE_STRING_COLUMN] = peptide_site_values
    if resolved.site_sequence is not None:
        adapted[_ADAPTED_SITE_SEQUENCE_COLUMN] = [
            optional_text(value)
            for value in source.loc[:, resolved.site_sequence].tolist()
        ]

    diagnostics = {
        "resolved_columns": _resolved_columns_payload(resolved),
        "protein_group_rows_collapsed_to_first_accession": int(protein_group_rows),
        "gene_group_rows_collapsed_to_first_symbol": int(gene_group_rows),
        "multi_site_rows": int(
            sum(1 for site_value in site_values if len(site_value.split(",")) > 1)
        ),
    }
    warnings: list[str] = []
    if resolved.modified_peptide_sequence is None:
        warnings.append(
            "MaxQuant modified peptide sequence column was not found; peptide "
            "sequence was reused as modified_peptide_sequence"
        )
    if protein_group_rows:
        warnings.append(
            "MaxQuant protein-group rows were represented by the first listed "
            "protein accession for protein-scoped identity"
        )
    if gene_group_rows:
        warnings.append(
            "MaxQuant gene-name group rows were represented by the first listed "
            "gene symbol for display metadata"
        )
    return adapted, diagnostics, tuple(warnings)


def _augment_mapped_result(
    mapped_result: PhosphositeImportResult,
    *,
    adapted: pd.DataFrame,
    flags: pd.DataFrame,
    contaminant_policy: str,
    reverse_policy: str,
    resolved: _ResolvedMaxQuantColumns,
    filter_diagnostics: dict[str, object],
    adapter_diagnostics: dict[str, object],
    warnings: tuple[str, ...],
) -> PhosphositeImportResult:
    site_metadata = mapped_result.site_metadata_candidate
    peptide_evidence = mapped_result.peptide_evidence
    flag_values = flags.copy(deep=True)
    flag_values.index = pd.Index(
        adapted.loc[:, _ADAPTED_ROW_ID_COLUMN].astype(str).tolist(),
        name=_ADAPTED_ROW_ID_COLUMN,
    )
    if contaminant_policy == MAXQUANT_FLAG_POLICY_FLAG:
        site_metadata[_MAXQUANT_CONTAMINANT_OUTPUT_COLUMN] = (
            flag_values.loc[site_metadata.index, _MAXQUANT_CONTAMINANT_OUTPUT_COLUMN]
            .astype(bool)
            .tolist()
        )
        if peptide_evidence is not None:
            peptide_evidence[_MAXQUANT_CONTAMINANT_OUTPUT_COLUMN] = (
                flag_values.loc[
                    peptide_evidence.index, _MAXQUANT_CONTAMINANT_OUTPUT_COLUMN
                ]
                .astype(bool)
                .tolist()
            )
    if reverse_policy == MAXQUANT_FLAG_POLICY_FLAG:
        site_metadata[_MAXQUANT_REVERSE_OUTPUT_COLUMN] = (
            flag_values.loc[site_metadata.index, _MAXQUANT_REVERSE_OUTPUT_COLUMN]
            .astype(bool)
            .tolist()
        )
        if peptide_evidence is not None:
            peptide_evidence[_MAXQUANT_REVERSE_OUTPUT_COLUMN] = (
                flag_values.loc[peptide_evidence.index, _MAXQUANT_REVERSE_OUTPUT_COLUMN]
                .astype(bool)
                .tolist()
            )

    diagnostics = dict(mapped_result.diagnostics)
    diagnostics["maxquant"] = {
        "source_type": "maxquant_phosphosite",
        "resolved_columns": _resolved_columns_payload(resolved),
        "filtering": filter_diagnostics,
        "adaptation": adapter_diagnostics,
    }
    combined_warnings = tuple(dict.fromkeys((*mapped_result.warnings, *warnings)))
    quality_report = _augment_quality_report(
        mapped_result.quality_report,
        resolved=resolved,
        filter_diagnostics=filter_diagnostics,
        adapter_diagnostics=adapter_diagnostics,
        contaminant_policy=contaminant_policy,
        reverse_policy=reverse_policy,
        warnings=combined_warnings,
    )
    return PhosphositeImportResult(
        phospho_matrix_candidate=mapped_result.phospho_matrix_candidate,
        site_metadata_candidate=site_metadata,
        peptide_evidence=peptide_evidence,
        sample_column_mapping=mapped_result.sample_column_mapping,
        localisation_confidence_column=mapped_result.localisation_confidence_column,
        warnings=combined_warnings,
        diagnostics=diagnostics,
        source_name=mapped_result.source_name,
        quality_report=quality_report,
        _assume_owned=True,
    )


def _augment_quality_report(
    mapped_report: ImporterQualityReport,
    *,
    resolved: _ResolvedMaxQuantColumns,
    filter_diagnostics: dict[str, object],
    adapter_diagnostics: dict[str, object],
    contaminant_policy: str,
    reverse_policy: str,
    warnings: tuple[str, ...],
) -> ImporterQualityReport:
    format_specific = dict(mapped_report.format_specific)
    format_specific["maxquant"] = {
        "resolved_columns": _resolved_columns_payload(resolved),
        "filtering": dict(filter_diagnostics),
        "adaptation": dict(adapter_diagnostics),
    }
    return replace(
        mapped_report,
        rows_read=_diagnostic_int(filter_diagnostics, "input_row_count"),
        rows_retained=_diagnostic_int(filter_diagnostics, "retained_row_count"),
        rows_dropped=_diagnostic_int(filter_diagnostics, "removed_rows"),
        localisation_confidence=replace(
            mapped_report.localisation_confidence,
            source_column=resolved.localisation_confidence,
        ),
        flagged_rows=ImporterFlaggedRowSummary(
            contaminant=_maxquant_flag_quality_count(
                count=_diagnostic_int(
                    filter_diagnostics,
                    "potential_contaminant_rows",
                ),
                source_column=resolved.potential_contaminant,
                policy=contaminant_policy,
                missing_reason=("MaxQuant potential contaminant column was not found"),
            ),
            reverse=_maxquant_flag_quality_count(
                count=_diagnostic_int(filter_diagnostics, "reverse_rows"),
                source_column=resolved.reverse,
                policy=reverse_policy,
                missing_reason="MaxQuant reverse column was not found",
            ),
            decoy=ImporterQualityCount(
                status=IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
                reason="MaxQuant importer does not parse a separate decoy flag",
            ),
        ),
        format_specific=format_specific,
        warnings=warnings,
    )


def _maxquant_flag_quality_count(
    *,
    count: int,
    source_column: str | None,
    policy: str,
    missing_reason: str,
) -> ImporterQualityCount:
    if source_column is None:
        return ImporterQualityCount(
            status=IMPORTER_QUALITY_STATUS_NOT_APPLICABLE,
            reason=missing_reason,
        )
    return ImporterQualityCount(
        status=IMPORTER_QUALITY_STATUS_REPORTED,
        count=count,
        source_column=source_column,
        policy=policy,
    )


def _diagnostic_int(payload: dict[str, object], field_name: str) -> int:
    value = payload[field_name]
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(f"MaxQuant diagnostic {field_name} must be an int")
    return int(value)


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
        importer_label="MaxQuant",
        validate_column_name=validate_optional_maxquant_column_name,
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
        importer_label="MaxQuant",
        required=required,
        validate_column_name=validate_optional_maxquant_column_name,
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
        importer_label="MaxQuant",
        request_label="maxquant",
        mapping_class_name="MaxQuantColumnMapping",
        reject_duplicate_inferred_sample_ids=True,
    )


def _raise_for_forbidden_flags(
    values: pd.Series,
    *,
    policy: str,
    label: str,
) -> None:
    raise_for_forbidden_flags(
        values,
        policy=policy,
        error_policy=MAXQUANT_FLAG_POLICY_ERROR,
        importer_label="MaxQuant",
        label=label,
    )


def _resolve_row_site_tokens(
    row: Mapping[str, object],
    *,
    resolved: _ResolvedMaxQuantColumns,
    row_position: int,
) -> tuple[str, ...]:
    if resolved.modified_site is not None:
        return _parse_site_tokens(
            row[resolved.modified_site],
            field_name=f"MaxQuant {resolved.modified_site}",
            row_position=row_position,
        )
    if resolved.amino_acid is None or resolved.site_position is None:
        raise PhosPyInputError(
            "MaxQuant importer requires modified_site or amino_acid/site_position"
        )
    residue_tokens = split_multi_value(row[resolved.amino_acid])
    position_tokens = split_multi_value(row[resolved.site_position])
    protein_tokens = split_multi_value(row[resolved.protein_accession])
    if (
        len(residue_tokens) == 1
        and len(position_tokens) > 1
        and len(protein_tokens) > 1
    ):
        position_tokens = position_tokens[:1]
    if len(residue_tokens) == 1 and len(position_tokens) > 1:
        residue_tokens = residue_tokens * len(position_tokens)
    if len(residue_tokens) != len(position_tokens):
        raise PhosPyInputError(
            "MaxQuant amino_acid and site_position columns must describe the same "
            f"number of site tokens; row_position={row_position}, "
            f"amino_acid={row[resolved.amino_acid]!r}, "
            f"site_position={row[resolved.site_position]!r}"
        )
    raw_site = ";".join(
        f"{_normalise_residue(residue, row_position=row_position)}"
        f"{_normalise_position(position, row_position=row_position)}"
        for residue, position in zip(residue_tokens, position_tokens, strict=True)
    )
    return _parse_site_tokens(
        raw_site,
        field_name="MaxQuant amino_acid/site_position",
        row_position=row_position,
    )


def _parse_site_tokens(
    value: object,
    *,
    field_name: str,
    row_position: int,
) -> tuple[str, ...]:
    return tuple(
        token.token
        for token in parse_phospho_site_tokens(
            value,
            field_name=f"{field_name} row_position={row_position}",
        )
    )


def _normalise_residue(value: object, *, row_position: int) -> str:
    token = required_text(
        value,
        field_name="MaxQuant amino_acid",
        row_position=row_position,
    ).upper()
    if token not in {"S", "T", "Y"}:
        raise PhosPyInputError(
            "MaxQuant amino_acid values must be one of 'S', 'T', or 'Y'; "
            f"row_position={row_position}, value={value!r}"
        )
    return token


def _normalise_position(value: object, *, row_position: int) -> int:
    token = required_text(
        value,
        field_name="MaxQuant site_position",
        row_position=row_position,
    )
    try:
        position_value = float(token)
    except ValueError as exc:
        raise PhosPyInputError(
            "MaxQuant site_position values must be positive integers; "
            f"row_position={row_position}, value={value!r}"
        ) from exc
    if not math.isfinite(position_value) or not position_value.is_integer():
        raise PhosPyInputError(
            "MaxQuant site_position values must be positive integers; "
            f"row_position={row_position}, value={value!r}"
        )
    position = int(position_value)
    if position < 1:
        raise PhosPyInputError(
            "MaxQuant site_position values must be positive integers; "
            f"row_position={row_position}, value={value!r}"
        )
    return int(position)


def _normalise_localisation_source_values(values: pd.Series) -> list[object]:
    normalised: list[object] = []
    for value in values.tolist():
        parsed_probability_string = _parse_maxquant_probability_tokens(value)
        normalised.append(
            value if parsed_probability_string is None else parsed_probability_string
        )
    return normalised


def _parse_maxquant_probability_tokens(value: object) -> float | None:
    if not isinstance(value, str):
        return None
    matches = _LOCALISATION_PROBABILITY_TOKEN_PATTERN.findall(value)
    if not matches:
        return None
    probabilities: list[float] = []
    for match in matches:
        try:
            numeric = float(match)
        except ValueError:
            return None
        if not math.isfinite(numeric) or numeric < 0.0 or numeric > 1.0:
            return None
        probabilities.append(numeric)
    if not probabilities:
        return None
    return float(min(probabilities))


def _resolve_modified_peptide_sequences(
    source: pd.DataFrame,
    *,
    resolved: _ResolvedMaxQuantColumns,
) -> list[str]:
    if resolved.modified_peptide_sequence is None:
        return [
            required_text(
                value,
                field_name=f"MaxQuant {resolved.peptide_sequence}",
                row_position=position,
            )
            for position, value in enumerate(source.loc[:, resolved.peptide_sequence])
        ]
    return [
        required_text(
            value,
            field_name=f"MaxQuant {resolved.modified_peptide_sequence}",
            row_position=position,
        )
        for position, value in enumerate(
            source.loc[:, resolved.modified_peptide_sequence]
        )
    ]


def _build_row_ids(
    *,
    source: pd.DataFrame,
    resolved: _ResolvedMaxQuantColumns,
    protein_values: list[str],
    site_values: list[str],
    source_row_numbers: list[int],
) -> list[str]:
    return build_row_ids(
        source=source,
        explicit_column=resolved.row_id,
        protein_values=protein_values,
        site_values=site_values,
        source_row_numbers=source_row_numbers,
        importer_label="MaxQuant",
        generated_prefix="maxquant",
    )


def _build_unique_feature_ids(
    *,
    source: pd.DataFrame,
    resolved: _ResolvedMaxQuantColumns,
    source_row_numbers: list[int],
) -> list[str]:
    return build_unique_feature_ids(
        source=source,
        explicit_column=resolved.unique_feature_id,
        source_row_numbers=source_row_numbers,
        importer_label="MaxQuant",
        generated_prefix="maxquant",
    )


def _resolved_columns_payload(resolved: _ResolvedMaxQuantColumns) -> dict[str, object]:
    return {
        "protein_accession": resolved.protein_accession,
        "gene_symbol": resolved.gene_symbol,
        "modified_site": resolved.modified_site,
        "amino_acid": resolved.amino_acid,
        "site_position": resolved.site_position,
        "localisation_confidence": resolved.localisation_confidence,
        "peptide_sequence": resolved.peptide_sequence,
        "modified_peptide_sequence": resolved.modified_peptide_sequence,
        "intensity_columns": dict(resolved.intensity_columns),
        "potential_contaminant": resolved.potential_contaminant,
        "reverse": resolved.reverse,
        "row_id": resolved.row_id,
        "unique_feature_id": resolved.unique_feature_id,
        "site_sequence": resolved.site_sequence,
    }


def _require_non_empty_unique_columns(source: pd.DataFrame) -> None:
    require_non_empty_unique_columns(source, importer_label="MaxQuant")


__all__ = [
    "MaxQuantColumnMapping",
    "MaxQuantPhosphositeImportRequest",
    "MaxQuantPhosphositeImporter",
]
