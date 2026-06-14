"""FragPipe/Philosopher/PTMProphet phosphosite importer."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import pandas as pd

from phospy.contracts.requests import PhosphositeImportRequest
from phospy.contracts.results import PhosphositeImportResult
from phospy.errors.input import PhosPyInputError
from phospy.io.readers.importers import (
    MappedPhosphositeTableImporter,
    _read_upstream_table,
)
from phospy.science.evidence.modified_peptides import (
    ModifiedResidue,
    parse_modified_peptide_sequence,
)
from phospy.science.evidence.multi_site import parse_phospho_site_tokens
from phospy.validation.datasets.fragpipe import (
    FRAGPIPE_FLAG_POLICY_ERROR,
    FRAGPIPE_FLAG_POLICY_FLAG,
    FRAGPIPE_FLAG_POLICY_REMOVE,
    FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PEPTIDE,
    FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PROTEIN,
    validate_fragpipe_flag_policy,
    validate_optional_fragpipe_column_name,
    validate_ptmprophet_position_reference,
)
from phospy.validation.datasets.importers import normalise_sample_column_mapping

_ADAPTED_ROW_ID_COLUMN = "__phospy_fragpipe_row_id"
_ADAPTED_PROTEIN_ACCESSION_COLUMN = "__phospy_fragpipe_protein_accession"
_ADAPTED_PROTEIN_ID_COLUMN = "__phospy_fragpipe_protein_id"
_ADAPTED_GENE_SYMBOL_COLUMN = "__phospy_fragpipe_gene_symbol"
_ADAPTED_SITE_COLUMN = "__phospy_fragpipe_site"
_ADAPTED_SITE_SEQUENCE_COLUMN = "__phospy_fragpipe_site_sequence"
_ADAPTED_LOCALISATION_COLUMN = "__phospy_fragpipe_localisation_confidence"
_ADAPTED_UNIQUE_FEATURE_ID_COLUMN = "__phospy_fragpipe_feature_id"
_ADAPTED_PEPTIDE_SEQUENCE_COLUMN = "__phospy_fragpipe_peptide_sequence"
_ADAPTED_MODIFIED_PEPTIDE_SEQUENCE_COLUMN = (
    "__phospy_fragpipe_modified_peptide_sequence"
)
_ADAPTED_PEPTIDE_SITE_STRING_COLUMN = "__phospy_fragpipe_peptide_site_string"
_ADAPTED_CANDIDATE_SITES_COLUMN = "fragpipe_ptmprophet_candidate_sites"
_ADAPTED_SITE_PROBABILITIES_COLUMN = "fragpipe_ptmprophet_site_probabilities"
_ADAPTED_AMBIGUOUS_COLUMN = "fragpipe_ptmprophet_ambiguous"
_ADAPTED_MODIFIED_PHOSPHO_COUNT_COLUMN = "fragpipe_modified_peptide_phospho_count"
_FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN = "fragpipe_contaminant"
_FRAGPIPE_DECOY_OUTPUT_COLUMN = "fragpipe_decoy"
_DEFAULT_INTENSITY_PREFIXES = (
    "Intensity ",
    "LFQ Intensity ",
    "MaxLFQ Intensity ",
    "Abundance ",
    "Area ",
)
_PROTEIN_ACCESSION_CANDIDATES = (
    "Protein",
    "Protein ID",
    "Protein IDs",
    "Protein ID(s)",
    "Proteins",
    "Mapped Proteins",
    "Protein Accession",
    "Protein Accession(s)",
)
_GENE_SYMBOL_CANDIDATES = (
    "Gene",
    "Gene Name",
    "Gene names",
    "Genes",
    "Mapped Genes",
)
_PEPTIDE_SEQUENCE_CANDIDATES = (
    "Peptide",
    "Peptide Sequence",
    "Sequence",
)
_MODIFIED_PEPTIDE_SEQUENCE_CANDIDATES = (
    "Modified Peptide",
    "Modified Peptide Sequence",
    "Modified Sequence",
    "Modified Peptide Sequence With Flanking AAs",
)
_PTMPROPHET_PROBABILITY_CANDIDATES = (
    "PTMProphet Probability",
    "PTMProphet Probabilities",
    "PTMProphet Site Probabilities",
    "PTMProphet Localization",
    "PTMProphet Localisation",
    "Localization Probability",
    "Localisation Probability",
    "Best Localization",
    "Best Localisation",
)
_PROTEIN_START_CANDIDATES = (
    "Protein Start",
    "Protein Start Position",
    "Start",
    "Start Position",
    "Peptide Start",
    "Mapped Start",
)
_SITE_CANDIDATES = (
    "Site",
    "Phosphosite",
    "Phospho Site",
    "Modified Site",
)
_SITE_SEQUENCE_CANDIDATES = (
    "Sequence Window",
    "Sequence window",
    "Window Sequence",
    "Site Sequence",
)
_UNIQUE_FEATURE_CANDIDATES = (
    "Spectrum",
    "Spectrum ID",
    "PSM ID",
    "Peptide ID",
    "Index",
    "id",
    "ID",
)
_CONTAMINANT_CANDIDATES = (
    "Contaminant",
    "Potential contaminant",
    "Potential Contaminant",
    "Is Contaminant",
)
_DECOY_CANDIDATES = (
    "Decoy",
    "Reverse",
    "Is Decoy",
    "Protein Decoy",
)
_ROW_ID_CANDIDATES: tuple[str, ...] = ()
_MULTI_VALUE_SPLIT_PATTERN = re.compile(r"\s*[,;]\s*")
_NUMERIC_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
_POSITIONED_LOCALISATION_PATTERNS = (
    re.compile(
        rf"(?<![A-Za-z])([STYsty])\s*([1-9][0-9]*)\s*\(\s*({_NUMERIC_PATTERN})\s*\)"
    ),
    re.compile(
        rf"(?<![A-Za-z])([STYsty])\s*([1-9][0-9]*)\s*[:=]\s*({_NUMERIC_PATTERN})"
    ),
    re.compile(
        rf"(?<![A-Za-z])([STYsty])\s*@\s*([1-9][0-9]*)\s*[:=]\s*({_NUMERIC_PATTERN})"
    ),
    re.compile(
        rf"(?<![A-Za-z])([STYsty])\s*\(\s*([1-9][0-9]*)\s*[,;:]\s*({_NUMERIC_PATTERN})\s*\)"
    ),
)
_RESIDUE_ONLY_LOCALISATION_PATTERN = re.compile(
    rf"(?<![A-Za-z0-9])([STYsty])\s*\(\s*({_NUMERIC_PATTERN})\s*\)"
)
_CONTAMINANT_PREFIXES = ("CON__", "CON_", "CONTAM_", "CONTAMINANT_")
_DECOY_PREFIXES = ("REV__", "REV_", "DECOY__", "DECOY_")


@dataclass(frozen=True, slots=True)
class FragPipeColumnMapping:
    """Optional source-column overrides for FragPipe/PTMProphet imports."""

    protein_accession: str | None = None
    gene_symbol: str | None = None
    peptide_sequence: str | None = None
    modified_peptide_sequence: str | None = None
    ptmprophet_probabilities: str | None = None
    protein_start: str | None = None
    site: str | None = None
    site_sequence: str | None = None
    intensity_columns: Mapping[str, str] | Sequence[str] | None = None
    contaminant: str | None = None
    decoy: str | None = None
    row_id: str | None = None
    unique_feature_id: str | None = None


@dataclass(frozen=True, slots=True)
class FragPipePTMProphetImportRequest:
    """Request for importing FragPipe/Philosopher/PTMProphet phosphosite output."""

    source: object
    column_mapping: FragPipeColumnMapping = field(default_factory=FragPipeColumnMapping)
    contaminant_policy: str = FRAGPIPE_FLAG_POLICY_REMOVE
    decoy_policy: str = FRAGPIPE_FLAG_POLICY_REMOVE
    ptmprophet_position_reference: str = FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PEPTIDE
    intensity_column_prefixes: Sequence[str] = _DEFAULT_INTENSITY_PREFIXES
    source_name: str = "fragpipe_ptmprophet"


@dataclass(frozen=True, slots=True)
class _ResolvedFragPipeColumns:
    protein_accession: str
    gene_symbol: str
    peptide_sequence: str
    modified_peptide_sequence: str
    ptmprophet_probabilities: str
    protein_start: str | None
    site: str | None
    site_sequence: str | None
    intensity_columns: dict[str, str]
    contaminant: str | None
    decoy: str | None
    row_id: str | None
    unique_feature_id: str | None


@dataclass(frozen=True, slots=True)
class _LocalisationCandidate:
    residue: str
    position: int | None
    probability: float | None

    @property
    def has_position(self) -> bool:
        return self.position is not None


@dataclass(frozen=True, slots=True)
class _ProteinSiteCandidate:
    residue: str
    protein_position: int
    probability: float | None

    @property
    def token(self) -> str:
        return f"{self.residue}{self.protein_position}"


@dataclass(frozen=True, slots=True)
class _SiteCall:
    site_tokens: tuple[str, ...]
    peptide_site_string: str
    localisation_confidence: object
    candidate_sites: str
    site_probabilities: str
    ambiguous: bool
    phospho_site_count: int


class FragPipePTMProphetImporter:
    """Import FragPipe/PTMProphet phosphosite output into PhosPy candidates."""

    def __init__(
        self,
        *,
        mapped_importer: MappedPhosphositeTableImporter | None = None,
    ) -> None:
        self._mapped_importer = mapped_importer or MappedPhosphositeTableImporter()

    def run(
        self,
        request: FragPipePTMProphetImportRequest,
    ) -> PhosphositeImportResult:
        if not isinstance(request, FragPipePTMProphetImportRequest):
            raise PhosPyInputError(
                "FragPipe importer input must be a FragPipePTMProphetImportRequest"
            )
        contaminant_policy = validate_fragpipe_flag_policy(
            request.contaminant_policy,
            field_name="fragpipe import request contaminant_policy",
        )
        decoy_policy = validate_fragpipe_flag_policy(
            request.decoy_policy,
            field_name="fragpipe import request decoy_policy",
        )
        position_reference = validate_ptmprophet_position_reference(
            request.ptmprophet_position_reference,
            field_name="fragpipe import request ptmprophet_position_reference",
        )
        source = _read_upstream_table(request.source)
        _require_non_empty_unique_columns(source)
        resolved = _resolve_fragpipe_columns(
            source,
            request.column_mapping,
            intensity_column_prefixes=request.intensity_column_prefixes,
        )
        filtered, flags, filter_diagnostics, filter_warnings = _apply_flag_policies(
            source,
            resolved=resolved,
            contaminant_policy=contaminant_policy,
            decoy_policy=decoy_policy,
        )
        if filtered.empty:
            raise PhosPyInputError(
                "FragPipe importer removed all rows after contaminant/decoy filtering"
            )

        adapted, adapter_diagnostics, adapter_warnings = _adapt_fragpipe_source(
            filtered,
            resolved=resolved,
            ptmprophet_position_reference=position_reference,
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
                localisation_confidence_scale="probability",
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
            decoy_policy=decoy_policy,
            resolved=resolved,
            filter_diagnostics=filter_diagnostics,
            adapter_diagnostics=adapter_diagnostics,
            warnings=filter_warnings + adapter_warnings,
        )


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


def _apply_flag_policies(
    source: pd.DataFrame,
    *,
    resolved: _ResolvedFragPipeColumns,
    contaminant_policy: str,
    decoy_policy: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object], tuple[str, ...]]:
    protein_values = source.loc[:, resolved.protein_accession]
    contaminant_flags = _resolve_flag_series(
        source,
        column=resolved.contaminant,
        field_name="FragPipe contaminant flag",
    )
    decoy_flags = _resolve_flag_series(
        source,
        column=resolved.decoy,
        field_name="FragPipe decoy flag",
    )
    contaminant_prefix_flags = protein_values.map(
        lambda value: _has_any_prefixed_token(value, prefixes=_CONTAMINANT_PREFIXES)
    )
    decoy_prefix_flags = protein_values.map(
        lambda value: _has_any_prefixed_token(value, prefixes=_DECOY_PREFIXES)
    )

    flags = pd.DataFrame(index=source.index.copy())
    flags[_FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN] = (
        contaminant_prefix_flags.astype(bool)
        if contaminant_flags is None
        else (contaminant_flags.astype(bool) | contaminant_prefix_flags.astype(bool))
    ).tolist()
    flags[_FRAGPIPE_DECOY_OUTPUT_COLUMN] = (
        decoy_prefix_flags.astype(bool)
        if decoy_flags is None
        else (decoy_flags.astype(bool) | decoy_prefix_flags.astype(bool))
    ).tolist()

    _raise_for_forbidden_flags(
        flags[_FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN],
        policy=contaminant_policy,
        label="contaminant",
    )
    _raise_for_forbidden_flags(
        flags[_FRAGPIPE_DECOY_OUTPUT_COLUMN],
        policy=decoy_policy,
        label="decoy",
    )
    keep_mask = pd.Series(True, index=source.index.copy(), dtype=bool)
    if contaminant_policy == FRAGPIPE_FLAG_POLICY_REMOVE:
        keep_mask &= ~flags[_FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN].astype(bool)
    if decoy_policy == FRAGPIPE_FLAG_POLICY_REMOVE:
        keep_mask &= ~flags[_FRAGPIPE_DECOY_OUTPUT_COLUMN].astype(bool)

    diagnostics = {
        "input_row_count": int(source.shape[0]),
        "contaminant_column": resolved.contaminant,
        "decoy_column": resolved.decoy,
        "contaminant_policy": contaminant_policy,
        "decoy_policy": decoy_policy,
        "contaminant_rows": int(
            flags[_FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN].astype(bool).sum()
        ),
        "decoy_rows": int(flags[_FRAGPIPE_DECOY_OUTPUT_COLUMN].astype(bool).sum()),
        "removed_rows": int((~keep_mask).sum()),
        "retained_row_count": int(keep_mask.sum()),
        "contaminant_prefix_rows": int(contaminant_prefix_flags.astype(bool).sum()),
        "decoy_prefix_rows": int(decoy_prefix_flags.astype(bool).sum()),
    }
    warnings: list[str] = []
    if resolved.contaminant is None:
        warnings.append(
            "FragPipe contaminant column was not found; contaminant filtering used "
            "protein accession prefixes only"
        )
    if resolved.decoy is None:
        warnings.append(
            "FragPipe decoy column was not found; decoy filtering used protein "
            "accession prefixes only"
        )
    filtered = source.loc[keep_mask, :].copy(deep=True)
    filtered_flags = flags.loc[keep_mask, :].copy(deep=True)
    return filtered, filtered_flags, diagnostics, tuple(warnings)


def _adapt_fragpipe_source(
    source: pd.DataFrame,
    *,
    resolved: _ResolvedFragPipeColumns,
    ptmprophet_position_reference: str,
) -> tuple[pd.DataFrame, dict[str, object], tuple[str, ...]]:
    adapted = source.copy(deep=True)
    source_row_numbers = [position + 1 for position in range(int(source.shape[0]))]
    protein_values: list[str] = []
    gene_values: list[str] = []
    site_values: list[str] = []
    peptide_site_values: list[str] = []
    localisation_values: list[object] = []
    candidate_site_values: list[str] = []
    site_probability_values: list[str] = []
    ambiguous_values: list[bool] = []
    phospho_counts: list[int] = []
    protein_group_rows = 0
    peptide_sequence_mismatch_rows = 0

    columns = source.columns.astype(str).tolist()
    for position, row in enumerate(source.itertuples(index=False, name=None)):
        row_lookup = dict(zip(columns, row, strict=True))
        if _multi_value_count(row_lookup[resolved.protein_accession]) > 1:
            protein_group_rows += 1
        protein_values.append(
            _parse_protein_accession(
                row_lookup[resolved.protein_accession],
                field_name=f"FragPipe {resolved.protein_accession}",
                row_position=position,
            )
        )
        gene_values.append(
            _first_list_token(
                row_lookup[resolved.gene_symbol],
                field_name=f"FragPipe {resolved.gene_symbol}",
                row_position=position,
            )
        )
        parsed_modified = parse_modified_peptide_sequence(
            row_lookup[resolved.modified_peptide_sequence],
            field_name=(
                f"FragPipe {resolved.modified_peptide_sequence} row_position={position}"
            ),
        )
        peptide_sequence = _required_text(
            row_lookup[resolved.peptide_sequence],
            field_name=f"FragPipe {resolved.peptide_sequence}",
            row_position=position,
        )
        if peptide_sequence.strip().upper() != parsed_modified.sequence:
            peptide_sequence_mismatch_rows += 1
        site_call = _resolve_site_call(
            row_lookup,
            resolved=resolved,
            modified_phospho_sites=parsed_modified.phospho_sites,
            ptmprophet_position_reference=ptmprophet_position_reference,
            row_position=position,
        )
        site_values.append(",".join(site_call.site_tokens))
        peptide_site_values.append(site_call.peptide_site_string)
        localisation_values.append(site_call.localisation_confidence)
        candidate_site_values.append(site_call.candidate_sites)
        site_probability_values.append(site_call.site_probabilities)
        ambiguous_values.append(site_call.ambiguous)
        phospho_counts.append(site_call.phospho_site_count)

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
    adapted[_ADAPTED_LOCALISATION_COLUMN] = localisation_values
    adapted[_ADAPTED_UNIQUE_FEATURE_ID_COLUMN] = _build_unique_feature_ids(
        source=source,
        resolved=resolved,
        source_row_numbers=source_row_numbers,
    )
    adapted[_ADAPTED_PEPTIDE_SEQUENCE_COLUMN] = [
        _required_text(
            value,
            field_name=f"FragPipe {resolved.peptide_sequence}",
            row_position=position,
        )
        for position, value in enumerate(source.loc[:, resolved.peptide_sequence])
    ]
    adapted[_ADAPTED_MODIFIED_PEPTIDE_SEQUENCE_COLUMN] = [
        _required_text(
            value,
            field_name=f"FragPipe {resolved.modified_peptide_sequence}",
            row_position=position,
        )
        for position, value in enumerate(
            source.loc[:, resolved.modified_peptide_sequence]
        )
    ]
    adapted[_ADAPTED_PEPTIDE_SITE_STRING_COLUMN] = peptide_site_values
    adapted[_ADAPTED_CANDIDATE_SITES_COLUMN] = candidate_site_values
    adapted[_ADAPTED_SITE_PROBABILITIES_COLUMN] = site_probability_values
    adapted[_ADAPTED_AMBIGUOUS_COLUMN] = ambiguous_values
    adapted[_ADAPTED_MODIFIED_PHOSPHO_COUNT_COLUMN] = phospho_counts
    if resolved.site_sequence is not None:
        adapted[_ADAPTED_SITE_SEQUENCE_COLUMN] = [
            _optional_text(value)
            for value in source.loc[:, resolved.site_sequence].tolist()
        ]

    diagnostics = {
        "resolved_columns": _resolved_columns_payload(resolved),
        "ptmprophet_position_reference": ptmprophet_position_reference,
        "protein_group_rows_collapsed_to_first_accession": int(protein_group_rows),
        "peptide_sequence_mismatch_rows": int(peptide_sequence_mismatch_rows),
        "ambiguous_localisation_rows": int(sum(ambiguous_values)),
        "multi_site_rows": int(
            sum(1 for site_value in site_values if len(site_value.split(",")) > 1)
        ),
    }
    warnings: list[str] = []
    if protein_group_rows:
        warnings.append(
            "FragPipe protein-group rows were represented by the first listed "
            "protein accession for protein-scoped identity"
        )
    if peptide_sequence_mismatch_rows:
        warnings.append(
            "FragPipe peptide sequence and parsed modified peptide sequence differed "
            "for some rows; modified peptide parsing was retained for site evidence"
        )
    if any(ambiguous_values):
        warnings.append(
            "FragPipe/PTMProphet ambiguous localisation rows were retained as "
            "joint multi-site observations rather than selecting the first site"
        )
    return adapted, diagnostics, tuple(warnings)


def _resolve_site_call(
    row: Mapping[str, object],
    *,
    resolved: _ResolvedFragPipeColumns,
    modified_phospho_sites: tuple[ModifiedResidue, ...],
    ptmprophet_position_reference: str,
    row_position: int,
) -> _SiteCall:
    localisation_candidates = _parse_ptmprophet_localisation_candidates(
        row[resolved.ptmprophet_probabilities],
        modified_phospho_sites=modified_phospho_sites,
        field_name=f"FragPipe {resolved.ptmprophet_probabilities}",
        row_position=row_position,
    )
    protein_start: int | None = None
    protein_start_was_resolved = False

    def resolve_protein_start() -> int | None:
        nonlocal protein_start, protein_start_was_resolved
        if not protein_start_was_resolved:
            protein_start = _resolve_protein_start(
                row,
                resolved=resolved,
                row_position=row_position,
            )
            protein_start_was_resolved = True
        return protein_start

    if resolved.site is not None:
        protein_start_for_localisation = (
            resolve_protein_start()
            if localisation_candidates
            and ptmprophet_position_reference
            == FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PEPTIDE
            else None
        )
        selected = tuple(
            _ProteinSiteCandidate(
                residue=token.residue,
                protein_position=token.position,
                probability=_probability_for_site_token(
                    token.token,
                    candidates=localisation_candidates,
                    protein_start=protein_start_for_localisation,
                    ptmprophet_position_reference=ptmprophet_position_reference,
                    row_position=row_position,
                ),
            )
            for token in parse_phospho_site_tokens(
                row[resolved.site],
                field_name=f"FragPipe {resolved.site} row_position={row_position}",
            )
        )
        all_candidates = _convert_candidates_to_protein_sites(
            localisation_candidates,
            protein_start=protein_start_for_localisation,
            ptmprophet_position_reference=ptmprophet_position_reference,
            row_position=row_position,
        )
        if not all_candidates:
            all_candidates = selected
        _, ambiguous = _select_localised_sites(
            all_candidates,
            phospho_site_count=len(modified_phospho_sites),
        )
        return _site_call_from_candidates(
            selected,
            all_candidates=all_candidates,
            ambiguous=ambiguous,
            phospho_site_count=len(modified_phospho_sites),
        )

    if localisation_candidates:
        protein_start_for_localisation = (
            resolve_protein_start()
            if ptmprophet_position_reference
            == FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PEPTIDE
            else None
        )
        protein_candidates = _convert_candidates_to_protein_sites(
            localisation_candidates,
            protein_start=protein_start_for_localisation,
            ptmprophet_position_reference=ptmprophet_position_reference,
            row_position=row_position,
        )
        selected, ambiguous = _select_localised_sites(
            protein_candidates,
            phospho_site_count=len(modified_phospho_sites),
        )
        return _site_call_from_candidates(
            selected,
            all_candidates=protein_candidates,
            ambiguous=ambiguous,
            phospho_site_count=len(modified_phospho_sites),
        )

    if not modified_phospho_sites:
        raise PhosPyInputError(
            "FragPipe importer could not extract phosphosite evidence from "
            f"modified peptide or PTMProphet localisation; row_position={row_position}"
        )
    fallback_candidates = tuple(
        _protein_site_from_peptide_position(
            residue=site.residue,
            peptide_position=site.position,
            probability=None,
            protein_start=resolve_protein_start(),
            row_position=row_position,
        )
        for site in modified_phospho_sites
    )
    return _site_call_from_candidates(
        fallback_candidates,
        all_candidates=fallback_candidates,
        ambiguous=False,
        phospho_site_count=len(modified_phospho_sites),
    )


def _parse_ptmprophet_localisation_candidates(
    value: object,
    *,
    modified_phospho_sites: tuple[ModifiedResidue, ...],
    field_name: str,
    row_position: int,
) -> tuple[_LocalisationCandidate, ...]:
    if _is_missing(value):
        return ()
    if isinstance(value, int | float) and not isinstance(value, bool):
        probability = _parse_probability(
            value, field_name=field_name, row_position=row_position
        )
        return tuple(
            _LocalisationCandidate(
                residue=site.residue,
                position=site.position,
                probability=probability,
            )
            for site in modified_phospho_sites
        )
    if not isinstance(value, str):
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} must be a PTMProphet "
            "localisation string or probability"
        )
    token = value.strip()
    if token == "":
        return ()

    positioned: list[_LocalisationCandidate] = []
    consumed_spans: set[tuple[int, int]] = set()
    for pattern in _POSITIONED_LOCALISATION_PATTERNS:
        for match in pattern.finditer(token):
            span = match.span()
            if span in consumed_spans:
                continue
            consumed_spans.add(span)
            positioned.append(
                _LocalisationCandidate(
                    residue=match.group(1).upper(),
                    position=int(match.group(2)),
                    probability=_parse_probability(
                        match.group(3),
                        field_name=field_name,
                        row_position=row_position,
                    ),
                )
            )
    if positioned:
        return tuple(_deduplicate_localisation_candidates(positioned))

    residue_only_matches = list(_RESIDUE_ONLY_LOCALISATION_PATTERN.finditer(token))
    if residue_only_matches:
        if len(residue_only_matches) > len(modified_phospho_sites):
            raise PhosPyInputError(
                f"{field_name} row_position={row_position} has "
                "residue-only probabilities but more probability tokens than parsed "
                "modified phosphosites"
            )
        candidates: list[_LocalisationCandidate] = []
        search_start = 0
        for match in residue_only_matches:
            residue = match.group(1).upper()
            mapped_site = _next_modified_site_with_residue(
                modified_phospho_sites,
                residue=residue,
                start=search_start,
                row_position=row_position,
                field_name=field_name,
            )
            candidates.append(
                _LocalisationCandidate(
                    residue=residue,
                    position=mapped_site.position,
                    probability=_parse_probability(
                        match.group(2),
                        field_name=field_name,
                        row_position=row_position,
                    ),
                )
            )
            search_start = modified_phospho_sites.index(mapped_site) + 1
        return tuple(candidates)

    try:
        probability = _parse_probability(
            token,
            field_name=field_name,
            row_position=row_position,
        )
    except PhosPyInputError as exc:
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} contains malformed "
            f"PTMProphet localisation string {value!r}; expected tokens like "
            "'S3(0.95)' or 'S3:0.95'"
        ) from exc
    return tuple(
        _LocalisationCandidate(
            residue=site.residue,
            position=site.position,
            probability=probability,
        )
        for site in modified_phospho_sites
    )


def _convert_candidates_to_protein_sites(
    candidates: tuple[_LocalisationCandidate, ...],
    *,
    protein_start: int | None,
    ptmprophet_position_reference: str,
    row_position: int,
) -> tuple[_ProteinSiteCandidate, ...]:
    protein_candidates: list[_ProteinSiteCandidate] = []
    for candidate in candidates:
        if candidate.position is None:
            raise PhosPyInputError(
                "FragPipe PTMProphet localisation candidate is missing a position; "
                f"row_position={row_position}"
            )
        if (
            ptmprophet_position_reference
            == FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PROTEIN
        ):
            protein_candidates.append(
                _ProteinSiteCandidate(
                    residue=candidate.residue,
                    protein_position=candidate.position,
                    probability=candidate.probability,
                )
            )
            continue
        protein_candidates.append(
            _protein_site_from_peptide_position(
                residue=candidate.residue,
                peptide_position=candidate.position,
                probability=candidate.probability,
                protein_start=protein_start,
                row_position=row_position,
            )
        )
    return tuple(_deduplicate_protein_candidates(protein_candidates))


def _select_localised_sites(
    candidates: tuple[_ProteinSiteCandidate, ...],
    *,
    phospho_site_count: int,
) -> tuple[tuple[_ProteinSiteCandidate, ...], bool]:
    if not candidates:
        return (), False
    if phospho_site_count <= 1 and len(candidates) > 1:
        probabilities = [candidate.probability for candidate in candidates]
        if any(probability is None for probability in probabilities):
            return candidates, True
        numeric_probabilities = tuple(
            probability for probability in probabilities if probability is not None
        )
        max_probability = max(numeric_probabilities)
        top_candidates = tuple(
            candidate
            for candidate in candidates
            if candidate.probability is not None
            and math.isclose(
                float(candidate.probability), max_probability, abs_tol=1e-12
            )
        )
        return top_candidates, len(top_candidates) > 1
    if phospho_site_count > 1 and len(candidates) > phospho_site_count:
        return candidates, True
    return candidates, False


def _site_call_from_candidates(
    selected: tuple[_ProteinSiteCandidate, ...],
    *,
    all_candidates: tuple[_ProteinSiteCandidate, ...],
    ambiguous: bool,
    phospho_site_count: int,
) -> _SiteCall:
    if not selected:
        raise PhosPyInputError("FragPipe importer resolved zero phosphosite candidates")
    ordered_selected = tuple(_deduplicate_protein_candidates(selected))
    site_tokens = tuple(candidate.token for candidate in ordered_selected)
    selected_probabilities = [candidate.probability for candidate in ordered_selected]
    if selected_probabilities and all(
        probability is not None for probability in selected_probabilities
    ):
        numeric_selected_probabilities = tuple(
            probability
            for probability in selected_probabilities
            if probability is not None
        )
        localisation_confidence: object = float(min(numeric_selected_probabilities))
    else:
        localisation_confidence = pd.NA
    return _SiteCall(
        site_tokens=site_tokens,
        peptide_site_string=";".join(site_tokens),
        localisation_confidence=localisation_confidence,
        candidate_sites=";".join(
            candidate.token
            for candidate in _deduplicate_protein_candidates(all_candidates)
        ),
        site_probabilities=_format_site_probabilities(all_candidates),
        ambiguous=bool(ambiguous),
        phospho_site_count=int(phospho_site_count),
    )


def _protein_site_from_peptide_position(
    *,
    residue: str,
    peptide_position: int,
    probability: float | None,
    protein_start: int | None,
    row_position: int,
) -> _ProteinSiteCandidate:
    if protein_start is None:
        raise PhosPyInputError(
            "FragPipe importer requires a protein_start column when PTMProphet "
            "positions are peptide-relative and no explicit site column is mapped; "
            f"row_position={row_position}"
        )
    protein_position = int(protein_start) + int(peptide_position) - 1
    if protein_position < 1:
        raise PhosPyInputError(
            f"FragPipe computed invalid protein position {protein_position}; "
            f"row_position={row_position}"
        )
    return _ProteinSiteCandidate(
        residue=residue.upper(),
        protein_position=protein_position,
        probability=probability,
    )


def _probability_for_site_token(
    token: str,
    *,
    candidates: tuple[_LocalisationCandidate, ...],
    protein_start: int | None,
    ptmprophet_position_reference: str,
    row_position: int,
) -> float | None:
    if not candidates:
        return None
    protein_candidates = _convert_candidates_to_protein_sites(
        candidates,
        protein_start=protein_start,
        ptmprophet_position_reference=ptmprophet_position_reference,
        row_position=row_position,
    )
    for candidate in protein_candidates:
        if candidate.token == token:
            return candidate.probability
    if len(protein_candidates) == 1:
        return protein_candidates[0].probability
    return None


def _resolve_protein_start(
    row: Mapping[str, object],
    *,
    resolved: _ResolvedFragPipeColumns,
    row_position: int,
) -> int | None:
    if resolved.protein_start is None:
        return None
    return _normalise_positive_int(
        row[resolved.protein_start],
        field_name=f"FragPipe {resolved.protein_start}",
        row_position=row_position,
    )


def _resolve_required_column(
    columns: pd.Index,
    *,
    explicit: str | None,
    candidates: tuple[str, ...],
    field_name: str,
) -> str:
    resolved = _resolve_column(
        columns,
        explicit=explicit,
        candidates=candidates,
        field_name=field_name,
        required=True,
    )
    if resolved is None:  # pragma: no cover - _resolve_column raises first.
        raise PhosPyInputError(
            f"FragPipe importer could not infer {field_name}; configure an explicit "
            "column mapping."
        )
    return resolved


def _resolve_column(
    columns: pd.Index,
    *,
    explicit: str | None,
    candidates: tuple[str, ...],
    field_name: str,
    required: bool,
) -> str | None:
    override = validate_optional_fragpipe_column_name(
        explicit,
        field_name=field_name,
    )
    if override is not None:
        if override in columns:
            return override
        raise PhosPyInputError(f"{field_name}={override!r} is not present in source")
    for candidate in candidates:
        match = _find_column(columns, candidate)
        if match is not None:
            return match
    if not required:
        return None
    accepted = ", ".join(repr(candidate) for candidate in candidates)
    raise PhosPyInputError(
        f"FragPipe importer could not infer {field_name}; configure an explicit "
        f"column mapping. tried: {accepted}"
    )


def _find_column(columns: pd.Index, wanted: str) -> str | None:
    if wanted in columns:
        return str(wanted)
    normalised_wanted = _normalise_column_label(wanted)
    matches = [
        str(column)
        for column in columns.tolist()
        if _normalise_column_label(str(column)) == normalised_wanted
    ]
    if len(matches) > 1:
        joined = ", ".join(repr(item) for item in matches)
        raise PhosPyInputError(
            f"FragPipe importer found ambiguous source columns for {wanted!r}: {joined}"
        )
    if matches:
        return matches[0]
    return None


def _resolve_intensity_columns(
    source: pd.DataFrame,
    value: Mapping[str, str] | Sequence[str] | None,
    *,
    intensity_column_prefixes: Sequence[str],
) -> dict[str, str]:
    if value is not None:
        mapping = normalise_sample_column_mapping(value)
        missing = [column for column in mapping if column not in source.columns]
        if missing:
            joined = ", ".join(missing)
            raise PhosPyInputError(
                f"FragPipe intensity column mapping includes missing columns: {joined}"
            )
        return mapping
    if isinstance(intensity_column_prefixes, str) or not isinstance(
        intensity_column_prefixes,
        Sequence,
    ):
        raise PhosPyInputError(
            "fragpipe import request intensity_column_prefixes must be a sequence "
            "of strings"
        )
    prefixes = tuple(str(prefix) for prefix in intensity_column_prefixes)
    if not prefixes or any(prefix.strip() == "" for prefix in prefixes):
        raise PhosPyInputError(
            "fragpipe import request intensity_column_prefixes must contain "
            "non-empty strings"
        )
    mapping: dict[str, str] = {}
    for column in source.columns.astype(str).tolist():
        for prefix in prefixes:
            if not column.startswith(prefix):
                continue
            sample_id = column[len(prefix) :].strip() or column
            mapping[column] = sample_id
            break
    if not mapping:
        accepted = ", ".join(repr(prefix) for prefix in prefixes)
        raise PhosPyInputError(
            "FragPipe importer could not infer intensity columns. Configure "
            "FragPipeColumnMapping.intensity_columns or adjust "
            f"intensity_column_prefixes. tried prefixes: {accepted}"
        )
    return normalise_sample_column_mapping(mapping)


def _augment_mapped_result(
    mapped_result: PhosphositeImportResult,
    *,
    adapted: pd.DataFrame,
    flags: pd.DataFrame,
    contaminant_policy: str,
    decoy_policy: str,
    resolved: _ResolvedFragPipeColumns,
    filter_diagnostics: dict[str, object],
    adapter_diagnostics: dict[str, object],
    warnings: tuple[str, ...],
) -> PhosphositeImportResult:
    site_metadata = mapped_result.site_metadata_candidate
    peptide_evidence = mapped_result.peptide_evidence
    row_ids = adapted.loc[:, _ADAPTED_ROW_ID_COLUMN].astype(str).tolist()
    sidecars = adapted.loc[
        :,
        [
            _ADAPTED_CANDIDATE_SITES_COLUMN,
            _ADAPTED_SITE_PROBABILITIES_COLUMN,
            _ADAPTED_AMBIGUOUS_COLUMN,
            _ADAPTED_MODIFIED_PHOSPHO_COUNT_COLUMN,
        ],
    ].copy(deep=True)
    sidecars.index = pd.Index(row_ids, name=_ADAPTED_ROW_ID_COLUMN)
    for column_name in sidecars.columns.astype(str).tolist():
        site_metadata[column_name] = sidecars.loc[
            site_metadata.index,
            column_name,
        ].tolist()
        if peptide_evidence is not None:
            peptide_evidence[column_name] = sidecars.loc[
                peptide_evidence.index,
                column_name,
            ].tolist()

    flag_values = flags.copy(deep=True)
    flag_values.index = pd.Index(row_ids, name=_ADAPTED_ROW_ID_COLUMN)
    if contaminant_policy == FRAGPIPE_FLAG_POLICY_FLAG:
        site_metadata[_FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN] = (
            flag_values.loc[site_metadata.index, _FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN]
            .astype(bool)
            .tolist()
        )
        if peptide_evidence is not None:
            peptide_evidence[_FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN] = (
                flag_values.loc[
                    peptide_evidence.index,
                    _FRAGPIPE_CONTAMINANT_OUTPUT_COLUMN,
                ]
                .astype(bool)
                .tolist()
            )
    if decoy_policy == FRAGPIPE_FLAG_POLICY_FLAG:
        site_metadata[_FRAGPIPE_DECOY_OUTPUT_COLUMN] = (
            flag_values.loc[site_metadata.index, _FRAGPIPE_DECOY_OUTPUT_COLUMN]
            .astype(bool)
            .tolist()
        )
        if peptide_evidence is not None:
            peptide_evidence[_FRAGPIPE_DECOY_OUTPUT_COLUMN] = (
                flag_values.loc[peptide_evidence.index, _FRAGPIPE_DECOY_OUTPUT_COLUMN]
                .astype(bool)
                .tolist()
            )

    diagnostics = dict(mapped_result.diagnostics)
    diagnostics["fragpipe"] = {
        "source_type": "fragpipe_ptmprophet_phosphosite",
        "resolved_columns": _resolved_columns_payload(resolved),
        "filtering": filter_diagnostics,
        "adaptation": adapter_diagnostics,
    }
    return PhosphositeImportResult(
        phospho_matrix_candidate=mapped_result.phospho_matrix_candidate,
        site_metadata_candidate=site_metadata,
        peptide_evidence=peptide_evidence,
        sample_column_mapping=mapped_result.sample_column_mapping,
        localisation_confidence_column=mapped_result.localisation_confidence_column,
        warnings=tuple(dict.fromkeys((*mapped_result.warnings, *warnings))),
        diagnostics=diagnostics,
        source_name=mapped_result.source_name,
        _assume_owned=True,
    )


def _parse_probability(value: object, *, field_name: str, row_position: int) -> float:
    if isinstance(value, bool):
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} contains boolean probability"
        )
    if isinstance(value, str):
        try:
            probability = float(value.strip())
        except ValueError as exc:
            raise PhosPyInputError(
                f"{field_name} row_position={row_position} contains non-numeric "
                f"probability {value!r}"
            ) from exc
    elif isinstance(value, (int, float)):
        probability = float(value)
    else:
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} contains non-numeric "
            f"probability {value!r}"
        )
    if not math.isfinite(probability) or probability < 0.0 or probability > 1.0:
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} probability must be in "
            f"[0.0, 1.0]; value={value!r}"
        )
    return float(probability)


def _deduplicate_localisation_candidates(
    candidates: Sequence[_LocalisationCandidate],
) -> tuple[_LocalisationCandidate, ...]:
    values: dict[tuple[str, int | None], _LocalisationCandidate] = {}
    for candidate in candidates:
        key = (candidate.residue.upper(), candidate.position)
        current = values.get(key)
        if current is None:
            values[key] = candidate
            continue
        if current.probability is None:
            values[key] = candidate
            continue
        if candidate.probability is not None and candidate.probability > (
            current.probability
        ):
            values[key] = candidate
    return tuple(values.values())


def _deduplicate_protein_candidates(
    candidates: Sequence[_ProteinSiteCandidate],
) -> tuple[_ProteinSiteCandidate, ...]:
    values: dict[str, _ProteinSiteCandidate] = {}
    for candidate in candidates:
        current = values.get(candidate.token)
        if current is None:
            values[candidate.token] = candidate
            continue
        if current.probability is None:
            values[candidate.token] = candidate
            continue
        if candidate.probability is not None and candidate.probability > (
            current.probability
        ):
            values[candidate.token] = candidate
    return tuple(values.values())


def _format_site_probabilities(
    candidates: tuple[_ProteinSiteCandidate, ...],
) -> str:
    parts: list[str] = []
    for candidate in _deduplicate_protein_candidates(candidates):
        if candidate.probability is None:
            parts.append(f"{candidate.token}:NA")
            continue
        parts.append(f"{candidate.token}:{float(candidate.probability):.6g}")
    return ";".join(parts)


def _next_modified_site_with_residue(
    sites: tuple[ModifiedResidue, ...],
    *,
    residue: str,
    start: int,
    row_position: int,
    field_name: str,
) -> ModifiedResidue:
    residue = residue.upper()
    for site in sites[start:]:
        if site.residue == residue:
            return site
    raise PhosPyInputError(
        f"{field_name} row_position={row_position} could not align residue-only "
        f"PTMProphet probability for residue {residue!r} to the modified peptide"
    )


def _resolve_flag_series(
    source: pd.DataFrame,
    *,
    column: str | None,
    field_name: str,
) -> pd.Series | None:
    if column is None:
        return None
    return pd.Series(
        [
            _parse_fragpipe_flag(
                value,
                field_name=field_name,
                row_position=position,
            )
            for position, value in enumerate(source.loc[:, column].tolist())
        ],
        index=source.index.copy(),
        dtype=bool,
    )


def _raise_for_forbidden_flags(
    values: pd.Series,
    *,
    policy: str,
    label: str,
) -> None:
    if policy != FRAGPIPE_FLAG_POLICY_ERROR:
        return
    count = int(values.astype(bool).sum())
    if not count:
        return
    raise PhosPyInputError(
        f"FragPipe importer encountered {count} {label} row(s) under policy='error'"
    )


def _parse_fragpipe_flag(
    value: object,
    *,
    field_name: str,
    row_position: int,
) -> bool:
    if _is_missing(value):
        return False
    if isinstance(value, bool):
        return bool(value)
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        numeric = float(value)
        if numeric == 0.0:
            return False
        if numeric == 1.0:
            return True
    token = str(value).strip().lower()
    if token in {"", "-", "false", "f", "no", "n", "0"}:
        return False
    if token in {"+", "true", "t", "yes", "y", "1"}:
        return True
    raise PhosPyInputError(
        f"{field_name} row_position={row_position} contains unsupported flag "
        f"value {value!r}; expected '+', blank, boolean, or 0/1 style values"
    )


def _parse_protein_accession(
    value: object,
    *,
    field_name: str,
    row_position: int,
) -> str:
    token = _first_list_token(value, field_name=field_name, row_position=row_position)
    cleaned = _strip_protein_prefixes(token)
    parts = cleaned.split("|")
    if len(parts) >= 3 and parts[0].strip().lower() in {"sp", "tr"}:
        accession = parts[1].strip()
    else:
        accession = cleaned.split()[0].strip()
    if accession == "":
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} did not contain a protein "
            "accession after parsing"
        )
    return accession


def _strip_protein_prefixes(value: str) -> str:
    token = value.strip()
    changed = True
    while changed:
        changed = False
        upper = token.upper()
        for prefix in (*_CONTAMINANT_PREFIXES, *_DECOY_PREFIXES):
            if upper.startswith(prefix):
                token = token[len(prefix) :].strip()
                changed = True
                break
    return token


def _has_any_prefixed_token(value: object, *, prefixes: tuple[str, ...]) -> bool:
    for token in _split_multi_value(value):
        upper = token.upper()
        if any(upper.startswith(prefix) for prefix in prefixes):
            return True
    return False


def _multi_value_count(value: object) -> int:
    return len(_split_multi_value(value))


def _split_multi_value(value: object) -> list[str]:
    if _is_missing(value):
        return []
    if isinstance(value, str):
        stripped = value.strip()
        if stripped == "":
            return []
        return [
            token.strip()
            for token in _MULTI_VALUE_SPLIT_PATTERN.split(stripped)
            if token.strip()
        ]
    return [str(value).strip()]


def _first_list_token(
    value: object,
    *,
    field_name: str,
    row_position: int,
) -> str:
    tokens = _split_multi_value(value)
    if not tokens:
        raise PhosPyInputError(
            f"{field_name} must contain non-empty values; row_position={row_position}"
        )
    return tokens[0]


def _required_text(
    value: object,
    *,
    field_name: str,
    row_position: int,
) -> str:
    if _is_missing(value):
        raise PhosPyInputError(
            f"{field_name} must not contain missing values; row_position={row_position}"
        )
    token = str(value).strip()
    if token == "":
        raise PhosPyInputError(
            f"{field_name} must contain non-empty values; row_position={row_position}"
        )
    return token


def _optional_text(value: object) -> str | None:
    if _is_missing(value):
        return None
    token = str(value).strip()
    if token == "":
        return None
    return token


def _normalise_positive_int(
    value: object,
    *,
    field_name: str,
    row_position: int,
) -> int:
    token = _required_text(value, field_name=field_name, row_position=row_position)
    try:
        numeric = float(token)
    except ValueError as exc:
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} must be a positive integer; "
            f"value={value!r}"
        ) from exc
    if not math.isfinite(numeric) or not numeric.is_integer():
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} must be a positive integer; "
            f"value={value!r}"
        )
    integer = int(numeric)
    if integer < 1:
        raise PhosPyInputError(
            f"{field_name} row_position={row_position} must be a positive integer; "
            f"value={value!r}"
        )
    return integer


def _build_row_ids(
    *,
    source: pd.DataFrame,
    resolved: _ResolvedFragPipeColumns,
    protein_values: list[str],
    site_values: list[str],
    source_row_numbers: list[int],
) -> list[str]:
    if resolved.row_id is not None:
        return [
            _required_text(
                value,
                field_name=f"FragPipe {resolved.row_id}",
                row_position=position,
            )
            for position, value in enumerate(source.loc[:, resolved.row_id])
        ]
    return [
        f"fragpipe:{protein}:{site}:row{row_number}"
        for protein, site, row_number in zip(
            protein_values,
            site_values,
            source_row_numbers,
            strict=True,
        )
    ]


def _build_unique_feature_ids(
    *,
    source: pd.DataFrame,
    resolved: _ResolvedFragPipeColumns,
    source_row_numbers: list[int],
) -> list[str]:
    if resolved.unique_feature_id is not None:
        return [
            _required_text(
                value,
                field_name=f"FragPipe {resolved.unique_feature_id}",
                row_position=position,
            )
            for position, value in enumerate(source.loc[:, resolved.unique_feature_id])
        ]
    return [f"fragpipe_feature_{row_number}" for row_number in source_row_numbers]


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


def _require_non_empty_unique_columns(source: pd.DataFrame) -> None:
    if not isinstance(source, pd.DataFrame):
        raise PhosPyInputError("FragPipe import source must be a pandas DataFrame")
    if source.empty:
        raise PhosPyInputError("FragPipe import source must not be empty")
    column_labels = [str(column) for column in source.columns.tolist()]
    if len(set(column_labels)) != len(column_labels):
        raise PhosPyInputError("FragPipe import source columns must be unique")


def _normalise_column_label(value: str) -> str:
    return " ".join(str(value).strip().lower().split())


def _is_missing(value: object) -> bool:
    try:
        return bool(pd.Series((value,), dtype="object").isna().iat[0])
    except (TypeError, ValueError):
        return False


__all__ = [
    "FragPipeColumnMapping",
    "FragPipePTMProphetImportRequest",
    "FragPipePTMProphetImporter",
]
