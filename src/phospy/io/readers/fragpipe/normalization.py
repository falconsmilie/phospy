"""FragPipe-specific normalization into mapped-importer input columns."""
# pyright: reportMissingTypeStubs=false

from __future__ import annotations

import math

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.io.readers._table_parsing import (
    build_row_ids,
    build_unique_feature_ids,
    first_list_token,
    is_missing,
    multi_value_count,
    optional_text,
    required_text,
)
from phospy.io.readers.fragpipe.columns import _resolved_columns_payload
from phospy.io.readers.fragpipe.constants import (
    _ADAPTED_AMBIGUOUS_COLUMN,
    _ADAPTED_CANDIDATE_SITES_COLUMN,
    _ADAPTED_GENE_SYMBOL_COLUMN,
    _ADAPTED_LOCALISATION_COLUMN,
    _ADAPTED_MODIFIED_PEPTIDE_SEQUENCE_COLUMN,
    _ADAPTED_MODIFIED_PHOSPHO_COUNT_COLUMN,
    _ADAPTED_PEPTIDE_SEQUENCE_COLUMN,
    _ADAPTED_PEPTIDE_SITE_STRING_COLUMN,
    _ADAPTED_PROTEIN_ACCESSION_COLUMN,
    _ADAPTED_PROTEIN_ID_COLUMN,
    _ADAPTED_QUANTIFICATION_DEPTH_COLUMN,
    _ADAPTED_ROW_ID_COLUMN,
    _ADAPTED_SITE_COLUMN,
    _ADAPTED_SITE_PROBABILITIES_COLUMN,
    _ADAPTED_SITE_SEQUENCE_COLUMN,
    _ADAPTED_UNIQUE_FEATURE_ID_COLUMN,
)
from phospy.io.readers.fragpipe.conversion import (
    _parse_protein_accession,
    _resolve_site_call,
)
from phospy.io.readers.fragpipe.models import _ResolvedFragPipeColumns
from phospy.science.differential.quantification_depth import (
    QUANTIFICATION_DEPTH_COLUMN,
    QUANTIFICATION_DEPTH_INTEGER_TOLERANCE,
)
from phospy.science.evidence.modified_peptides import parse_modified_peptide_sequence


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
        if multi_value_count(row_lookup[resolved.protein_accession]) > 1:
            protein_group_rows += 1
        protein_values.append(
            _parse_protein_accession(
                row_lookup[resolved.protein_accession],
                field_name=f"FragPipe {resolved.protein_accession}",
                row_position=position,
            )
        )
        gene_values.append(
            first_list_token(
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
        peptide_sequence = required_text(
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
    adapted[_ADAPTED_LOCALISATION_COLUMN] = pd.Series(
        localisation_values,
        index=adapted.index,
        dtype=object,
    )
    adapted[_ADAPTED_UNIQUE_FEATURE_ID_COLUMN] = _build_unique_feature_ids(
        source=source,
        resolved=resolved,
        source_row_numbers=source_row_numbers,
    )
    adapted[_ADAPTED_PEPTIDE_SEQUENCE_COLUMN] = [
        required_text(
            value,
            field_name=f"FragPipe {resolved.peptide_sequence}",
            row_position=position,
        )
        for position, value in enumerate(source.loc[:, resolved.peptide_sequence])
    ]
    modified_peptide_sequences = pd.Series(
        [
            required_text(
                value,
                field_name=f"FragPipe {resolved.modified_peptide_sequence}",
                row_position=position,
            )
            for position, value in enumerate(
                source.loc[:, resolved.modified_peptide_sequence]
            )
        ],
        index=adapted.index,
        dtype=object,
    )
    adapted[_ADAPTED_MODIFIED_PEPTIDE_SEQUENCE_COLUMN] = modified_peptide_sequences
    adapted[_ADAPTED_PEPTIDE_SITE_STRING_COLUMN] = peptide_site_values
    adapted[_ADAPTED_CANDIDATE_SITES_COLUMN] = candidate_site_values
    adapted[_ADAPTED_SITE_PROBABILITIES_COLUMN] = site_probability_values
    adapted[_ADAPTED_AMBIGUOUS_COLUMN] = ambiguous_values
    adapted[_ADAPTED_MODIFIED_PHOSPHO_COUNT_COLUMN] = phospho_counts
    if resolved.site_sequence is not None:
        adapted[_ADAPTED_SITE_SEQUENCE_COLUMN] = [
            optional_text(value)
            for value in source.loc[:, resolved.site_sequence].tolist()
        ]
    quantification_depth_diagnostics = _adapt_quantification_depth(
        adapted=adapted,
        source=source,
        resolved=resolved,
    )

    diagnostics = {
        "resolved_columns": _resolved_columns_payload(resolved),
        "ptmprophet_position_reference": ptmprophet_position_reference,
        "quantification_depth": quantification_depth_diagnostics,
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


def _build_row_ids(
    *,
    source: pd.DataFrame,
    resolved: _ResolvedFragPipeColumns,
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
        importer_label="FragPipe",
        generated_prefix="fragpipe",
    )


def _adapt_quantification_depth(
    *,
    adapted: pd.DataFrame,
    source: pd.DataFrame,
    resolved: _ResolvedFragPipeColumns,
) -> dict[str, object]:
    if resolved.quantification_depth is None:
        return {
            "status": "not_mapped",
            "source_column": None,
            "output_column": QUANTIFICATION_DEPTH_COLUMN,
            "quantification_depth_kind": None,
        }
    values = _normalise_quantification_depth_source_values(
        source.loc[:, resolved.quantification_depth],
        source_column=resolved.quantification_depth,
    )
    adapted.loc[:, _ADAPTED_QUANTIFICATION_DEPTH_COLUMN] = pd.Series(
        values,
        index=adapted.index.copy(),
        dtype=float,
    )
    return {
        "status": "reported",
        "source_column": resolved.quantification_depth,
        "adapted_column": _ADAPTED_QUANTIFICATION_DEPTH_COLUMN,
        "output_column": QUANTIFICATION_DEPTH_COLUMN,
        "quantification_depth_kind": resolved.quantification_depth_kind,
        "row_count": int(len(values)),
        "missing_count": 0,
        "minimum_count": float(min(values)),
        "maximum_count": float(max(values)),
    }


def _normalise_quantification_depth_source_values(
    series: pd.Series,
    *,
    source_column: str,
) -> list[float]:
    return [
        _normalise_quantification_depth_value(
            value,
            source_column=source_column,
            row_position=position,
        )
        for position, value in enumerate(series.tolist())
    ]


def _normalise_quantification_depth_value(
    value: object,
    *,
    source_column: str,
    row_position: int,
) -> float:
    field_name = f"FragPipe {source_column} quantification_depth"
    if _is_missing_quantification_depth_value(value):
        raise PhosPyInputError(
            f"{field_name} must not contain missing values; row_position={row_position}"
        )
    if isinstance(value, bool) or type(value).__name__ == "bool_":
        raise PhosPyInputError(
            f"{field_name} must contain numeric count values, not booleans; "
            f"row_position={row_position}, offending_value={value!r}"
        )
    try:
        numeric = float(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise PhosPyInputError(
            f"{field_name} must contain numeric count values; "
            f"row_position={row_position}, offending_value={value!r}"
        ) from exc
    if not math.isfinite(numeric):
        raise PhosPyInputError(
            f"{field_name} must contain finite numeric count values; "
            f"row_position={row_position}, offending_value={value!r}"
        )
    if numeric < 1.0:
        raise PhosPyInputError(
            f"{field_name} values must be >= 1; "
            f"row_position={row_position}, offending_value={value!r}"
        )
    if not math.isclose(
        numeric,
        round(numeric),
        rel_tol=0.0,
        abs_tol=QUANTIFICATION_DEPTH_INTEGER_TOLERANCE,
    ):
        raise PhosPyInputError(
            f"{field_name} count values must be integer-valued within tolerance "
            f"{QUANTIFICATION_DEPTH_INTEGER_TOLERANCE:g}; "
            f"row_position={row_position}, offending_value={value!r}"
        )
    return float(round(numeric))


def _is_missing_quantification_depth_value(value: object) -> bool:
    if is_missing(value):
        return True
    if isinstance(value, str):
        return value.strip().lower() in {"", "na", "n/a", "nan", "null"}
    return False


def _build_unique_feature_ids(
    *,
    source: pd.DataFrame,
    resolved: _ResolvedFragPipeColumns,
    source_row_numbers: list[int],
) -> list[str]:
    return build_unique_feature_ids(
        source=source,
        explicit_column=resolved.unique_feature_id,
        source_row_numbers=source_row_numbers,
        importer_label="FragPipe",
        generated_prefix="fragpipe",
    )


__all__ = ["_adapt_fragpipe_source"]
