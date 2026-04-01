from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

from .errors import InputCompatibilityError
from .normalization import normalize_identifier_series
from .requests import KinaseWorkflowRequest
from .tables import PredMatSchema, SiteMatrixSchema


def validate_core_column_alignment(
    total_cols: Sequence[str],
    phospho_cols: Sequence[str],
    corrected_cols: Sequence[str] | None = None,
    *,
    context: str = "Core preprocessing inputs",
) -> None:
    if len(total_cols) != len(phospho_cols):
        msg = f"{context} require the same number of total and phospho value columns"
        raise InputCompatibilityError(msg)
    if corrected_cols is not None and len(corrected_cols) != len(total_cols):
        msg = (
            f"{context} require corrected value columns to align with total and "
            "phospho value columns"
        )
        raise InputCompatibilityError(msg)


def validate_protein_correction_inputs(
    phospho_df: pd.DataFrame,
    total_df: pd.DataFrame,
    *,
    phospho_gene_col: str,
    total_gene_col: str,
    phospho_cols: Sequence[str],
    protein_cols: Sequence[str],
    max_unmatched_fraction: float = 0.0,
    context: str = "Protein correction inputs",
) -> None:
    if len(phospho_cols) != len(protein_cols):
        msg = f"{context} require the same number of phospho and protein columns"
        raise InputCompatibilityError(msg)
    if phospho_df.empty:
        msg = f"{context} contain no phosphosite rows after filtering"
        raise InputCompatibilityError(msg)
    if total_df.empty:
        msg = f"{context} contain no protein rows after filtering"
        raise InputCompatibilityError(msg)
    if phospho_gene_col not in phospho_df.columns:
        msg = f"{context} is missing phospho gene column: {phospho_gene_col}"
        raise InputCompatibilityError(msg)
    if total_gene_col not in total_df.columns:
        msg = f"{context} is missing total gene column: {total_gene_col}"
        raise InputCompatibilityError(msg)

    total_gene_series = normalize_identifier_series(total_df[total_gene_col])
    if total_gene_series.duplicated().any():
        msg = (
            f"{context} require unique values in {total_gene_col} before protein "
            "correction to avoid duplicating phosphosite rows"
        )
        raise InputCompatibilityError(msg)

    phospho_genes = normalize_identifier_series(phospho_df[phospho_gene_col])
    total_gene_values = set(total_gene_series)
    matched_mask = phospho_genes.isin(total_gene_values)
    matched_rows = int(matched_mask.sum())

    if matched_rows == 0:
        msg = (
            f"{context} have no overlapping gene identifiers between "
            f"{phospho_gene_col} and {total_gene_col}"
        )
        raise InputCompatibilityError(msg)

    unmatched_rows = len(phospho_df) - matched_rows
    if unmatched_rows == 0:
        return

    unmatched_fraction = unmatched_rows / len(phospho_df)
    if unmatched_fraction > max_unmatched_fraction:
        unmatched_genes = pd.unique(phospho_genes.loc[~matched_mask].dropna())
        unmatched_preview = ", ".join(str(gene) for gene in unmatched_genes[:5])
        percent = unmatched_fraction * 100.0
        msg = (
            f"{context} would drop {unmatched_rows} of {len(phospho_df)} phosphosite "
            f"rows ({percent:.1f}%) due to missing protein matches in "
            f"{total_gene_col}: {unmatched_preview}"
        )
        raise InputCompatibilityError(msg)


def validate_kinase_activity_inputs(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    *,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
    min_overlap: int = 1,
    min_fraction: float = 0.1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    validated_pred_mat = PredMatSchema.validate(pred_mat, context=pred_context)
    validated_matrix = SiteMatrixSchema.validate(phospho_matrix, context=matrix_context)
    validate_pred_mat_overlap(
        validated_pred_mat,
        validated_matrix,
        pred_context=pred_context,
        matrix_context=matrix_context,
        min_overlap=min_overlap,
        min_fraction=min_fraction,
    )
    return validated_pred_mat, validated_matrix


def validate_workflow_request(
    request: KinaseWorkflowRequest,
    *,
    context: str = "Kinase workflow inputs",
) -> pd.DataFrame:
    return validate_workflow_inputs(
        request.phospho_matrix,
        request.substrate_map,
        request.site_sequences,
        request.motif_sequences,
        context=context,
    )


def validate_pred_mat_overlap(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    *,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
    min_overlap: int = 1,
    min_fraction: float = 0.1,
) -> None:
    overlap = pred_mat.index.intersection(phospho_matrix.index)
    overlap_count = len(overlap)
    if overlap_count == 0:
        msg = f"{pred_context} and {matrix_context} have no overlapping phosphosite IDs"
        raise InputCompatibilityError(msg)

    matrix_rows = len(phospho_matrix.index)
    overlap_fraction = overlap_count / max(matrix_rows, 1)
    if overlap_count < min_overlap or overlap_fraction < min_fraction:
        percent = overlap_fraction * 100.0
        msg = (
            f"{pred_context} and {matrix_context} have insufficient overlapping "
            f"phosphosite IDs: {overlap_count} row(s) ({percent:.1f}%)"
        )
        raise InputCompatibilityError(msg)


def validate_workflow_inputs(
    phospho_matrix: pd.DataFrame,
    substrate_map: Mapping[str, Sequence[str]],
    site_sequences: Mapping[str, str] | pd.Series | None,
    motif_sequences: Mapping[str, Sequence[str]] | None,
    *,
    context: str = "Kinase workflow inputs",
) -> pd.DataFrame:
    validated_matrix = SiteMatrixSchema.validate(
        phospho_matrix,
        context="phospho_matrix",
    )

    overlapping_sites = {
        site
        for sites in substrate_map.values()
        for site in sites
        if site in validated_matrix.index
    }
    if not overlapping_sites:
        msg = f"{context} contain no overlap between substrate_map and phospho_matrix"
        raise InputCompatibilityError(msg)

    if site_sequences is not None:
        sequence_index = _extract_sequence_index(site_sequences)
        missing = [
            site for site in validated_matrix.index if site not in sequence_index
        ]
        if missing:
            missing_preview = ", ".join(missing[:5])
            msg = (
                f"site_sequences is missing entries for phosphosites: {missing_preview}"
            )
            raise InputCompatibilityError(msg)

    if motif_sequences is not None:
        widths: set[int] = set()
        for kinase, sequences in motif_sequences.items():
            if not sequences:
                continue
            kinase_widths = {len(str(sequence)) for sequence in sequences}
            if len(kinase_widths) != 1:
                msg = (
                    f"motif_sequences for kinase {kinase} must use a consistent "
                    "sequence width"
                )
                raise InputCompatibilityError(msg)
            widths.update(kinase_widths)
        if len(widths) > 1:
            msg = "motif_sequences must use the same sequence width across kinases"
            raise InputCompatibilityError(msg)

    return validated_matrix


def _extract_sequence_index(
    site_sequences: Mapping[str, str] | pd.Series,
) -> set[str]:
    if isinstance(site_sequences, pd.Series):
        return {str(value) for value in site_sequences.index}
    if isinstance(site_sequences, Mapping):
        return {str(value) for value in site_sequences}
    msg = (
        "site_sequences must be provided as a mapping keyed by phosphosite ID "
        "or as a pandas Series with an explicit phosphosite index; plain "
        "sequences are not supported"
    )
    raise InputCompatibilityError(msg)
