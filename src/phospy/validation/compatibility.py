from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .errors import InputCompatibilityError
from .frames import require_columns
from .identifiers import normalize_identifier_series
from .scalars import validate_fraction
from .schemas import PredictionScoreMatrixSchema, PredMatSchema, SiteMatrixSchema


@dataclass(frozen=True, slots=True)
class ProteinCorrectionMatchSummary:
    """Describe phosphosite-to-protein matching before correction."""

    input_rows: int
    matched_rows: int
    unmatched_rows: int
    unmatched_fraction: float
    unmatched_gene_preview: tuple[str, ...]


def validate_core_column_alignment(
    total_cols: Sequence[str],
    phospho_cols: Sequence[str],
    corrected_cols: Sequence[str] | None = None,
    *,
    context: str = "Core preprocessing inputs",
) -> None:
    """Validate that paired preprocessing column groups align by width."""

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
) -> ProteinCorrectionMatchSummary:
    """Validate phosphosite/protein correction inputs and matching coverage."""

    resolved_max_unmatched_fraction = validate_fraction(
        max_unmatched_fraction,
        name="max_unmatched_fraction",
    )
    if len(phospho_cols) != len(protein_cols):
        msg = f"{context} require the same number of phospho and protein columns"
        raise InputCompatibilityError(msg)
    if phospho_df.empty:
        msg = f"{context} contain no phosphosite rows after filtering"
        raise InputCompatibilityError(msg)
    if total_df.empty:
        msg = f"{context} contain no protein rows after filtering"
        raise InputCompatibilityError(msg)

    require_columns(
        phospho_df,
        required_columns=[phospho_gene_col, *phospho_cols],
        context=f"{context} phospho input",
    )
    require_columns(
        total_df,
        required_columns=[total_gene_col, *protein_cols],
        context=f"{context} total input",
    )

    total_gene_series = normalize_identifier_series(total_df[total_gene_col])
    if total_gene_series.duplicated().any():
        msg = (
            f"{context}: {total_gene_col} must be unique before protein "
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

    input_rows = int(len(phospho_df))
    unmatched_rows = input_rows - matched_rows
    unmatched_fraction = unmatched_rows / input_rows
    unmatched_genes = pd.unique(phospho_genes.loc[~matched_mask].dropna())
    unmatched_gene_preview = tuple(str(gene) for gene in unmatched_genes[:5])

    if unmatched_rows > 0 and unmatched_fraction > resolved_max_unmatched_fraction:
        unmatched_preview = ", ".join(unmatched_gene_preview)
        percent = unmatched_fraction * 100.0
        msg = (
            f"{context} would drop {unmatched_rows} of {input_rows} phosphosite "
            f"rows ({percent:.1f}%) due to missing protein matches in "
            f"{total_gene_col}: {unmatched_preview}"
        )
        raise InputCompatibilityError(msg)

    return ProteinCorrectionMatchSummary(
        input_rows=input_rows,
        matched_rows=matched_rows,
        unmatched_rows=unmatched_rows,
        unmatched_fraction=unmatched_fraction,
        unmatched_gene_preview=unmatched_gene_preview,
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
    """Validate overlap between a prediction matrix and a phosphosite matrix."""

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


def validate_workflow_matrix_inputs(
    phospho_matrix: pd.DataFrame,
    substrate_map: Mapping[str, Sequence[str]],
    site_sequences: Mapping[str, str] | pd.Series | None,
    *,
    require_site_sequences_for_prediction: bool,
    context: str,
) -> tuple[pd.DataFrame, tuple[str, ...]]:
    """Validate workflow matrices and resolve the scoring subset."""

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

    scoring_site_index = tuple(str(site) for site in validated_matrix.index)

    if require_site_sequences_for_prediction:
        if site_sequences is None:
            msg = "site_sequences are required when motif_sequences are provided"
            raise InputCompatibilityError(msg)

        sequence_index = _extract_sequence_index(site_sequences)
        scoring_site_index = tuple(
            site for site in validated_matrix.index if str(site) in sequence_index
        )
        if not scoring_site_index:
            msg = (
                f"{context} contain no phosphosites with sequence coverage "
                "required for scoring and prediction"
            )
            raise InputCompatibilityError(msg)

    return validated_matrix, scoring_site_index


def validate_signalome_alignment(
    *,
    scoring_matrix: pd.DataFrame,
    pred_mat: pd.DataFrame,
    expression_matrix: pd.DataFrame,
    kinases_of_interest: Sequence[str],
    module_count: int | None,
    scoring_context: str,
    pred_mat_context: str,
    expression_context: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, tuple[str, ...]]:
    """Validate and align signalome inputs across sites and kinase columns."""

    validated_scoring_matrix = PredictionScoreMatrixSchema.validate(
        scoring_matrix,
        context=scoring_context,
    )
    validated_pred_mat = PredMatSchema.validate(
        pred_mat,
        context=pred_mat_context,
    )
    _require_finite_pred_mat(validated_pred_mat, context=pred_mat_context)
    validated_expression_matrix = SiteMatrixSchema.validate(
        expression_matrix,
        context=expression_context,
    )

    common_sites = tuple(
        site_id
        for site_id in validated_scoring_matrix.index.astype(str)
        if site_id in validated_pred_mat.index
        and site_id in validated_expression_matrix.index
    )
    if not common_sites:
        msg = (
            f"{scoring_context}, {pred_mat_context}, and {expression_context} "
            "must share at least one phosphosite row"
        )
        raise InputCompatibilityError(msg)

    common_kinases = tuple(
        kinase
        for kinase in validated_scoring_matrix.columns.astype(str)
        if kinase in validated_pred_mat.columns
    )
    if not common_kinases:
        msg = (
            f"{scoring_context} and {pred_mat_context} must share at least one "
            "kinase column"
        )
        raise InputCompatibilityError(msg)

    missing_koi = [
        kinase for kinase in kinases_of_interest if kinase not in common_kinases
    ]
    if missing_koi:
        missing = ", ".join(missing_koi)
        msg = (
            "kinases_of_interest are not available in the aligned signalome "
            f"inputs: {missing}"
        )
        raise InputCompatibilityError(msg)

    if module_count is not None and module_count > len(common_sites):
        msg = "module_count cannot exceed the number of aligned phosphosite rows"
        raise InputCompatibilityError(msg)

    return (
        validated_scoring_matrix.loc[list(common_sites), list(common_kinases)],
        validated_pred_mat.loc[list(common_sites), list(common_kinases)],
        validated_expression_matrix.loc[list(common_sites)],
        common_sites,
    )


def _extract_sequence_index(site_sequences: Mapping[str, str] | pd.Series) -> set[str]:
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


def _require_finite_pred_mat(
    frame: pd.DataFrame,
    *,
    context: str,
) -> None:
    failures: list[str] = []
    for column in frame.columns.astype(str):
        series = frame.loc[:, column]
        invalid_mask = ~np.isfinite(series.to_numpy(dtype=float))
        if invalid_mask.any():
            sample_values = series.loc[invalid_mask].astype(str).unique()[:3]
            sample_preview = ", ".join(str(value) for value in sample_values)
            failures.append(f"{column} ({sample_preview})")
    if failures:
        failures_str = "; ".join(failures)
        msg = f"{context} contains non-finite values in numeric columns: {failures_str}"
        raise InputCompatibilityError(msg)


__all__ = [
    "ProteinCorrectionMatchSummary",
    "validate_core_column_alignment",
    "validate_pred_mat_overlap",
    "validate_protein_correction_inputs",
    "validate_signalome_alignment",
    "validate_workflow_matrix_inputs",
]
