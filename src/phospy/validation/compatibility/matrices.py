from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd

from ...errors import InputCompatibilityError, format_overlap_failure_message
from ..schema.tables import (
    PredictionScoreMatrixSchema,
    PredMatForSignalomeSchema,
    SiteMatrixSchema,
)

DEFAULT_MIN_PRED_MAT_OVERLAP = 1
DEFAULT_MIN_PRED_MAT_OVERLAP_FRACTION = 0.5


@dataclass(frozen=True, slots=True)
class PredMatOverlapSummary:
    """Resolved overlap details between a predMat and a phosphosite matrix."""

    overlap_count: int
    matrix_rows: int
    pred_mat_rows: int
    pred_context: str
    matrix_context: str

    @property
    def matrix_fraction(self) -> float:
        return self.overlap_count / max(self.matrix_rows, 1)

    @property
    def pred_mat_fraction(self) -> float:
        return self.overlap_count / max(self.pred_mat_rows, 1)

    @property
    def is_partial(self) -> bool:
        return (
            self.overlap_count != self.matrix_rows
            or self.overlap_count != self.pred_mat_rows
        )

    @property
    def message(self) -> str:
        return (
            f"{self.pred_context} and {self.matrix_context} partially overlap: using "
            f"{self.overlap_count} shared phosphosite row(s) "
            f"({self.matrix_fraction * 100.0:.1f}% of {self.matrix_context}; "
            f"{self.pred_mat_fraction * 100.0:.1f}% of {self.pred_context})"
        )


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


def validate_pred_mat_overlap(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    *,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
    min_overlap: int = DEFAULT_MIN_PRED_MAT_OVERLAP,
    min_fraction: float = DEFAULT_MIN_PRED_MAT_OVERLAP_FRACTION,
) -> PredMatOverlapSummary:
    """Validate overlap between a prediction matrix and a phosphosite matrix."""

    matrix_rows = len(phospho_matrix.index)
    pred_mat_rows = len(pred_mat.index)
    overlap = pred_mat.index.intersection(phospho_matrix.index)
    overlap_count = len(overlap)
    if overlap_count == 0:
        msg = format_overlap_failure_message(
            pred_context=pred_context,
            matrix_context=matrix_context,
            overlap_count=overlap_count,
            pred_mat_rows=pred_mat_rows,
            matrix_rows=matrix_rows,
            min_overlap=min_overlap,
            min_fraction=min_fraction,
        )
        raise InputCompatibilityError(msg)

    overlap_fraction = overlap_count / max(matrix_rows, 1)
    if overlap_count < min_overlap or overlap_fraction < min_fraction:
        msg = format_overlap_failure_message(
            pred_context=pred_context,
            matrix_context=matrix_context,
            overlap_count=overlap_count,
            pred_mat_rows=pred_mat_rows,
            matrix_rows=matrix_rows,
            min_overlap=min_overlap,
            min_fraction=min_fraction,
        )
        raise InputCompatibilityError(msg)

    return PredMatOverlapSummary(
        overlap_count=overlap_count,
        matrix_rows=matrix_rows,
        pred_mat_rows=pred_mat_rows,
        pred_context=pred_context,
        matrix_context=matrix_context,
    )


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

    unique_substrate_sites = {
        site for sites in substrate_map.values() for site in sites
    }
    overlapping_sites = {
        site for site in unique_substrate_sites if site in validated_matrix.index
    }
    if not overlapping_sites:
        msg = f"{context} contain no overlap between substrate_map and phospho_matrix"
        msg = (
            f"{msg} at the substrate-map alignment seam "
            f"(shared=0, substrate_map sites={len(unique_substrate_sites)}, "
            f"phospho_matrix rows={len(validated_matrix.index)}). "
            "Use the same phosphosite IDs on both inputs."
        )
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
                "required for scoring and prediction "
                f"(sequence IDs={len(sequence_index)}, phospho_matrix rows="
                f"{len(validated_matrix.index)})."
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
    validated_pred_mat = PredMatForSignalomeSchema.validate(
        pred_mat,
        context=pred_mat_context,
    )
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
            "must share at least one phosphosite row "
            "(shared=0 at the signalome site-alignment seam; "
            f"{scoring_context} rows={len(validated_scoring_matrix.index)}, "
            f"{pred_mat_context} rows={len(validated_pred_mat.index)}, "
            f"{expression_context} rows={len(validated_expression_matrix.index)})."
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
            "kinase column "
            "(shared=0 at the signalome kinase-alignment seam; "
            f"{scoring_context} columns={len(validated_scoring_matrix.columns)}, "
            f"{pred_mat_context} columns={len(validated_pred_mat.columns)})."
        )
        raise InputCompatibilityError(msg)

    missing_koi = [
        kinase for kinase in kinases_of_interest if kinase not in common_kinases
    ]
    if missing_koi:
        missing = ", ".join(missing_koi)
        available = ", ".join(common_kinases) if common_kinases else "<none>"
        msg = (
            "kinases_of_interest are not available in the aligned signalome "
            f"inputs: {missing}. Available aligned kinases: {available}"
        )
        raise InputCompatibilityError(msg)

    if module_count is not None and module_count > len(common_sites):
        msg = (
            "module_count cannot exceed the number of aligned phosphosite rows "
            f"(module_count={module_count}, aligned_rows={len(common_sites)})"
        )
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


__all__ = [
    "DEFAULT_MIN_PRED_MAT_OVERLAP",
    "DEFAULT_MIN_PRED_MAT_OVERLAP_FRACTION",
    "PredMatOverlapSummary",
    "validate_core_column_alignment",
    "validate_pred_mat_overlap",
    "validate_signalome_alignment",
    "validate_workflow_matrix_inputs",
]
