from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from numbers import Real
from typing import TYPE_CHECKING

import pandas as pd

from phospy.motifs import KinaseMotifScorer

from .errors import InputCompatibilityError, PhospyValidationError, TableSchemaError
from .normalization import normalize_identifier_series
from .tables import PredMatSchema, SiteMatrixSchema

if TYPE_CHECKING:
    from .requests import KinaseWorkflowRequest


@dataclass(frozen=True, slots=True)
class ProteinCorrectionMatchSummary:
    """Describe phosphosite-to-protein matching before correction."""

    input_rows: int
    matched_rows: int
    unmatched_rows: int
    unmatched_fraction: float
    unmatched_gene_preview: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ValidatedKinaseActivityInputs:
    pred_mat: pd.DataFrame
    phospho_matrix: pd.DataFrame


@dataclass(frozen=True, slots=True)
class ValidatedKinaseWorkflowInputs:
    request: KinaseWorkflowRequest
    phospho_matrix: pd.DataFrame
    motif_scorer: KinaseMotifScorer | None


def _validate_fraction(
    value: float,
    *,
    name: str,
) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise PhospyValidationError(
            f"{name} must be a finite numeric value between 0 and 1"
        )

    resolved = float(value)
    if (
        not 0.0 <= resolved <= 1.0
        or resolved != resolved
        or resolved in {float("inf"), float("-inf")}
    ):
        raise PhospyValidationError(
            f"{name} must be a finite numeric value between 0 and 1"
        )
    return resolved


def _require_columns(
    frame: pd.DataFrame,
    *,
    required_columns: Sequence[str],
    context: str,
) -> None:
    missing_columns = [
        column for column in required_columns if column not in frame.columns
    ]
    if missing_columns:
        joined_columns = ", ".join(missing_columns)
        raise TableSchemaError(
            f"{context} is missing required columns: {joined_columns}"
        )


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
) -> ProteinCorrectionMatchSummary:
    resolved_max_unmatched_fraction = _validate_fraction(
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

    _require_columns(
        phospho_df,
        required_columns=[phospho_gene_col, *phospho_cols],
        context=f"{context} phospho input",
    )
    _require_columns(
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


def build_kinase_activity_inputs(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    *,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
    min_overlap: int = 1,
    min_fraction: float = 0.1,
) -> ValidatedKinaseActivityInputs:
    validated_pred_mat = PredMatSchema.validate(pred_mat, context=pred_context)
    return build_loaded_kinase_activity_inputs(
        validated_pred_mat=validated_pred_mat,
        phospho_matrix=phospho_matrix,
        pred_context=pred_context,
        matrix_context=matrix_context,
        min_overlap=min_overlap,
        min_fraction=min_fraction,
    )


def build_loaded_kinase_activity_inputs(
    *,
    validated_pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
    min_overlap: int = 1,
    min_fraction: float = 0.1,
) -> ValidatedKinaseActivityInputs:
    validated_matrix = SiteMatrixSchema.validate(phospho_matrix, context=matrix_context)
    validate_pred_mat_overlap(
        validated_pred_mat,
        validated_matrix,
        pred_context=pred_context,
        matrix_context=matrix_context,
        min_overlap=min_overlap,
        min_fraction=min_fraction,
    )
    return ValidatedKinaseActivityInputs(
        pred_mat=validated_pred_mat,
        phospho_matrix=validated_matrix,
    )


def validate_kinase_activity_inputs(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    *,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
    min_overlap: int = 1,
    min_fraction: float = 0.1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    validated_inputs = build_kinase_activity_inputs(
        pred_mat=pred_mat,
        phospho_matrix=phospho_matrix,
        pred_context=pred_context,
        matrix_context=matrix_context,
        min_overlap=min_overlap,
        min_fraction=min_fraction,
    )
    return validated_inputs.pred_mat, validated_inputs.phospho_matrix


def build_workflow_request_inputs(
    request: KinaseWorkflowRequest,
    *,
    flank_size: int,
    context: str = "Kinase workflow inputs",
) -> ValidatedKinaseWorkflowInputs:
    validated_matrix = _validate_workflow_matrix_inputs(
        request.phospho_matrix,
        request.substrate_map,
        request.site_sequences,
        context=context,
    )
    motif_scorer = (
        None
        if request.motif_sequences is None
        else KinaseMotifScorer.from_substrate_sequences(
            motif_sequences=request.motif_sequences,
            flank_size=flank_size,
        )
    )
    return ValidatedKinaseWorkflowInputs(
        request=request,
        phospho_matrix=validated_matrix,
        motif_scorer=motif_scorer,
    )


def validate_workflow_request(
    request: KinaseWorkflowRequest,
    *,
    flank_size: int = 7,
    context: str = "Kinase workflow inputs",
) -> pd.DataFrame:
    return build_workflow_request_inputs(
        request,
        flank_size=flank_size,
        context=context,
    ).phospho_matrix


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
    flank_size: int = 7,
    context: str = "Kinase workflow inputs",
) -> pd.DataFrame:
    validated_matrix = _validate_workflow_matrix_inputs(
        phospho_matrix,
        substrate_map,
        site_sequences,
        context=context,
    )
    if motif_sequences is not None:
        KinaseMotifScorer.from_substrate_sequences(
            motif_sequences=motif_sequences,
            flank_size=flank_size,
        )
    return validated_matrix


def _validate_workflow_matrix_inputs(
    phospho_matrix: pd.DataFrame,
    substrate_map: Mapping[str, Sequence[str]],
    site_sequences: Mapping[str, str] | pd.Series | None,
    *,
    context: str,
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
