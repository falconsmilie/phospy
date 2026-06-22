"""Analysis-ready site-sequence boundary validation for dataset builder."""

from __future__ import annotations

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    PreprocessingStageExecution,
)
from phospy.science.sites.identifiers import (
    canonicalize_site_components,
    canonicalize_site_identifier,
)

_SITE_SEQUENCE_RESOLUTION_STAGE = DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION
_SITE_SEQUENCE_FAILURE_CATEGORY_MISSING_REFERENCE_SUPPORT = "missing_reference_support"
_SITE_SEQUENCE_FAILURE_CATEGORY_AMBIGUOUS_MAPPING = "ambiguous_mapping"
_SITE_SEQUENCE_FAILURE_CATEGORY_INVALID_METADATA = "invalid_metadata"


class AnalysisReadySiteSequenceValidator:
    """Enforce `site_metadata.site_sequence` validity before dataset construction."""

    def run(
        self,
        *,
        site_metadata: pd.DataFrame,
        preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    ) -> pd.DataFrame:
        if "site_sequence" not in site_metadata.columns:
            _raise_site_sequence_boundary_error(
                site_metadata=site_metadata,
                invalid_index=list(range(int(len(site_metadata.index)))),
                preprocessing_trace=preprocessing_trace,
                reason_override="site_sequence column is missing after preprocessing",
            )
        sequence_column = site_metadata.loc[:, "site_sequence"]
        missing_or_blank_mask = sequence_column.isna()
        missing_or_blank_mask = missing_or_blank_mask | (
            sequence_column.astype("string").str.strip().isna()
        )
        missing_or_blank_mask = missing_or_blank_mask | (
            sequence_column.astype("string").str.strip() == ""
        )
        non_string_mask = ~sequence_column.map(lambda value: isinstance(value, str))
        non_string_mask = non_string_mask & ~sequence_column.isna()
        invalid_mask = missing_or_blank_mask | non_string_mask
        if not bool(invalid_mask.any()):
            return site_metadata
        invalid_positions = [
            position
            for position, flagged in enumerate(invalid_mask.tolist())
            if bool(flagged)
        ]
        _raise_site_sequence_boundary_error(
            site_metadata=site_metadata,
            invalid_index=invalid_positions,
            preprocessing_trace=preprocessing_trace,
            reason_override=None,
        )
        return site_metadata


def _raise_site_sequence_boundary_error(
    *,
    site_metadata: pd.DataFrame,
    invalid_index: list[int],
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    reason_override: str | None,
) -> None:
    context_by_row_id = _resolve_site_sequence_failure_context(preprocessing_trace)
    row_messages: list[str] = []
    for position in invalid_index[:8]:
        row_id = str(site_metadata.index[position])
        row = site_metadata.iloc[position]
        gene_symbol = _optional_row_value(row, "gene_symbol")
        site_value = _optional_row_value(row, "site")
        protein_id = _optional_row_value(row, "protein_id")
        category, identity_reason = _classify_identity_failure(
            row_id=row_id,
            gene_symbol=gene_symbol,
            site=site_value,
        )
        context = context_by_row_id.get(row_id)
        if context is None:
            reason = (
                identity_reason
                if reason_override is None
                else f"{reason_override}; {identity_reason}"
            )
            reason_category = category
        else:
            context_reason, context_category = context
            if reason_override is None:
                reason = context_reason
            else:
                reason = f"{reason_override}; {context_reason}"
            reason_category = context_category
        row_messages.append(
            f"site_id={row_id!r}, gene_symbol={gene_symbol!r}, site={site_value!r}, "
            f"protein_id={protein_id!r}, reason={reason!r}, failure_category={reason_category!r}"
        )
    suffix = "" if len(invalid_index) <= 8 else f"; +{len(invalid_index) - 8} more rows"
    details = " | ".join(row_messages)
    raise PhosPyInputError(
        "dataset builder cannot construct AnalysisReadyPhosphoDataset because "
        "dataset.site_metadata.site_sequence is missing, blank, or invalid after "
        f"builder enrichment; unresolved_rows={len(invalid_index)}; row_context={details}{suffix}"
    )


def _resolve_site_sequence_failure_context(
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
) -> dict[str, tuple[str, str]]:
    if preprocessing_trace is None:
        return {}
    contexts: dict[str, tuple[str, str]] = {}
    for stage in preprocessing_trace:
        if stage.stage != _SITE_SEQUENCE_RESOLUTION_STAGE:
            continue
        diagnostics = stage.diagnostics if stage.diagnostics is not None else {}
        row_diagnostics = diagnostics.get("row_diagnostics")
        if not isinstance(row_diagnostics, list):
            continue
        for row in row_diagnostics:
            if not isinstance(row, dict):
                continue
            row_id_value = row.get("row_id")
            if row_id_value is None:
                continue
            row_id = str(row_id_value)
            resolved_value = row.get("resolved_site_sequence")
            if (
                resolved_value is not None
                and str(resolved_value).strip() != ""
                and str(resolved_value).strip().lower() != "none"
            ):
                continue
            status = str(row.get("status", "")).strip().lower()
            reason = str(row.get("reason", "")).strip()
            category = _classify_stage_failure_category(status, reason)
            resolved_reason = (
                reason
                if reason
                else "site_sequence resolution failed during FASTA enrichment"
            )
            contexts[row_id] = (resolved_reason, category)
    return contexts


def _classify_stage_failure_category(status: str, reason: str) -> str:
    combined = f"{status} {reason}".strip().lower()
    if "ambiguous" in combined:
        return _SITE_SEQUENCE_FAILURE_CATEGORY_AMBIGUOUS_MAPPING
    invalid_tokens = (
        "missing_accession",
        "missing_existing_sequence",
        "invalid",
        "blank",
        "empty",
    )
    if any(token in combined for token in invalid_tokens):
        return _SITE_SEQUENCE_FAILURE_CATEGORY_INVALID_METADATA
    return _SITE_SEQUENCE_FAILURE_CATEGORY_MISSING_REFERENCE_SUPPORT


def _classify_identity_failure(
    *,
    row_id: str,
    gene_symbol: str | None,
    site: str | None,
) -> tuple[str, str]:
    metadata_site_id: str | None = None
    index_site_id: str | None = None
    metadata_error: str | None = None
    index_error: str | None = None
    if gene_symbol is not None and site is not None:
        try:
            metadata_site_id = canonicalize_site_components(
                gene_symbol,
                site,
                field_name=(
                    "dataset build request site_metadata.gene_symbol/site for "
                    "analysis-ready sequence enforcement"
                ),
                error_type=PhosPyInputError,
            )
        except PhosPyInputError as exc:
            metadata_error = str(exc)
    else:
        metadata_error = "missing gene_symbol/site metadata"
    try:
        if ";" in str(row_id):
            index_site_id = canonicalize_site_identifier(
                row_id,
                field_name=(
                    "dataset build request site_metadata.index for analysis-ready "
                    "sequence enforcement"
                ),
                error_type=PhosPyInputError,
            )
    except PhosPyInputError as exc:
        index_error = str(exc)
    if metadata_site_id is not None and index_site_id is not None:
        if metadata_site_id != index_site_id:
            return (
                _SITE_SEQUENCE_FAILURE_CATEGORY_AMBIGUOUS_MAPPING,
                "gene/site metadata and site index map to different site identities",
            )
        return (
            _SITE_SEQUENCE_FAILURE_CATEGORY_MISSING_REFERENCE_SUPPORT,
            "site_sequence is unresolved for the expected site identity despite "
            "coherent metadata",
        )
    if metadata_site_id is not None or index_site_id is not None:
        return (
            _SITE_SEQUENCE_FAILURE_CATEGORY_MISSING_REFERENCE_SUPPORT,
            "site_sequence is unresolved for the expected site identity despite "
            "partial identity metadata",
        )
    detail_chunks = [chunk for chunk in (metadata_error, index_error) if chunk]
    detail = (
        "; ".join(detail_chunks) if detail_chunks else "invalid site identity metadata"
    )
    return (_SITE_SEQUENCE_FAILURE_CATEGORY_INVALID_METADATA, detail)


def _optional_row_value(row: pd.Series, column_name: str) -> str | None:
    if column_name not in row.index:
        return None
    value = row[column_name]
    if pd.isna(value):
        return None
    if isinstance(value, str):
        stripped = value.strip()
        return stripped if stripped else None
    return str(value)
