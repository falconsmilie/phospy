from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from numbers import Integral, Real
from typing import Literal, overload

import numpy as np
import pandas as pd

from .constants import ComparisonSpec
from .dataset_schema import DatasetSchema
from .validation.compatibility import (
    validate_core_column_alignment,
    validate_protein_correction_inputs,
)
from .validation.errors import (
    InputCompatibilityError,
    PhospyValidationError,
    TableSchemaError,
)
from .validation.normalization import normalize_identifier_series
from .validation.primitives import validate_non_negative_int

"""Standalone preprocessing helpers.

`PhosphoDataset.preprocessing.run()` is the preferred public entrypoint for the
end-to-end core preprocessing path. This module contains targeted helpers for
advanced use when a caller only wants a specific preprocessing step.
"""


@dataclass(frozen=True, slots=True)
class LocalizationFilterSummary:
    """Small summary describing a phosphosite localisation filter pass."""

    input_rows: int
    retained_rows: int
    removed_rows: int
    retention_fraction: float
    threshold: float
    localization_col: str


@dataclass(frozen=True, slots=True)
class LocalizationFilterResult:
    """Filtered phosphosites together with a localisation filter summary."""

    filtered: pd.DataFrame
    summary: LocalizationFilterSummary


@dataclass(frozen=True, slots=True)
class CoverageFilterSummary:
    """Small summary describing a phosphosite coverage filter pass."""

    input_rows: int
    retained_rows: int
    removed_rows: int
    retention_fraction: float
    min_coverage: float
    required_observed_count: int
    value_columns: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class CoverageFilterResult:
    """Filtered phosphosites together with a coverage filter summary."""

    filtered: pd.DataFrame
    summary: CoverageFilterSummary


@dataclass(frozen=True, slots=True)
class ProteinCorrectionSummary:
    """Describe a phospho-to-protein correction pass."""

    input_rows: int
    matched_rows: int
    unmatched_rows: int
    unmatched_fraction: float
    phospho_gene_col: str
    total_gene_col: str
    unmatched_gene_preview: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ProteinCorrectionResult:
    """Corrected phosphosite rows together with protein-match metadata."""

    corrected: pd.DataFrame
    summary: ProteinCorrectionSummary


def _require_columns(
    df: pd.DataFrame,
    *,
    required_columns: Sequence[str],
    context: str,
) -> None:
    missing_columns = [
        column for column in required_columns if column not in df.columns
    ]
    if missing_columns:
        joined_columns = ", ".join(missing_columns)
        msg = f"{context} is missing required columns: {joined_columns}"
        raise TableSchemaError(msg)


def _validate_probability_threshold(value: float, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise PhospyValidationError(
            f"{name} must be a finite numeric value between 0 and 1"
        )

    resolved = float(value)
    if not math.isfinite(resolved) or not 0.0 <= resolved <= 1.0:
        raise PhospyValidationError(
            f"{name} must be a finite numeric value between 0 and 1"
        )

    return resolved


def _require_numeric_series(
    values: pd.Series,
    *,
    column: str,
    context: str,
) -> pd.Series:
    try:
        numeric_values = pd.to_numeric(values, errors="raise")
    except (TypeError, ValueError) as exc:
        raise TableSchemaError(
            f"{context} requires numeric values in column '{column}'"
        ) from exc

    return pd.Series(numeric_values, index=values.index, copy=False)


def _resolve_required_columns(
    columns: Iterable[str],
    *,
    argument_name: str,
    context: str,
) -> list[str]:
    resolved_columns = list(columns)
    if not resolved_columns:
        raise PhospyValidationError(
            f"{context} requires at least one column name in '{argument_name}'"
        )
    return resolved_columns


def _validate_non_negative_integer(value: int, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise PhospyValidationError(f"{name} must be a non-negative integer")

    resolved = int(value)
    try:
        validate_non_negative_int(resolved, name)
    except PhospyValidationError as exc:
        raise PhospyValidationError(f"{name} must be a non-negative integer") from exc
    return resolved


def _require_numeric_columns(
    df: pd.DataFrame,
    *,
    columns: Sequence[str],
    context: str,
) -> None:
    for column in columns:
        df[column] = _require_numeric_series(
            df[column],
            column=column,
            context=context,
        )


@overload
def filter_localized_sites(
    df: pd.DataFrame,
    *,
    localization_col: str = "localization_prob",
    threshold: float = 0.75,
    return_summary: Literal[False] = False,
) -> pd.DataFrame: ...


@overload
def filter_localized_sites(
    df: pd.DataFrame,
    *,
    localization_col: str = "localization_prob",
    threshold: float = 0.75,
    return_summary: Literal[True],
) -> LocalizationFilterResult: ...


def filter_localized_sites(
    df: pd.DataFrame,
    *,
    localization_col: str = "localization_prob",
    threshold: float = 0.75,
    return_summary: bool = False,
) -> pd.DataFrame | LocalizationFilterResult:
    """Filter phosphosites by localisation probability.

    The helper keeps rows whose localisation probability is greater than or
    equal to ``threshold`` and returns a copy of the retained rows.
    """

    _require_columns(
        df,
        required_columns=[localization_col],
        context="filter_localized_sites() input",
    )
    resolved_threshold = _validate_probability_threshold(
        threshold,
        name="threshold",
    )
    localization_values = _require_numeric_series(
        df[localization_col],
        column=localization_col,
        context="filter_localized_sites()",
    )

    filtered = _filter_localized_sites_without_copy(
        df,
        localization_col=localization_col,
        threshold=resolved_threshold,
        localization_values=localization_values,
    ).copy()
    if not return_summary:
        return filtered

    input_rows = int(len(df))
    retained_rows = int(len(filtered))
    removed_rows = input_rows - retained_rows
    summary = LocalizationFilterSummary(
        input_rows=input_rows,
        retained_rows=retained_rows,
        removed_rows=removed_rows,
        retention_fraction=(retained_rows / input_rows) if input_rows else 0.0,
        threshold=resolved_threshold,
        localization_col=localization_col,
    )
    return LocalizationFilterResult(filtered=filtered, summary=summary)


def _filter_localized_sites_without_copy(
    df: pd.DataFrame,
    *,
    localization_col: str,
    threshold: float,
    localization_values: pd.Series | None = None,
) -> pd.DataFrame:
    values = localization_values
    if values is None:
        values = _require_numeric_series(
            df[localization_col],
            column=localization_col,
            context="filter_localized_sites()",
        )
    return df.loc[values >= threshold]


@overload
def filter_sites_by_coverage(
    df: pd.DataFrame,
    *,
    columns: Sequence[str],
    min_coverage: float = 0.0,
    return_summary: Literal[False] = False,
) -> pd.DataFrame: ...


@overload
def filter_sites_by_coverage(
    df: pd.DataFrame,
    *,
    columns: Sequence[str],
    min_coverage: float = 0.0,
    return_summary: Literal[True],
) -> CoverageFilterResult: ...


def filter_sites_by_coverage(
    df: pd.DataFrame,
    *,
    columns: Sequence[str],
    min_coverage: float = 0.0,
    return_summary: bool = False,
) -> pd.DataFrame | CoverageFilterResult:
    """Filter phosphosites by the minimum observed proportion across samples.

    The helper keeps rows whose observed-value proportion across ``columns`` is
    greater than or equal to ``min_coverage``. Coverage is evaluated as the
    fraction of non-missing values among the selected sample columns.
    """

    resolved_columns = _resolve_required_columns(
        columns,
        argument_name="columns",
        context="filter_sites_by_coverage()",
    )
    _require_columns(
        df,
        required_columns=resolved_columns,
        context="filter_sites_by_coverage() input",
    )
    resolved_min_coverage = _validate_probability_threshold(
        min_coverage,
        name="min_coverage",
    )
    required_observed_count = math.ceil(resolved_min_coverage * len(resolved_columns))

    filtered = filter_min_observed(
        df,
        resolved_columns,
        min_observed=required_observed_count,
    )
    if not return_summary:
        return filtered

    input_rows = int(len(df))
    retained_rows = int(len(filtered))
    removed_rows = input_rows - retained_rows
    summary = CoverageFilterSummary(
        input_rows=input_rows,
        retained_rows=retained_rows,
        removed_rows=removed_rows,
        retention_fraction=(retained_rows / input_rows) if input_rows else 0.0,
        min_coverage=resolved_min_coverage,
        required_observed_count=required_observed_count,
        value_columns=tuple(resolved_columns),
    )
    return CoverageFilterResult(filtered=filtered, summary=summary)


def _replace_sentinel_with_nan_in_place(
    df: pd.DataFrame,
    columns: Iterable[str],
    sentinel: float | int,
) -> pd.DataFrame:
    cols = list(columns)
    for col in cols:
        df[col] = df[col].astype(float).replace(sentinel, np.nan)
    return df


def replace_sentinel_with_nan(
    df: pd.DataFrame,
    columns: Iterable[str],
    sentinel: float | int,
) -> pd.DataFrame:
    resolved_columns = _resolve_required_columns(
        columns,
        argument_name="columns",
        context="replace_sentinel_with_nan()",
    )
    _require_columns(
        df,
        required_columns=resolved_columns,
        context="replace_sentinel_with_nan() input",
    )
    result = df.copy()
    _require_numeric_columns(
        result,
        columns=resolved_columns,
        context="replace_sentinel_with_nan()",
    )
    return _replace_sentinel_with_nan_in_place(
        result,
        resolved_columns,
        sentinel,
    )


def _filter_min_observed_without_copy(
    df: pd.DataFrame,
    columns: Sequence[str],
    min_observed: int,
) -> pd.DataFrame:
    mask = df.loc[:, list(columns)].notna().sum(axis=1) >= min_observed
    return df.loc[mask]


def filter_min_observed(
    df: pd.DataFrame,
    columns: Sequence[str],
    min_observed: int,
) -> pd.DataFrame:
    resolved_columns = _resolve_required_columns(
        columns,
        argument_name="columns",
        context="filter_min_observed()",
    )
    _require_columns(
        df,
        required_columns=resolved_columns,
        context="filter_min_observed() input",
    )
    resolved_min_observed = _validate_non_negative_integer(
        min_observed,
        name="min_observed",
    )
    return _filter_min_observed_without_copy(
        df,
        resolved_columns,
        resolved_min_observed,
    ).copy()


def _collapse_duplicate_genes_owned(
    df: pd.DataFrame,
    gene_col: str,
    value_cols: Sequence[str],
    uppercase: bool = True,
) -> pd.DataFrame:
    """Collapse duplicate gene rows using an explicit ranking policy.

    Gene identifiers are converted to pandas ``string`` values before ranking.
    When ``uppercase`` is ``True``, identifiers are normalised to uppercase
    before duplicate grouping so mixed-case duplicates collapse together.

    Rows are ranked within each gene group by:
    1. highest observed-value count across ``value_cols``
    2. highest mean signal across ``value_cols``
    3. earliest original row order as a stable tie-breaker

    Rows with zero observed values across ``value_cols`` are dropped before the
    winning row is selected so all-missing groups do not survive deduplication.
    The top-ranked remaining row for each gene is retained.
    """
    _require_columns(
        df,
        required_columns=[gene_col, *value_cols],
        context="collapse_duplicate_genes() input",
    )

    df[gene_col] = df[gene_col].astype("string")
    if uppercase:
        df[gene_col] = df[gene_col].str.upper()
    ranked_cols = list(value_cols)
    df["__observed_count"] = df.loc[:, ranked_cols].notna().sum(axis=1)
    df["__mean_signal"] = df.loc[:, ranked_cols].mean(axis=1, skipna=True)
    df["__original_order"] = np.arange(len(df), dtype=int)

    ranked = df.sort_values(
        by=[gene_col, "__observed_count", "__mean_signal", "__original_order"],
        ascending=[True, False, False, True],
        kind="mergesort",
    )
    ranked = ranked.loc[ranked["__observed_count"] > 0]
    result = ranked.drop_duplicates(subset=[gene_col], keep="first").drop(
        columns=["__observed_count", "__mean_signal", "__original_order"]
    )

    return result.reset_index(drop=True)


def collapse_duplicate_genes(
    df: pd.DataFrame,
    gene_col: str,
    value_cols: Sequence[str],
    uppercase: bool = True,
) -> pd.DataFrame:
    work = df.copy()
    return _collapse_duplicate_genes_owned(
        work,
        gene_col=gene_col,
        value_cols=value_cols,
        uppercase=uppercase,
    )


@overload
def correct_phospho_to_protein(
    df_phospho: pd.DataFrame,
    df_total: pd.DataFrame,
    phospho_gene_col: str,
    total_gene_col: str,
    phospho_cols: Sequence[str],
    protein_cols: Sequence[str],
    corrected_cols: Sequence[str] | None = None,
    output_prefix: str = "phospho_corrected_",
    *,
    max_unmatched_fraction: float = 1.0,
    return_summary: Literal[False] = False,
) -> pd.DataFrame: ...


@overload
def correct_phospho_to_protein(
    df_phospho: pd.DataFrame,
    df_total: pd.DataFrame,
    phospho_gene_col: str,
    total_gene_col: str,
    phospho_cols: Sequence[str],
    protein_cols: Sequence[str],
    corrected_cols: Sequence[str] | None = None,
    output_prefix: str = "phospho_corrected_",
    *,
    max_unmatched_fraction: float = 1.0,
    return_summary: Literal[True],
) -> ProteinCorrectionResult: ...


def correct_phospho_to_protein(
    df_phospho: pd.DataFrame,
    df_total: pd.DataFrame,
    phospho_gene_col: str,
    total_gene_col: str,
    phospho_cols: Sequence[str],
    protein_cols: Sequence[str],
    corrected_cols: Sequence[str] | None = None,
    output_prefix: str = "phospho_corrected_",
    *,
    max_unmatched_fraction: float = 1.0,
    return_summary: bool = False,
) -> pd.DataFrame | ProteinCorrectionResult:
    resolved_phospho_cols = _resolve_required_columns(
        phospho_cols,
        argument_name="phospho_cols",
        context="correct_phospho_to_protein()",
    )
    resolved_protein_cols = _resolve_required_columns(
        protein_cols,
        argument_name="protein_cols",
        context="correct_phospho_to_protein()",
    )
    resolved_corrected_cols = (
        list(corrected_cols)
        if corrected_cols is not None
        else [
            f"{output_prefix}{idx}" for idx in range(1, len(resolved_phospho_cols) + 1)
        ]
    )

    validate_core_column_alignment(
        resolved_protein_cols,
        resolved_phospho_cols,
        resolved_corrected_cols,
        context="correct_phospho_to_protein() inputs",
    )
    _require_columns(
        df_phospho,
        required_columns=[phospho_gene_col, *resolved_phospho_cols],
        context="correct_phospho_to_protein() phospho input",
    )
    _require_columns(
        df_total,
        required_columns=[total_gene_col, *resolved_protein_cols],
        context="correct_phospho_to_protein() total input",
    )

    phospho_join_col = "__phospy_normalized_phospho_gene_key"
    total_join_col = "__phospy_normalized_total_gene_key"

    phospho_work = df_phospho.copy()
    total_work = df_total.copy()
    _require_numeric_columns(
        phospho_work,
        columns=resolved_phospho_cols,
        context="correct_phospho_to_protein() phospho input",
    )
    _require_numeric_columns(
        total_work,
        columns=resolved_protein_cols,
        context="correct_phospho_to_protein() total input",
    )
    match_summary = validate_protein_correction_inputs(
        phospho_work,
        total_work,
        phospho_gene_col=phospho_gene_col,
        total_gene_col=total_gene_col,
        phospho_cols=resolved_phospho_cols,
        protein_cols=resolved_protein_cols,
        max_unmatched_fraction=max_unmatched_fraction,
        context="correct_phospho_to_protein() inputs",
    )

    phospho_work[phospho_join_col] = normalize_identifier_series(
        phospho_work[phospho_gene_col]
    )
    total_work[total_join_col] = normalize_identifier_series(total_work[total_gene_col])

    if total_work[total_join_col].duplicated().any():
        msg = (
            f"{total_gene_col} must be unique before protein correction to avoid "
            "duplicating phosphosite rows during the merge"
        )
        raise InputCompatibilityError(msg)

    merged = phospho_work.merge(
        total_work[[total_join_col, total_gene_col, *resolved_protein_cols]],
        left_on=phospho_join_col,
        right_on=total_join_col,
        how="inner",
    )

    drop_columns: list[str] = [phospho_join_col, total_join_col]
    if total_gene_col != phospho_gene_col and total_gene_col in merged.columns:
        drop_columns.append(total_gene_col)
    merged = merged.drop(columns=drop_columns, errors="ignore")

    for corrected_col, p_col, t_col in zip(
        resolved_corrected_cols,
        resolved_phospho_cols,
        resolved_protein_cols,
        strict=True,
    ):
        merged[corrected_col] = merged[p_col] - merged[t_col]

    if not return_summary:
        return merged

    summary = ProteinCorrectionSummary(
        input_rows=match_summary.input_rows,
        matched_rows=match_summary.matched_rows,
        unmatched_rows=match_summary.unmatched_rows,
        unmatched_fraction=match_summary.unmatched_fraction,
        phospho_gene_col=phospho_gene_col,
        total_gene_col=total_gene_col,
        unmatched_gene_preview=match_summary.unmatched_gene_preview,
    )
    return ProteinCorrectionResult(corrected=merged, summary=summary)


def _add_pairwise_comparisons_in_place(
    df: pd.DataFrame,
    comparisons: Sequence[ComparisonSpec],
    group_to_corrected_col: dict[str, str] | None = None,
    output_prefix: str = "p_",
    schema: DatasetSchema | None = None,
) -> pd.DataFrame:
    if group_to_corrected_col is None:
        resolved_schema = schema or DatasetSchema()
        resolved_schema.validate_comparisons(
            comparisons,
            context="Pairwise comparison configuration",
        )
        group_to_corrected_col = dict(resolved_schema.group_to_corrected_col)
    else:
        resolved = tuple(comparisons)
        valid_groups = frozenset(group_to_corrected_col)
        seen: set[ComparisonSpec] = set()
        for left, right in resolved:
            if left not in valid_groups or right not in valid_groups:
                raise KeyError(f"Missing group mapping for comparison: {(left, right)}")
            pair = (left, right)
            if pair in seen:
                raise ValueError(f"Duplicate comparison pair: {left!r}, {right!r}")
            seen.add(pair)

    for left, right in comparisons:
        if left not in group_to_corrected_col or right not in group_to_corrected_col:
            raise KeyError(f"Missing group mapping for comparison: {(left, right)}")
        df[f"{output_prefix}{left}_{right}"] = (
            df[group_to_corrected_col[left]] - df[group_to_corrected_col[right]]
        )

    return df


def add_pairwise_comparisons(
    df: pd.DataFrame,
    comparisons: Sequence[ComparisonSpec],
    group_to_corrected_col: dict[str, str] | None = None,
    output_prefix: str = "p_",
    schema: DatasetSchema | None = None,
) -> pd.DataFrame:
    result = df.copy()
    return _add_pairwise_comparisons_in_place(
        result,
        comparisons,
        group_to_corrected_col=group_to_corrected_col,
        output_prefix=output_prefix,
        schema=schema,
    )


__all__ = [
    "CoverageFilterResult",
    "CoverageFilterSummary",
    "LocalizationFilterResult",
    "LocalizationFilterSummary",
    "ProteinCorrectionResult",
    "ProteinCorrectionSummary",
    "add_pairwise_comparisons",
    "collapse_duplicate_genes",
    "correct_phospho_to_protein",
    "filter_localized_sites",
    "filter_min_observed",
    "filter_sites_by_coverage",
    "replace_sentinel_with_nan",
]
