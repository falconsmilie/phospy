from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from numbers import Integral
from typing import Literal, overload

import pandas as pd

from ..datasets.schema import DatasetSchema
from ..errors import PhospyValidationError
from ..internal.constants import LOCALIZATION_PROB_COLUMN, ComparisonSpec
from ..validation.schema.frames import (
    require_columns,
    require_numeric_columns,
    require_numeric_series,
)
from ..validation.values.collections import resolve_required_columns
from ..validation.values.numeric import validate_fraction, validate_non_negative_int
from .primitives import (
    _add_pairwise_comparisons_in_place,
    _collapse_duplicate_genes_owned,
    _filter_localized_sites_without_copy,
    _filter_min_observed_without_copy,
    _replace_sentinel_with_nan_in_place,
)
from .protein_correction import (
    ProteinCorrectionResult,
    ProteinCorrectionSummary,
    run_protein_correction,
)

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


@dataclass(slots=True)
class LocalizationFilterResult:
    """Localisation filter result bundle with the filtered phosphosite table.

    The contained DataFrame is mutable pandas state; this wrapper does not imply
    immutable snapshot semantics.
    """

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


@dataclass(slots=True)
class CoverageFilterResult:
    """Coverage filter result bundle with the filtered phosphosite table.

    The contained DataFrame is mutable pandas state; this wrapper does not imply
    immutable snapshot semantics.
    """

    filtered: pd.DataFrame
    summary: CoverageFilterSummary


def _validate_non_negative_integer(value: int, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise PhospyValidationError(f"{name} must be a non-negative integer")

    resolved = int(value)
    try:
        validate_non_negative_int(resolved, name)
    except PhospyValidationError as exc:
        raise PhospyValidationError(f"{name} must be a non-negative integer") from exc
    return resolved


@overload
def filter_localized_sites(
    df: pd.DataFrame,
    *,
    localization_col: str = LOCALIZATION_PROB_COLUMN,
    threshold: float = 0.75,
    return_summary: Literal[False] = False,
) -> pd.DataFrame: ...


@overload
def filter_localized_sites(
    df: pd.DataFrame,
    *,
    localization_col: str = LOCALIZATION_PROB_COLUMN,
    threshold: float = 0.75,
    return_summary: Literal[True],
) -> LocalizationFilterResult: ...


def filter_localized_sites(
    df: pd.DataFrame,
    *,
    localization_col: str = LOCALIZATION_PROB_COLUMN,
    threshold: float = 0.75,
    return_summary: bool = False,
) -> pd.DataFrame | LocalizationFilterResult:
    """Filter phosphosites by localisation probability.

    The helper keeps rows whose localisation probability is greater than or
    equal to ``threshold`` and returns a copy of the retained rows.
    """

    require_columns(
        df,
        required_columns=[localization_col],
        context="filter_localized_sites() input",
    )
    resolved_threshold = validate_fraction(
        threshold,
        name="threshold",
    )
    localization_values = require_numeric_series(
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

    resolved_columns = resolve_required_columns(
        columns,
        argument_name="columns",
        context="filter_sites_by_coverage()",
    )
    require_columns(
        df,
        required_columns=resolved_columns,
        context="filter_sites_by_coverage() input",
    )
    resolved_min_coverage = validate_fraction(
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


def replace_sentinel_with_nan(
    df: pd.DataFrame,
    columns: Iterable[str],
    sentinel: float | int,
) -> pd.DataFrame:
    resolved_columns = resolve_required_columns(
        columns,
        argument_name="columns",
        context="replace_sentinel_with_nan()",
    )
    require_columns(
        df,
        required_columns=resolved_columns,
        context="replace_sentinel_with_nan() input",
    )
    result = df.copy()
    require_numeric_columns(
        result,
        columns=resolved_columns,
        context="replace_sentinel_with_nan()",
    )
    return _replace_sentinel_with_nan_in_place(
        result,
        resolved_columns,
        sentinel,
    )


def filter_min_observed(
    df: pd.DataFrame,
    columns: Sequence[str],
    min_observed: int,
) -> pd.DataFrame:
    resolved_columns = resolve_required_columns(
        columns,
        argument_name="columns",
        context="filter_min_observed()",
    )
    require_columns(
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
    return run_protein_correction(
        df_phospho=df_phospho,
        df_total=df_total,
        phospho_gene_col=phospho_gene_col,
        total_gene_col=total_gene_col,
        phospho_cols=phospho_cols,
        protein_cols=protein_cols,
        corrected_cols=corrected_cols,
        output_prefix=output_prefix,
        max_unmatched_fraction=max_unmatched_fraction,
        return_summary=return_summary,
    )


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
