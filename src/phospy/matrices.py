from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd

from .errors import TableSchemaError
from .internal.constants import (
    PHOSPHO_GENE_COLUMN,
    PHOSPHO_UID_COLUMN,
    SITE_MATRIX_GENE_COLUMN,
    SITE_MATRIX_ID_COLUMN,
    SITE_MATRIX_P_SITE_COLUMN,
)
from .internal.types import DuplicateSiteStrategy
from .validation.schema.tables import SiteMatrixSourceSchema
from .validation.values.enums import validate_duplicate_site_strategy

__all__ = [
    "DEFAULT_SITE_MATRIX_POLICY",
    "SiteMatrixPolicy",
    "build_site_matrix",
    "format_row_drop_diagnostics",
]


@dataclass(frozen=True, slots=True)
class SiteMatrixPolicy:
    """Explicit policy for duplicate phosphosite rows during site-matrix building.

    ``duplicate_site_strategy`` controls how multiple corrected phosphosite rows
    that resolve to the same site identifier are collapsed:

    - ``"max_mean_signal"`` keeps the row with the largest mean corrected signal
    - ``"first"`` keeps the first encountered complete row
    - ``"aggregate_mean"`` averages corrected values across duplicate rows
    - ``"aggregate_median"`` takes the column-wise median across duplicate rows
    - ``"error"`` rejects duplicate site IDs instead of collapsing them
    """

    duplicate_site_strategy: DuplicateSiteStrategy = "max_mean_signal"

    def __post_init__(self) -> None:
        validate_duplicate_site_strategy(self.duplicate_site_strategy)

    @classmethod
    def from_value(cls, value: object) -> SiteMatrixPolicy:
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if isinstance(value, Mapping):
            return cls(**dict(value))
        msg = "site_matrix_policy must be a SiteMatrixPolicy or mapping"
        raise TypeError(msg)


DEFAULT_SITE_MATRIX_POLICY = SiteMatrixPolicy()


def format_row_drop_diagnostics(row_drop_stats: Mapping[str, int | str]) -> str:
    """Format human-readable site-matrix row-drop diagnostics."""

    stats = {
        "input_rows": int(row_drop_stats.get("input_rows", 0)),
        "dropped_missing_sequence": int(
            row_drop_stats.get("dropped_missing_sequence", 0)
        ),
        "dropped_incomplete_values": int(
            row_drop_stats.get("dropped_incomplete_values", 0)
        ),
        "deduplicated_site_rows": int(row_drop_stats.get("deduplicated_site_rows", 0)),
        "retained_rows": int(row_drop_stats.get("retained_rows", 0)),
        "duplicate_site_strategy": str(
            row_drop_stats.get(
                "duplicate_site_strategy",
                DEFAULT_SITE_MATRIX_POLICY.duplicate_site_strategy,
            )
        ),
    }
    known_drops = (
        stats["dropped_missing_sequence"]
        + stats["dropped_incomplete_values"]
        + stats["deduplicated_site_rows"]
    )
    other_dropped_rows = max(
        stats["input_rows"] - stats["retained_rows"] - known_drops, 0
    )
    return (
        "row-drop diagnostics: "
        f"input_rows={stats['input_rows']}, "
        f"dropped_missing_sequence={stats['dropped_missing_sequence']}, "
        f"dropped_incomplete_values={stats['dropped_incomplete_values']}, "
        f"deduplicated_site_rows={stats['deduplicated_site_rows']}, "
        f"duplicate_site_strategy={stats['duplicate_site_strategy']}, "
        f"other_dropped_rows={other_dropped_rows}, "
        f"retained_rows={stats['retained_rows']}"
    )


def build_site_matrix(
    df: pd.DataFrame,
    gene_p_site_col: str,
    sequence_col: str,
    value_cols: Sequence[str],
    gene_col_name: str = SITE_MATRIX_GENE_COLUMN,
    p_site_col_name: str = SITE_MATRIX_P_SITE_COLUMN,
    *,
    policy: SiteMatrixPolicy | None = None,
    duplicate_site_strategy: DuplicateSiteStrategy | None = None,
    copy_frame: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    resolved_policy = _resolve_site_matrix_policy(
        policy=policy,
        duplicate_site_strategy=duplicate_site_strategy,
    )
    work = SiteMatrixSourceSchema.validate(
        df,
        gene_p_site_col=gene_p_site_col,
        sequence_col=sequence_col,
        value_cols=value_cols,
        context="site-matrix source table",
        copy_frame=copy_frame,
    )

    source_cols = [
        col
        for col in [
            PHOSPHO_GENE_COLUMN,
            PHOSPHO_UID_COLUMN,
            gene_p_site_col,
            sequence_col,
            *list(value_cols),
        ]
        if col in work.columns
    ]
    site_work = work.loc[:, source_cols].copy(deep=True)

    split_cols = (
        site_work[gene_p_site_col]
        .astype("string")
        .str.split(
            "_",
            n=1,
            expand=True,
        )
    )
    site_work[gene_col_name] = split_cols[0].astype("string")
    site_work[p_site_col_name] = split_cols[1].astype("string")
    site_work[SITE_MATRIX_ID_COLUMN] = (
        site_work[gene_col_name].str.upper()
        + ";"
        + site_work[p_site_col_name].str.upper()
        + ";"
    )

    base_cols = [
        col
        for col in [
            PHOSPHO_GENE_COLUMN,
            gene_col_name,
            p_site_col_name,
            PHOSPHO_UID_COLUMN,
            sequence_col,
            SITE_MATRIX_ID_COLUMN,
        ]
        if col in site_work.columns
    ]
    keep_cols: list[str] = [*base_cols, *list(value_cols)]
    base_input = site_work.loc[:, keep_cols]

    total_rows = len(base_input)
    with_sequence = base_input.loc[base_input[sequence_col].notna()]
    dropped_missing_sequence = total_rows - len(with_sequence)

    complete_cases = with_sequence.loc[
        with_sequence.loc[:, list(value_cols)].notna().all(axis=1)
    ]
    dropped_incomplete_values = len(with_sequence) - len(complete_cases)

    phosr_input = _apply_duplicate_site_policy(
        complete_cases=complete_cases,
        value_cols=value_cols,
        policy=resolved_policy,
    )
    phosr_input = phosr_input.sort_values(
        SITE_MATRIX_ID_COLUMN,
        ascending=True,
        kind="stable",
    ).reset_index(drop=True)
    deduplicated_site_rows = len(complete_cases) - len(phosr_input)

    row_drop_stats: dict[str, int | str] = {
        "input_rows": total_rows,
        "dropped_missing_sequence": dropped_missing_sequence,
        "dropped_incomplete_values": dropped_incomplete_values,
        "deduplicated_site_rows": deduplicated_site_rows,
        "duplicate_site_strategy": resolved_policy.duplicate_site_strategy,
        "retained_rows": len(phosr_input),
    }

    matrix = phosr_input.loc[:, [SITE_MATRIX_ID_COLUMN, *value_cols]].set_index(
        SITE_MATRIX_ID_COLUMN
    )
    matrix.attrs["row_drop_stats"] = row_drop_stats.copy()
    sequences = phosr_input.set_index(SITE_MATRIX_ID_COLUMN)[sequence_col]
    sequences.attrs["row_drop_stats"] = row_drop_stats.copy()
    phosr_input.attrs["row_drop_stats"] = row_drop_stats.copy()
    return phosr_input, matrix, sequences


def _resolve_site_matrix_policy(
    *,
    policy: SiteMatrixPolicy | None,
    duplicate_site_strategy: DuplicateSiteStrategy | None,
) -> SiteMatrixPolicy:
    if policy is not None and duplicate_site_strategy is not None:
        msg = (
            "build_site_matrix() accepts either policy or duplicate_site_strategy, "
            "not both"
        )
        raise ValueError(msg)
    if policy is not None:
        return SiteMatrixPolicy.from_value(policy)
    if duplicate_site_strategy is not None:
        return SiteMatrixPolicy(duplicate_site_strategy=duplicate_site_strategy)
    return DEFAULT_SITE_MATRIX_POLICY


def _apply_duplicate_site_policy(
    *,
    complete_cases: pd.DataFrame,
    value_cols: Sequence[str],
    policy: SiteMatrixPolicy,
) -> pd.DataFrame:
    if complete_cases.empty:
        return complete_cases.copy().reset_index(drop=True)

    duplicate_mask = complete_cases.duplicated(SITE_MATRIX_ID_COLUMN, keep=False)
    if not bool(duplicate_mask.any()):
        return complete_cases.copy().reset_index(drop=True)

    if policy.duplicate_site_strategy == "error":
        duplicate_sites = (
            complete_cases.loc[duplicate_mask, SITE_MATRIX_ID_COLUMN]
            .astype(str)
            .drop_duplicates()
            .tolist()
        )
        preview = ", ".join(duplicate_sites[:3])
        msg = (
            "site-matrix source table contains duplicate phosphosite identifiers "
            f"after correction and duplicate_site_strategy='error': {preview}"
        )
        raise TableSchemaError(msg)

    if policy.duplicate_site_strategy == "max_mean_signal":
        work = complete_cases.copy()
        work["__mean_signal"] = work.loc[:, list(value_cols)].mean(axis=1, skipna=True)
        idx = work.groupby(SITE_MATRIX_ID_COLUMN, sort=False)["__mean_signal"].idxmax()
        return work.loc[idx].drop(columns="__mean_signal").copy().reset_index(drop=True)

    if policy.duplicate_site_strategy == "first":
        return (
            complete_cases.drop_duplicates(SITE_MATRIX_ID_COLUMN, keep="first")
            .copy()
            .reset_index(drop=True)
        )

    if policy.duplicate_site_strategy in {"aggregate_mean", "aggregate_median"}:
        return _aggregate_duplicate_sites(
            complete_cases=complete_cases,
            value_cols=value_cols,
            strategy=policy.duplicate_site_strategy,
        )

    msg = f"Unsupported duplicate_site_strategy: {policy.duplicate_site_strategy}"
    raise TableSchemaError(msg)


def _aggregate_duplicate_sites(
    *,
    complete_cases: pd.DataFrame,
    value_cols: Sequence[str],
    strategy: DuplicateSiteStrategy,
) -> pd.DataFrame:
    value_col_list = list(value_cols)
    ordered_columns = list(complete_cases.columns)
    metadata_columns = [
        column
        for column in ordered_columns
        if column not in {SITE_MATRIX_ID_COLUMN, *value_col_list}
    ]

    metadata = (
        complete_cases.loc[:, [SITE_MATRIX_ID_COLUMN, *metadata_columns]]
        .groupby(SITE_MATRIX_ID_COLUMN, sort=False)
        .first()
    )
    grouped_values = complete_cases.groupby(SITE_MATRIX_ID_COLUMN, sort=False)[
        value_col_list
    ]
    if strategy == "aggregate_mean":
        numeric_values = grouped_values.mean()
    else:
        numeric_values = grouped_values.median()

    aggregated = metadata.join(numeric_values)
    aggregated.index.name = SITE_MATRIX_ID_COLUMN
    aggregated = aggregated.reset_index()
    return aggregated.loc[:, ordered_columns].copy().reset_index(drop=True)
