from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd

from ..errors import TableSchemaError
from ..internal.constants import (
    PHOSPHO_GENE_COLUMN,
    PHOSPHO_UID_COLUMN,
    ROW_DROP_DEDUPLICATED_SITE_ROWS_KEY,
    ROW_DROP_DROPPED_INCOMPLETE_VALUES_KEY,
    ROW_DROP_DROPPED_MISSING_SEQUENCE_KEY,
    ROW_DROP_DUPLICATE_SITE_STRATEGY_KEY,
    ROW_DROP_INPUT_ROWS_KEY,
    ROW_DROP_MISSING_DATA_POLICY_KEY,
    ROW_DROP_REQUIRED_OBSERVED_COUNT_KEY,
    ROW_DROP_RETAINED_ROWS_KEY,
    ROW_DROP_STATS_ATTR,
    SITE_MATRIX_GENE_COLUMN,
    SITE_MATRIX_ID_COLUMN,
    SITE_MATRIX_P_SITE_COLUMN,
)
from ..internal.types import (
    DUPLICATE_SITE_STRATEGY_AGGREGATE_MEAN,
    DUPLICATE_SITE_STRATEGY_AGGREGATE_MEDIAN,
    DUPLICATE_SITE_STRATEGY_ERROR,
    DUPLICATE_SITE_STRATEGY_FIRST,
    DUPLICATE_SITE_STRATEGY_MAX_MEAN_SIGNAL,
    SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
    SITE_MATRIX_MISSING_DATA_POLICY_REQUIRE_MIN_OBSERVED_VALUES,
    SITE_MATRIX_MISSING_DATA_POLICY_RETAIN_MISSING,
    DuplicateSiteStrategy,
    SiteMatrixMissingDataPolicy,
)
from ..validation.schema.tables import SiteMatrixSourceSchema
from ..validation.values.enums import (
    validate_duplicate_site_strategy,
    validate_site_matrix_missing_data_policy,
)
from ..validation.values.identifiers import build_canonical_site_id
from ..validation.values.numeric import validate_positive_int

__all__ = [
    "DEFAULT_SITE_MATRIX_POLICY",
    "SiteMatrixPolicy",
    "build_site_matrix",
    "format_row_drop_diagnostics",
]


@dataclass(frozen=True, slots=True)
class SiteMatrixPolicy:
    """Explicit policy for duplicate phosphosite rows during site-matrix building.

    ``missing_data_policy`` controls how missing corrected values are handled:

    - ``"drop_any_missing"`` keeps only complete rows across all value columns
    - ``"retain_missing"`` keeps rows even when some corrected values are missing
    - ``"require_min_observed_values"`` keeps rows with at least
      ``minimum_observed_values`` non-missing corrected values

    ``duplicate_site_strategy`` controls how multiple corrected phosphosite rows
    that resolve to the same site identifier are collapsed:

    - ``"max_mean_signal"`` keeps the row with the largest mean corrected signal
    - ``"first"`` keeps the first encountered complete row
    - ``"aggregate_mean"`` averages corrected values across duplicate rows
    - ``"aggregate_median"`` takes the column-wise median across duplicate rows
    - ``"error"`` rejects duplicate site IDs instead of collapsing them
    """

    duplicate_site_strategy: DuplicateSiteStrategy = (
        DUPLICATE_SITE_STRATEGY_MAX_MEAN_SIGNAL
    )
    missing_data_policy: SiteMatrixMissingDataPolicy = (
        SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING
    )
    minimum_observed_values: int | None = None

    def __post_init__(self) -> None:
        validate_duplicate_site_strategy(self.duplicate_site_strategy)
        validate_site_matrix_missing_data_policy(self.missing_data_policy)

        if (
            self.missing_data_policy
            == SITE_MATRIX_MISSING_DATA_POLICY_REQUIRE_MIN_OBSERVED_VALUES
        ):
            if self.minimum_observed_values is None:
                msg = (
                    "minimum_observed_values is required when "
                    "missing_data_policy='require_min_observed_values'"
                )
                raise ValueError(msg)
            validate_positive_int(
                self.minimum_observed_values,
                name="minimum_observed_values",
            )
            return

        if self.minimum_observed_values is not None:
            msg = (
                "minimum_observed_values may only be provided when "
                "missing_data_policy='require_min_observed_values'"
            )
            raise ValueError(msg)

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
        ROW_DROP_INPUT_ROWS_KEY: int(row_drop_stats.get(ROW_DROP_INPUT_ROWS_KEY, 0)),
        ROW_DROP_DROPPED_MISSING_SEQUENCE_KEY: int(
            row_drop_stats.get(ROW_DROP_DROPPED_MISSING_SEQUENCE_KEY, 0)
        ),
        ROW_DROP_DROPPED_INCOMPLETE_VALUES_KEY: int(
            row_drop_stats.get(ROW_DROP_DROPPED_INCOMPLETE_VALUES_KEY, 0)
        ),
        ROW_DROP_MISSING_DATA_POLICY_KEY: str(
            row_drop_stats.get(
                ROW_DROP_MISSING_DATA_POLICY_KEY,
                DEFAULT_SITE_MATRIX_POLICY.missing_data_policy,
            )
        ),
        ROW_DROP_REQUIRED_OBSERVED_COUNT_KEY: int(
            row_drop_stats.get(ROW_DROP_REQUIRED_OBSERVED_COUNT_KEY, 0)
        ),
        ROW_DROP_DEDUPLICATED_SITE_ROWS_KEY: int(
            row_drop_stats.get(ROW_DROP_DEDUPLICATED_SITE_ROWS_KEY, 0)
        ),
        ROW_DROP_RETAINED_ROWS_KEY: int(
            row_drop_stats.get(ROW_DROP_RETAINED_ROWS_KEY, 0)
        ),
        ROW_DROP_DUPLICATE_SITE_STRATEGY_KEY: str(
            row_drop_stats.get(
                ROW_DROP_DUPLICATE_SITE_STRATEGY_KEY,
                DEFAULT_SITE_MATRIX_POLICY.duplicate_site_strategy,
            )
        ),
    }
    known_drops = (
        stats[ROW_DROP_DROPPED_MISSING_SEQUENCE_KEY]
        + stats[ROW_DROP_DROPPED_INCOMPLETE_VALUES_KEY]
        + stats[ROW_DROP_DEDUPLICATED_SITE_ROWS_KEY]
    )
    other_dropped_rows = max(
        stats[ROW_DROP_INPUT_ROWS_KEY]
        - stats[ROW_DROP_RETAINED_ROWS_KEY]
        - known_drops,
        0,
    )
    return (
        "row-drop diagnostics: "
        f"{ROW_DROP_INPUT_ROWS_KEY}={stats[ROW_DROP_INPUT_ROWS_KEY]}, "
        f"{ROW_DROP_DROPPED_MISSING_SEQUENCE_KEY}="
        f"{stats[ROW_DROP_DROPPED_MISSING_SEQUENCE_KEY]}, "
        f"{ROW_DROP_DROPPED_INCOMPLETE_VALUES_KEY}="
        f"{stats[ROW_DROP_DROPPED_INCOMPLETE_VALUES_KEY]}, "
        f"{ROW_DROP_MISSING_DATA_POLICY_KEY}={stats[ROW_DROP_MISSING_DATA_POLICY_KEY]}, "
        f"{ROW_DROP_REQUIRED_OBSERVED_COUNT_KEY}="
        f"{stats[ROW_DROP_REQUIRED_OBSERVED_COUNT_KEY]}, "
        f"{ROW_DROP_DEDUPLICATED_SITE_ROWS_KEY}={stats[ROW_DROP_DEDUPLICATED_SITE_ROWS_KEY]}, "
        f"{ROW_DROP_DUPLICATE_SITE_STRATEGY_KEY}="
        f"{stats[ROW_DROP_DUPLICATE_SITE_STRATEGY_KEY]}, "
        f"other_dropped_rows={other_dropped_rows}, "
        f"{ROW_DROP_RETAINED_ROWS_KEY}={stats[ROW_DROP_RETAINED_ROWS_KEY]}"
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
    missing_data_policy: SiteMatrixMissingDataPolicy | None = None,
    minimum_observed_values: int | None = None,
    copy_frame: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    resolved_policy = _resolve_site_matrix_policy(
        policy=policy,
        duplicate_site_strategy=duplicate_site_strategy,
        missing_data_policy=missing_data_policy,
        minimum_observed_values=minimum_observed_values,
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
    site_work[SITE_MATRIX_ID_COLUMN] = [
        build_canonical_site_id(
            entity=gene,
            site_token=site,
            context="site-matrix source table",
        )
        for gene, site in zip(
            site_work[gene_col_name],
            site_work[p_site_col_name],
            strict=True,
        )
    ]

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

    policy_filtered_rows, dropped_incomplete_values, required_observed_count = (
        _apply_missing_data_policy(
            with_sequence=with_sequence,
            value_cols=value_cols,
            policy=resolved_policy,
        )
    )

    phosr_input = _apply_duplicate_site_policy(
        source_rows=policy_filtered_rows,
        value_cols=value_cols,
        policy=resolved_policy,
    )
    phosr_input = phosr_input.sort_values(
        SITE_MATRIX_ID_COLUMN,
        ascending=True,
        kind="stable",
    ).reset_index(drop=True)
    deduplicated_site_rows = len(policy_filtered_rows) - len(phosr_input)

    row_drop_stats: dict[str, int | str] = {
        ROW_DROP_INPUT_ROWS_KEY: total_rows,
        ROW_DROP_DROPPED_MISSING_SEQUENCE_KEY: dropped_missing_sequence,
        ROW_DROP_DROPPED_INCOMPLETE_VALUES_KEY: dropped_incomplete_values,
        ROW_DROP_MISSING_DATA_POLICY_KEY: resolved_policy.missing_data_policy,
        ROW_DROP_REQUIRED_OBSERVED_COUNT_KEY: required_observed_count,
        ROW_DROP_DEDUPLICATED_SITE_ROWS_KEY: deduplicated_site_rows,
        ROW_DROP_DUPLICATE_SITE_STRATEGY_KEY: resolved_policy.duplicate_site_strategy,
        ROW_DROP_RETAINED_ROWS_KEY: len(phosr_input),
    }

    matrix = phosr_input.loc[:, [SITE_MATRIX_ID_COLUMN, *value_cols]].set_index(
        SITE_MATRIX_ID_COLUMN
    )
    matrix.attrs[ROW_DROP_STATS_ATTR] = row_drop_stats.copy()
    sequences = phosr_input.set_index(SITE_MATRIX_ID_COLUMN)[sequence_col]
    sequences.attrs[ROW_DROP_STATS_ATTR] = row_drop_stats.copy()
    phosr_input.attrs[ROW_DROP_STATS_ATTR] = row_drop_stats.copy()
    return phosr_input, matrix, sequences


def _resolve_site_matrix_policy(
    *,
    policy: SiteMatrixPolicy | None,
    duplicate_site_strategy: DuplicateSiteStrategy | None,
    missing_data_policy: SiteMatrixMissingDataPolicy | None,
    minimum_observed_values: int | None,
) -> SiteMatrixPolicy:
    has_policy_overrides = any(
        (
            duplicate_site_strategy is not None,
            missing_data_policy is not None,
            minimum_observed_values is not None,
        )
    )
    if policy is not None and has_policy_overrides:
        msg = (
            "build_site_matrix() accepts either policy or individual site-matrix "
            "policy overrides, not both"
        )
        raise ValueError(msg)
    if policy is not None:
        return SiteMatrixPolicy.from_value(policy)

    if not has_policy_overrides:
        return DEFAULT_SITE_MATRIX_POLICY

    resolved_missing_data_policy = missing_data_policy
    if resolved_missing_data_policy is None and minimum_observed_values is not None:
        resolved_missing_data_policy = (
            SITE_MATRIX_MISSING_DATA_POLICY_REQUIRE_MIN_OBSERVED_VALUES
        )

    return SiteMatrixPolicy(
        duplicate_site_strategy=(
            duplicate_site_strategy
            or DEFAULT_SITE_MATRIX_POLICY.duplicate_site_strategy
        ),
        missing_data_policy=(
            resolved_missing_data_policy
            or DEFAULT_SITE_MATRIX_POLICY.missing_data_policy
        ),
        minimum_observed_values=minimum_observed_values,
    )


def _apply_missing_data_policy(
    *,
    with_sequence: pd.DataFrame,
    value_cols: Sequence[str],
    policy: SiteMatrixPolicy,
) -> tuple[pd.DataFrame, int, int]:
    value_col_list = list(value_cols)
    required_observed_count = _resolve_required_observed_count(
        policy=policy,
        value_column_count=len(value_col_list),
    )

    if policy.missing_data_policy == SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING:
        filtered_rows = with_sequence.loc[
            with_sequence.loc[:, value_col_list].notna().all(axis=1)
        ]
    elif policy.missing_data_policy == SITE_MATRIX_MISSING_DATA_POLICY_RETAIN_MISSING:
        filtered_rows = with_sequence
    elif (
        policy.missing_data_policy
        == SITE_MATRIX_MISSING_DATA_POLICY_REQUIRE_MIN_OBSERVED_VALUES
    ):
        observed_counts = with_sequence.loc[:, value_col_list].notna().sum(axis=1)
        filtered_rows = with_sequence.loc[observed_counts >= required_observed_count]
    else:
        msg = f"Unsupported missing_data_policy: {policy.missing_data_policy}"
        raise TableSchemaError(msg)

    dropped_rows = len(with_sequence) - len(filtered_rows)
    return (
        filtered_rows.copy().reset_index(drop=True),
        dropped_rows,
        required_observed_count,
    )


def _resolve_required_observed_count(
    *,
    policy: SiteMatrixPolicy,
    value_column_count: int,
) -> int:
    if policy.missing_data_policy == SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING:
        return value_column_count
    if policy.missing_data_policy == SITE_MATRIX_MISSING_DATA_POLICY_RETAIN_MISSING:
        return 0
    if (
        policy.missing_data_policy
        != SITE_MATRIX_MISSING_DATA_POLICY_REQUIRE_MIN_OBSERVED_VALUES
    ):
        msg = f"Unsupported missing_data_policy: {policy.missing_data_policy}"
        raise TableSchemaError(msg)

    required = policy.minimum_observed_values
    if required is None:
        msg = (
            "minimum_observed_values is required when "
            "missing_data_policy='require_min_observed_values'"
        )
        raise ValueError(msg)
    if required > value_column_count:
        msg = (
            "minimum_observed_values cannot exceed the number of site-matrix "
            f"value columns ({value_column_count})"
        )
        raise ValueError(msg)
    return required


def _apply_duplicate_site_policy(
    *,
    source_rows: pd.DataFrame,
    value_cols: Sequence[str],
    policy: SiteMatrixPolicy,
) -> pd.DataFrame:
    if source_rows.empty:
        return source_rows.copy().reset_index(drop=True)

    duplicate_mask = source_rows.duplicated(SITE_MATRIX_ID_COLUMN, keep=False)
    if not bool(duplicate_mask.any()):
        return source_rows.copy().reset_index(drop=True)

    if policy.duplicate_site_strategy == DUPLICATE_SITE_STRATEGY_ERROR:
        duplicate_sites = (
            source_rows.loc[duplicate_mask, SITE_MATRIX_ID_COLUMN]
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

    if policy.duplicate_site_strategy == DUPLICATE_SITE_STRATEGY_MAX_MEAN_SIGNAL:
        work = source_rows.copy()
        work["__observed_values"] = work.loc[:, list(value_cols)].notna().sum(axis=1)
        work["__mean_signal"] = work.loc[:, list(value_cols)].mean(axis=1, skipna=True)
        work["__row_order"] = range(len(work))
        selected = work.sort_values(
            [
                SITE_MATRIX_ID_COLUMN,
                "__observed_values",
                "__mean_signal",
                "__row_order",
            ],
            ascending=[True, False, False, True],
            kind="stable",
            na_position="last",
        ).drop_duplicates(SITE_MATRIX_ID_COLUMN, keep="first")
        return (
            selected.drop(columns=["__observed_values", "__mean_signal", "__row_order"])
            .copy()
            .reset_index(drop=True)
        )

    if policy.duplicate_site_strategy == DUPLICATE_SITE_STRATEGY_FIRST:
        return (
            source_rows.drop_duplicates(SITE_MATRIX_ID_COLUMN, keep="first")
            .copy()
            .reset_index(drop=True)
        )

    if policy.duplicate_site_strategy in {
        DUPLICATE_SITE_STRATEGY_AGGREGATE_MEAN,
        DUPLICATE_SITE_STRATEGY_AGGREGATE_MEDIAN,
    }:
        return _aggregate_duplicate_sites(
            source_rows=source_rows,
            value_cols=value_cols,
            strategy=policy.duplicate_site_strategy,
        )

    msg = f"Unsupported duplicate_site_strategy: {policy.duplicate_site_strategy}"
    raise TableSchemaError(msg)


def _aggregate_duplicate_sites(
    *,
    source_rows: pd.DataFrame,
    value_cols: Sequence[str],
    strategy: DuplicateSiteStrategy,
) -> pd.DataFrame:
    value_col_list = list(value_cols)
    ordered_columns = list(source_rows.columns)
    metadata_columns = [
        column
        for column in ordered_columns
        if column not in {SITE_MATRIX_ID_COLUMN, *value_col_list}
    ]

    metadata = (
        source_rows.loc[:, [SITE_MATRIX_ID_COLUMN, *metadata_columns]]
        .groupby(SITE_MATRIX_ID_COLUMN, sort=False)
        .first()
    )
    grouped_values = source_rows.groupby(SITE_MATRIX_ID_COLUMN, sort=False)[
        value_col_list
    ]
    if strategy == DUPLICATE_SITE_STRATEGY_AGGREGATE_MEAN:
        numeric_values = grouped_values.mean()
    else:
        numeric_values = grouped_values.median()

    aggregated = metadata.join(numeric_values)
    aggregated.index.name = SITE_MATRIX_ID_COLUMN
    aggregated = aggregated.reset_index()
    return aggregated.loc[:, ordered_columns].copy().reset_index(drop=True)
