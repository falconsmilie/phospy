"""Site-matrix construction stage for dataset preprocessing."""

from __future__ import annotations

import re
from dataclasses import replace

import pandas as pd

from phospy.api.configs import (
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEAN,
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEDIAN,
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_ERROR,
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_FIRST,
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_MAX_MEAN_SIGNAL,
    DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
    DATASET_SITE_MATRIX_POLICY_AS_INPUT,
    DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA,
)
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    PreprocessingState,
)
from phospy.errors.input import PhosPyInputError

_GENE_SYMBOL_COLUMN = "gene_symbol"
_SITE_COLUMN = "site"
_SITE_SEQUENCE_COLUMN = "site_sequence"
_SITE_ID_COLUMN = "site_id"
_REQUIRED_SITE_METADATA_COLUMNS = (
    _GENE_SYMBOL_COLUMN,
    _SITE_COLUMN,
)
_SITE_TOKEN_PATTERN = re.compile(r"^[A-Za-z]+\d+$")
_GENE_TOKEN_PATTERN = re.compile(r"^[^;\s]+$")
_ROW_DROP_STATS_ATTR = "site_matrix_row_drop_stats"
_SITE_MATRIX_POLICY_ATTR = "site_matrix_policy"
_INTERNAL_SITE_MATRIX_MISSING_DATA_POLICY_RETAIN_MISSING = "retain_missing"
_INTERNAL_SITE_MATRIX_MISSING_DATA_POLICY_REQUIRE_MIN_OBSERVED_VALUES = (
    "require_min_observed_values"
)
_SUPPORTED_DUPLICATE_SITE_STRATEGIES = {
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_MAX_MEAN_SIGNAL,
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_FIRST,
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEAN,
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEDIAN,
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_ERROR,
}


class SiteMatrixStage:
    """Build site-matrix-ready phospho rows from site metadata when requested.

    This stage ports the legacy site-matrix policy surface behind
    `site_matrix.policy='build_from_metadata'`.
    """

    stage_key = DATASET_PREPROCESSING_STAGE_SITE_MATRIX

    def run(self, state: PreprocessingState) -> PreprocessingState:
        policy = state.plan.site_matrix_policy
        if policy == DATASET_SITE_MATRIX_POLICY_AS_INPUT:
            return state
        if policy != DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA:
            raise PhosPyInputError(
                "dataset build request preprocessing_config contains an unsupported "
                "site_matrix.policy"
            )

        self._require_site_metadata_columns(state.site_metadata)

        gene_symbol = _resolve_required_string_column(
            state.site_metadata,
            column_name=_GENE_SYMBOL_COLUMN,
        )
        site = _resolve_required_string_column(
            state.site_metadata,
            column_name=_SITE_COLUMN,
        )
        if _SITE_SEQUENCE_COLUMN in state.site_metadata.columns:
            site_sequence = _resolve_optional_string_column(
                state.site_metadata,
                column_name=_SITE_SEQUENCE_COLUMN,
            )
        else:
            site_sequence = pd.Series(
                pd.NA,
                index=state.site_metadata.index.copy(),
                dtype="string",
                name=_SITE_SEQUENCE_COLUMN,
            )
        constructed_site_id = _build_site_identifier(gene_symbol=gene_symbol, site=site)

        has_sequence = site_sequence.notna()
        with_sequence_phospho = state.phospho.loc[has_sequence]
        with_sequence_site_metadata = state.site_metadata.loc[has_sequence]
        with_sequence_site_id = constructed_site_id.loc[has_sequence]
        dropped_missing_sequence = int((~has_sequence).sum())

        (
            policy_filtered_phospho,
            dropped_incomplete_values,
            required_observed_count,
        ) = _apply_missing_data_policy(
            phospho=with_sequence_phospho,
            missing_data_policy=state.plan.site_matrix_missing_data_policy,
            minimum_observed_values=state.plan.site_matrix_minimum_observed_values,
        )
        policy_filtered_site_metadata = with_sequence_site_metadata.loc[
            policy_filtered_phospho.index
        ]
        policy_filtered_site_id = with_sequence_site_id.loc[
            policy_filtered_phospho.index
        ]

        (
            deduplicated_phospho,
            deduplicated_site_metadata,
            deduplicated_site_rows,
        ) = _apply_duplicate_site_policy(
            phospho=policy_filtered_phospho,
            site_metadata=policy_filtered_site_metadata,
            constructed_site_id=policy_filtered_site_id,
            duplicate_site_strategy=state.plan.site_matrix_duplicate_site_strategy,
        )

        final_phospho = deduplicated_phospho.sort_index(kind="stable")
        final_site_metadata = deduplicated_site_metadata.reindex(final_phospho.index)
        final_site_index = pd.Index(
            final_phospho.index.tolist(), name=state.phospho.index.name
        )
        final_phospho.index = final_site_index
        final_site_metadata.index = final_site_index.copy()

        row_drop_stats = {
            "input_rows": int(len(state.phospho.index)),
            "dropped_missing_sequence": dropped_missing_sequence,
            "dropped_incomplete_values": dropped_incomplete_values,
            "missing_data_policy": state.plan.site_matrix_missing_data_policy,
            "required_observed_count": required_observed_count,
            "deduplicated_site_rows": deduplicated_site_rows,
            "duplicate_site_strategy": state.plan.site_matrix_duplicate_site_strategy,
            "retained_rows": int(len(final_phospho.index)),
        }
        if final_phospho.empty:
            diagnostics = _format_row_drop_diagnostics(row_drop_stats)
            raise PhosPyInputError(
                "dataset build request preprocessing site-matrix construction "
                f"produced no retained rows after filtering; {diagnostics}"
            )

        final_phospho.attrs[_ROW_DROP_STATS_ATTR] = row_drop_stats.copy()
        final_site_metadata.attrs[_ROW_DROP_STATS_ATTR] = row_drop_stats.copy()
        final_phospho.attrs[_SITE_MATRIX_POLICY_ATTR] = policy
        final_site_metadata.attrs[_SITE_MATRIX_POLICY_ATTR] = policy

        return replace(
            state,
            phospho=final_phospho,
            site_metadata=final_site_metadata,
        )

    @staticmethod
    def _require_site_metadata_columns(site_metadata: pd.DataFrame) -> None:
        missing_columns = [
            column
            for column in _REQUIRED_SITE_METADATA_COLUMNS
            if column not in site_metadata.columns
        ]
        if missing_columns:
            joined_missing_columns = ", ".join(missing_columns)
            raise PhosPyInputError(
                "dataset build request preprocessing site-matrix construction "
                "requires site_metadata columns: "
                f"{joined_missing_columns}"
            )


def _resolve_required_string_column(
    site_metadata: pd.DataFrame,
    *,
    column_name: str,
) -> pd.Series:
    column = site_metadata.loc[:, column_name]
    normalized = column.astype("string").str.strip()
    invalid_mask = column.isna() | normalized.isna() | (normalized == "")
    if bool(invalid_mask.any()):
        raise PhosPyInputError(
            "dataset build request preprocessing site-matrix construction requires "
            f"site_metadata.{column_name} to contain non-empty values"
        )
    return normalized.astype(str)


def _resolve_optional_string_column(
    site_metadata: pd.DataFrame,
    *,
    column_name: str,
) -> pd.Series:
    column = site_metadata.loc[:, column_name]
    normalized = column.astype("string").str.strip()
    missing_mask = column.isna() | normalized.isna() | (normalized == "")
    return normalized.where(~missing_mask, other=pd.NA)


def _build_site_identifier(
    *,
    gene_symbol: pd.Series,
    site: pd.Series,
) -> pd.Series:
    normalized_gene_symbol = gene_symbol.astype(str).str.strip().str.upper()
    normalized_site = site.astype(str).str.strip().str.upper()

    invalid_gene_symbol = ~normalized_gene_symbol.str.fullmatch(_GENE_TOKEN_PATTERN)
    if bool(invalid_gene_symbol.any()):
        preview = ", ".join(
            normalized_gene_symbol.loc[invalid_gene_symbol].astype(str).head(3).tolist()
        )
        raise PhosPyInputError(
            "dataset build request preprocessing site-matrix construction requires "
            "site_metadata.gene_symbol values that normalize to canonical non-empty "
            f"tokens without whitespace/semicolons; example invalid values: {preview}"
        )

    invalid_site = ~normalized_site.str.fullmatch(_SITE_TOKEN_PATTERN)
    if bool(invalid_site.any()):
        preview = ", ".join(
            normalized_site.loc[invalid_site].astype(str).head(3).tolist()
        )
        raise PhosPyInputError(
            "dataset build request preprocessing site-matrix construction requires "
            "site_metadata.site values that normalize to site tokens like 'S123'; "
            f"example invalid values: {preview}"
        )

    site_id = normalized_gene_symbol + ";" + normalized_site + ";"
    site_id.index = gene_symbol.index.copy()
    site_id.name = _SITE_ID_COLUMN
    return site_id


def _apply_missing_data_policy(
    *,
    phospho: pd.DataFrame,
    missing_data_policy: str,
    minimum_observed_values: int | None,
) -> tuple[pd.DataFrame, int, int]:
    if missing_data_policy == DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING:
        retained_mask = phospho.notna().all(axis=1)
        required_observed_count = phospho.shape[1]
    elif (
        missing_data_policy == _INTERNAL_SITE_MATRIX_MISSING_DATA_POLICY_RETAIN_MISSING
    ):
        retained_mask = pd.Series(True, index=phospho.index)
        required_observed_count = 0
    elif (
        missing_data_policy
        == _INTERNAL_SITE_MATRIX_MISSING_DATA_POLICY_REQUIRE_MIN_OBSERVED_VALUES
    ):
        if minimum_observed_values is None:
            raise PhosPyInputError(
                "dataset build request preprocessing site-matrix construction requires "
                "minimum_observed_values when "
                "site_matrix.missing_data_policy='require_min_observed_values'"
            )
        if minimum_observed_values > phospho.shape[1]:
            raise PhosPyInputError(
                "dataset build request preprocessing site-matrix construction "
                "minimum_observed_values cannot exceed phospho sample count "
                f"({phospho.shape[1]})"
            )
        required_observed_count = minimum_observed_values
        retained_mask = phospho.notna().sum(axis=1) >= required_observed_count
    else:
        raise PhosPyInputError(
            "dataset build request preprocessing_config contains an unsupported "
            "site_matrix.missing_data_policy"
        )

    filtered = phospho.loc[retained_mask]
    dropped_rows = int(len(phospho.index) - len(filtered.index))
    return filtered, dropped_rows, required_observed_count


def _apply_duplicate_site_policy(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    constructed_site_id: pd.Series,
    duplicate_site_strategy: str,
) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    if duplicate_site_strategy not in _SUPPORTED_DUPLICATE_SITE_STRATEGIES:
        raise PhosPyInputError(
            "dataset build request preprocessing_config contains an unsupported "
            "site_matrix.duplicate_site_strategy"
        )

    if phospho.empty:
        empty_site_index = pd.Index([], name=_SITE_ID_COLUMN)
        empty_phospho = phospho.copy()
        empty_site_metadata = site_metadata.copy()
        empty_phospho.index = empty_site_index
        empty_site_metadata.index = empty_site_index.copy()
        return empty_phospho, empty_site_metadata, 0

    duplicate_mask = constructed_site_id.duplicated(keep=False)
    if not bool(duplicate_mask.any()):
        final_site_index = pd.Index(
            constructed_site_id.astype(str).tolist(), name=_SITE_ID_COLUMN
        )
        direct_phospho = phospho.copy()
        direct_site_metadata = site_metadata.copy()
        direct_phospho.index = final_site_index
        direct_site_metadata.index = final_site_index.copy()
        return direct_phospho, direct_site_metadata, 0

    if duplicate_site_strategy == DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_ERROR:
        duplicate_sites = (
            constructed_site_id.loc[duplicate_mask]
            .astype(str)
            .drop_duplicates()
            .head(3)
        )
        preview = ", ".join(duplicate_sites.tolist())
        raise PhosPyInputError(
            "dataset build request preprocessing site-matrix construction found "
            "duplicate constructed site identifiers and "
            "site_matrix.duplicate_site_strategy='error': "
            f"{preview}"
        )

    if duplicate_site_strategy == DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_FIRST:
        selected_rows = (
            pd.DataFrame({_SITE_ID_COLUMN: constructed_site_id}, index=phospho.index)
            .drop_duplicates(_SITE_ID_COLUMN, keep="first")
            .index
        )
        selected_phospho = phospho.loc[selected_rows].copy()
        selected_site_metadata = site_metadata.loc[selected_rows].copy()
        selected_site_ids = constructed_site_id.loc[selected_rows]
        final_site_index = pd.Index(
            selected_site_ids.astype(str).tolist(), name=_SITE_ID_COLUMN
        )
        selected_phospho.index = final_site_index
        selected_site_metadata.index = final_site_index.copy()
        return (
            selected_phospho,
            selected_site_metadata,
            int(len(phospho.index) - len(selected_phospho.index)),
        )

    if (
        duplicate_site_strategy
        == DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_MAX_MEAN_SIGNAL
    ):
        value_columns = list(phospho.columns)
        dedupe_work = pd.DataFrame(
            {
                _SITE_ID_COLUMN: constructed_site_id,
                "observed_values": phospho.loc[:, value_columns].notna().sum(axis=1),
                "mean_signal": phospho.loc[:, value_columns].mean(axis=1, skipna=True),
                "row_order": range(len(phospho)),
            },
            index=phospho.index,
        )
        selected_rows = (
            dedupe_work.sort_values(
                [_SITE_ID_COLUMN, "observed_values", "mean_signal", "row_order"],
                ascending=[True, False, False, True],
                kind="stable",
                na_position="last",
            )
            .drop_duplicates(_SITE_ID_COLUMN, keep="first")
            .index
        )
        selected_phospho = phospho.loc[selected_rows].copy()
        selected_site_metadata = site_metadata.loc[selected_rows].copy()
        selected_site_ids = constructed_site_id.loc[selected_rows]
        final_site_index = pd.Index(
            selected_site_ids.astype(str).tolist(), name=_SITE_ID_COLUMN
        )
        selected_phospho.index = final_site_index
        selected_site_metadata.index = final_site_index.copy()
        return (
            selected_phospho,
            selected_site_metadata,
            int(len(phospho.index) - len(selected_phospho.index)),
        )

    if duplicate_site_strategy in {
        DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEAN,
        DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEDIAN,
    }:
        metadata_columns = list(site_metadata.columns)
        grouped_metadata = (
            site_metadata.assign(**{_SITE_ID_COLUMN: constructed_site_id.to_numpy()})
            .groupby(_SITE_ID_COLUMN, sort=False)[metadata_columns]
            .first()
        )
        grouped_values = phospho.groupby(constructed_site_id, sort=False)
        if (
            duplicate_site_strategy
            == DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEAN
        ):
            grouped_phospho = grouped_values.mean()
        else:
            grouped_phospho = grouped_values.median()
        grouped_phospho.index = pd.Index(
            grouped_phospho.index.astype(str), name=_SITE_ID_COLUMN
        )
        grouped_metadata.index = pd.Index(
            grouped_metadata.index.astype(str), name=_SITE_ID_COLUMN
        )
        return (
            grouped_phospho,
            grouped_metadata,
            int(len(phospho.index) - len(grouped_phospho.index)),
        )

    raise RuntimeError("site-matrix duplicate strategy dispatch fell through")


def _format_row_drop_diagnostics(row_drop_stats: dict[str, int | str]) -> str:
    known_drops = (
        int(row_drop_stats.get("dropped_missing_sequence", 0))
        + int(row_drop_stats.get("dropped_incomplete_values", 0))
        + int(row_drop_stats.get("deduplicated_site_rows", 0))
    )
    input_rows = int(row_drop_stats.get("input_rows", 0))
    retained_rows = int(row_drop_stats.get("retained_rows", 0))
    other_dropped_rows = max(input_rows - retained_rows - known_drops, 0)
    return (
        f"input_rows={input_rows}, "
        "dropped_missing_sequence="
        f"{int(row_drop_stats.get('dropped_missing_sequence', 0))}, "
        "dropped_incomplete_values="
        f"{int(row_drop_stats.get('dropped_incomplete_values', 0))}, "
        "missing_data_policy="
        f"{str(row_drop_stats.get('missing_data_policy', 'drop_any_missing'))}, "
        "required_observed_count="
        f"{int(row_drop_stats.get('required_observed_count', 0))}, "
        "deduplicated_site_rows="
        f"{int(row_drop_stats.get('deduplicated_site_rows', 0))}, "
        "duplicate_site_strategy="
        f"{str(row_drop_stats.get('duplicate_site_strategy', 'max_mean_signal'))}, "
        f"other_dropped_rows={other_dropped_rows}, "
        f"retained_rows={retained_rows}"
    )


__all__ = ["SiteMatrixStage"]
