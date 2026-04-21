"""Site-matrix construction stage for dataset preprocessing."""

from __future__ import annotations

from dataclasses import replace

import pandas as pd

from phospy.api.configs import (
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
_REQUIRED_SITE_METADATA_COLUMNS = (
    _GENE_SYMBOL_COLUMN,
    _SITE_COLUMN,
    _SITE_SEQUENCE_COLUMN,
)


class SiteMatrixStage:
    """Build site-matrix-ready phospho rows from site metadata when requested."""

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
        site_sequence = _resolve_optional_string_column(
            state.site_metadata,
            column_name=_SITE_SEQUENCE_COLUMN,
        )
        constructed_site_id = _build_site_identifier(gene_symbol=gene_symbol, site=site)

        has_sequence = site_sequence.notna()
        has_complete_values = state.phospho.notna().all(axis=1)
        retained_mask = has_sequence & has_complete_values

        retained_rows = int(retained_mask.sum())
        if retained_rows == 0:
            dropped_missing_sequence = int((~has_sequence).sum())
            dropped_incomplete_values = int((~has_complete_values).sum())
            raise PhosPyInputError(
                "dataset build request preprocessing site-matrix construction "
                "produced no retained rows after filtering; "
                f"dropped_missing_sequence={dropped_missing_sequence}, "
                f"dropped_incomplete_values={dropped_incomplete_values}"
            )

        filtered_phospho = state.phospho.loc[retained_mask]
        filtered_site_metadata = state.site_metadata.loc[retained_mask]
        filtered_constructed_site_id = constructed_site_id.loc[retained_mask]

        selected_rows = _resolve_selected_rows(
            phospho=filtered_phospho,
            constructed_site_id=filtered_constructed_site_id,
        )

        deduplicated_phospho = filtered_phospho.loc[selected_rows]
        deduplicated_site_metadata = filtered_site_metadata.loc[selected_rows].copy()
        deduplicated_site_id = filtered_constructed_site_id.loc[selected_rows]

        ordered_site_id = deduplicated_site_id.sort_values(kind="stable")
        ordered_rows = ordered_site_id.index
        final_site_index = pd.Index(
            ordered_site_id.tolist(),
            name=state.phospho.index.name,
        )

        final_phospho = deduplicated_phospho.loc[ordered_rows]
        final_phospho.index = final_site_index

        final_site_metadata = deduplicated_site_metadata.loc[ordered_rows]
        final_site_metadata.index = final_site_index

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
    site_id = gene_symbol.astype(str) + ";" + site.astype(str) + ";"
    site_id.index = gene_symbol.index.copy()
    site_id.name = "site_id"
    return site_id


def _resolve_selected_rows(
    *,
    phospho: pd.DataFrame,
    constructed_site_id: pd.Series,
) -> pd.Index:
    if phospho.empty:
        return phospho.index.copy()
    duplicate_mask = constructed_site_id.duplicated(keep=False)
    if not bool(duplicate_mask.any()):
        return phospho.index.copy()

    value_columns = list(phospho.columns)
    dedupe_work = pd.DataFrame(
        {
            "site_id": constructed_site_id,
            "observed_values": phospho.loc[:, value_columns].notna().sum(axis=1),
            "mean_signal": phospho.loc[:, value_columns].mean(axis=1, skipna=True),
            "row_order": range(len(phospho)),
        },
        index=phospho.index,
    )
    selected = dedupe_work.sort_values(
        ["site_id", "observed_values", "mean_signal", "row_order"],
        ascending=[True, False, False, True],
        kind="stable",
        na_position="last",
    ).drop_duplicates("site_id", keep="first")
    return pd.Index(selected.index.copy())


__all__ = ["SiteMatrixStage"]
