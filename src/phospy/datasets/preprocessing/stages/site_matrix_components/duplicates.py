from __future__ import annotations

from typing import cast

import pandas as pd

from phospy.api.configs import (
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_ERROR,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL,
)
from phospy.datasets.preprocessing.models import DuplicateSiteResolutionResult
from phospy.datasets.preprocessing.report_schema import (
    DUPLICATE_SITE_RESOLUTION_COLUMNS,
    dataframe_from_duplicate_site_resolution_rows,
)
from phospy.datasets.preprocessing.stages.site_matrix_components.metadata import (
    MetadataConflictDetector,
    _empty_metadata_conflicts,
    _resolve_source_metadata_column,
)
from phospy.errors.input import PhosPyInputError

_SITE_ID_COLUMN = "site_id"
_SUPPORTED_DUPLICATE_SITE_POLICIES = {
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_ERROR,
}


class DuplicateSiteResolver:
    """Apply duplicate-site policy and build duplicate/conflict report tables."""

    def __init__(
        self,
        *,
        metadata_conflict_detector: MetadataConflictDetector | None = None,
    ) -> None:
        self._metadata_conflict_detector = (
            MetadataConflictDetector()
            if metadata_conflict_detector is None
            else metadata_conflict_detector
        )

    def resolve(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        constructed_site_id: pd.Series,
        duplicate_site_policy: str,
    ) -> DuplicateSiteResolutionResult:
        if duplicate_site_policy not in _SUPPORTED_DUPLICATE_SITE_POLICIES:
            raise PhosPyInputError(
                "dataset build request preprocessing_config contains an unsupported "
                "site_matrix.duplicate_site_policy"
            )

        if phospho.empty:
            empty_site_index = pd.Index([], name=_SITE_ID_COLUMN)
            empty_phospho = phospho.copy()
            empty_site_metadata = site_metadata.copy()
            empty_phospho.index = empty_site_index
            empty_site_metadata.index = empty_site_index.copy()
            return DuplicateSiteResolutionResult(
                phospho=empty_phospho,
                site_metadata=empty_site_metadata,
                dropped_row_count=0,
                duplicate_site_resolution=_empty_duplicate_site_resolution(),
                metadata_conflicts=_empty_metadata_conflicts(),
            )

        duplicate_mask = constructed_site_id.duplicated(keep=False)
        if not bool(duplicate_mask.any()):
            final_site_index = pd.Index(
                constructed_site_id.astype(str).tolist(), name=_SITE_ID_COLUMN
            )
            direct_phospho = phospho.copy()
            direct_site_metadata = site_metadata.copy()
            direct_phospho.index = final_site_index
            direct_site_metadata.index = final_site_index.copy()
            return DuplicateSiteResolutionResult(
                phospho=direct_phospho,
                site_metadata=direct_site_metadata,
                dropped_row_count=0,
                duplicate_site_resolution=_empty_duplicate_site_resolution(),
                metadata_conflicts=_empty_metadata_conflicts(),
            )

        dedupe_work = pd.DataFrame(
            {
                _SITE_ID_COLUMN: constructed_site_id.astype(str),
                "source_row_id": phospho.index.astype(str),
                "observed_values": phospho.notna().sum(axis=1),
                "mean_signal": phospho.mean(axis=1, skipna=True),
                "row_order": range(len(phospho.index)),
            },
            index=phospho.index,
        )
        duplicate_work = dedupe_work.loc[duplicate_mask].copy()
        duplicate_work.loc[:, "n_source_rows"] = (
            duplicate_work.groupby(_SITE_ID_COLUMN, sort=False)
            .size()
            .reindex(duplicate_work.loc[:, _SITE_ID_COLUMN])
            .to_numpy()
        )
        metadata_conflicts = self._metadata_conflict_detector.detect(
            site_metadata=site_metadata.loc[duplicate_mask],
            constructed_site_id=constructed_site_id.loc[duplicate_mask],
        )
        conflict_site_ids = set(metadata_conflicts.loc[:, "site_id"].astype(str))

        if duplicate_site_policy == DATASET_SITE_MATRIX_DUPLICATE_POLICY_ERROR:
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
                "site_matrix.duplicate_site_policy='error': "
                f"{preview}. Use a non-error duplicate policy to emit duplicate-site "
                "resolution and metadata-conflict diagnostics."
            )

        if duplicate_site_policy == DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST:
            selected_rows = (
                pd.DataFrame(
                    {_SITE_ID_COLUMN: constructed_site_id}, index=phospho.index
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
            duplicate_site_resolution = self._build_duplicate_site_resolution(
                duplicate_work=duplicate_work,
                site_metadata=site_metadata,
                selected_rows=selected_rows,
                duplicate_site_policy=duplicate_site_policy,
                retained_reason="selected first row by input order",
                dropped_reason=(
                    "dropped because another row was selected first by input order"
                ),
                conflict_site_ids=conflict_site_ids,
            )
            return DuplicateSiteResolutionResult(
                phospho=selected_phospho,
                site_metadata=selected_site_metadata,
                dropped_row_count=int(len(phospho.index) - len(selected_phospho.index)),
                duplicate_site_resolution=duplicate_site_resolution,
                metadata_conflicts=metadata_conflicts,
            )

        if (
            duplicate_site_policy
            == DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL
        ):
            value_columns = list(phospho.columns)
            dedupe_work = pd.DataFrame(
                {
                    _SITE_ID_COLUMN: constructed_site_id,
                    "observed_values": phospho.loc[:, value_columns]
                    .notna()
                    .sum(axis=1),
                    "mean_signal": phospho.loc[:, value_columns].mean(
                        axis=1, skipna=True
                    ),
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
            duplicate_site_resolution = self._build_duplicate_site_resolution(
                duplicate_work=duplicate_work,
                site_metadata=site_metadata,
                selected_rows=selected_rows,
                duplicate_site_policy=duplicate_site_policy,
                retained_reason=(
                    "selected row with highest mean signal under max_mean_signal criteria"
                ),
                dropped_reason=(
                    "dropped because another row ranked higher by max_mean_signal criteria"
                ),
                conflict_site_ids=conflict_site_ids,
            )
            return DuplicateSiteResolutionResult(
                phospho=selected_phospho,
                site_metadata=selected_site_metadata,
                dropped_row_count=int(len(phospho.index) - len(selected_phospho.index)),
                duplicate_site_resolution=duplicate_site_resolution,
                metadata_conflicts=metadata_conflicts,
            )

        if duplicate_site_policy in {
            DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
            DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN,
        }:
            metadata_columns = list(site_metadata.columns)
            grouped_metadata = (
                site_metadata.assign(
                    **{_SITE_ID_COLUMN: constructed_site_id.to_numpy()}
                )
                .groupby(_SITE_ID_COLUMN, sort=False)[metadata_columns]
                .first()
            )
            grouped_metadata = cast(pd.DataFrame, grouped_metadata)
            grouped_values = phospho.groupby(constructed_site_id, sort=False)
            if (
                duplicate_site_policy
                == DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN
            ):
                grouped_phospho = cast(pd.DataFrame, grouped_values.mean())
            else:
                grouped_phospho = cast(pd.DataFrame, grouped_values.median())
            grouped_phospho.index = pd.Index(
                grouped_phospho.index.astype(str), name=_SITE_ID_COLUMN
            )
            grouped_metadata.index = pd.Index(
                grouped_metadata.index.astype(str), name=_SITE_ID_COLUMN
            )
            duplicate_site_resolution = self._build_duplicate_site_resolution(
                duplicate_work=duplicate_work,
                site_metadata=site_metadata,
                selected_rows=duplicate_work.index,
                duplicate_site_policy=duplicate_site_policy,
                retained_reason=(
                    "contributed to site-level aggregate from duplicate source rows"
                ),
                dropped_reason=None,
                conflict_site_ids=conflict_site_ids,
                aggregated=True,
            )
            return DuplicateSiteResolutionResult(
                phospho=grouped_phospho,
                site_metadata=grouped_metadata,
                dropped_row_count=int(len(phospho.index) - len(grouped_phospho.index)),
                duplicate_site_resolution=duplicate_site_resolution,
                metadata_conflicts=metadata_conflicts,
            )

        raise RuntimeError("site-matrix duplicate policy dispatch fell through")

    @staticmethod
    def _build_duplicate_site_resolution(
        *,
        duplicate_work: pd.DataFrame,
        site_metadata: pd.DataFrame,
        selected_rows: pd.Index,
        duplicate_site_policy: str,
        retained_reason: str,
        dropped_reason: str | None,
        conflict_site_ids: set[str],
        aggregated: bool = False,
    ) -> pd.DataFrame:
        if duplicate_work.empty:
            return _empty_duplicate_site_resolution()

        selected_row_ids = set(selected_rows.astype(str).tolist())
        source_metadata = site_metadata.loc[duplicate_work.index]
        resolution = pd.DataFrame(
            {
                "site_id": duplicate_work.loc[:, _SITE_ID_COLUMN].astype(str).tolist(),
                "source_row_id": duplicate_work.loc[:, "source_row_id"]
                .astype(str)
                .tolist(),
                "retained": duplicate_work.index.astype(str).isin(selected_row_ids),
                "resolution_policy": duplicate_site_policy,
                "observed_values": duplicate_work.loc[:, "observed_values"].to_numpy(),
                "mean_signal": duplicate_work.loc[:, "mean_signal"].to_numpy(),
                "n_source_rows": duplicate_work.loc[:, "n_source_rows"].to_numpy(),
                "source_protein_id": _resolve_source_metadata_column(
                    source_metadata, "protein_id"
                ),
                "source_gene_symbol": _resolve_source_metadata_column(
                    source_metadata, "gene_symbol"
                ),
                "source_site": _resolve_source_metadata_column(source_metadata, "site"),
                "source_site_sequence": _resolve_source_metadata_column(
                    source_metadata, "site_sequence"
                ),
            },
            index=duplicate_work.index,
        )
        if aggregated:
            resolution.loc[:, "retained"] = True
            resolution.loc[:, "retained_reason"] = retained_reason
            resolution.loc[:, "dropped_reason"] = pd.NA
            resolution.loc[:, "n_aggregated_rows"] = resolution.loc[:, "n_source_rows"]
        else:
            retained_mask = resolution.loc[:, "retained"]
            resolution.loc[retained_mask, "retained_reason"] = retained_reason
            resolution.loc[~retained_mask, "retained_reason"] = pd.NA
            if dropped_reason is None:
                resolution.loc[~retained_mask, "dropped_reason"] = pd.NA
            else:
                resolution.loc[~retained_mask, "dropped_reason"] = dropped_reason
            resolution.loc[retained_mask, "dropped_reason"] = pd.NA
            resolution.loc[:, "n_aggregated_rows"] = pd.NA
        resolution.loc[:, "metadata_conflict_detected"] = (
            resolution.loc[:, "site_id"].astype(str).isin(conflict_site_ids)
        )
        return resolution.loc[:, list(DUPLICATE_SITE_RESOLUTION_COLUMNS)].reset_index(
            drop=True
        )


def _empty_duplicate_site_resolution() -> pd.DataFrame:
    return dataframe_from_duplicate_site_resolution_rows(())
