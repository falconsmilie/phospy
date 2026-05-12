from __future__ import annotations

from enum import Enum
from typing import cast

import pandas as pd

from phospy.datasets.preprocessing.models import DuplicateSiteResolutionResult
from phospy.datasets.preprocessing.report_schema import (
    DUPLICATE_SITE_RESOLUTION_COLUMNS,
    dataframe_from_duplicate_site_resolution_rows,
)
from phospy.datasets.preprocessing.stages.site_matrix_components.metadata import (
    MetadataConflictDetector,
    _empty_metadata_conflicts,
    _resolve_source_metadata_column,
    resolve_aggregate_site_metadata,
)
from phospy.errors.input import PhosPyInputError
from phospy.policy_models import SiteMatrixDuplicateSitePolicy
from phospy.sites.identity import validate_no_conflicting_identity_collisions

_SITE_ID_COLUMN = "site_id"
_AGGREGATE_DUPLICATE_SITE_POLICIES = {
    SiteMatrixDuplicateSitePolicy.AGGREGATE_MEAN,
    SiteMatrixDuplicateSitePolicy.AGGREGATE_MEDIAN,
}
_SUPPORTED_DUPLICATE_SITE_POLICIES = {
    SiteMatrixDuplicateSitePolicy.MAX_MEAN_SIGNAL,
    SiteMatrixDuplicateSitePolicy.FIRST,
    *_AGGREGATE_DUPLICATE_SITE_POLICIES,
    SiteMatrixDuplicateSitePolicy.ERROR,
}


class DuplicateMissingValuePolicy(str, Enum):
    SKIP_MISSING_VALUES = "skip_missing_values"
    REQUIRE_COMPLETE_OBSERVATIONS = "require_complete_observations"
    PROPAGATE_MISSING_IF_ANY_SOURCE_MISSING = "propagate_missing_if_any_source_missing"


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
        duplicate_site_policy: SiteMatrixDuplicateSitePolicy | str,
    ) -> DuplicateSiteResolutionResult:
        resolved_policy = SiteMatrixDuplicateSitePolicy.parse(
            duplicate_site_policy,
            field_name="site_matrix.duplicate_site_policy",
        )
        if resolved_policy not in _SUPPORTED_DUPLICATE_SITE_POLICIES:
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
                duplicate_aggregation_diagnostics=(
                    self._build_duplicate_aggregation_diagnostics(
                        input_phospho=phospho,
                        output_phospho=empty_phospho,
                        duplicate_site_resolution=_empty_duplicate_site_resolution(),
                        duplicate_site_policy=resolved_policy,
                    )
                ),
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
            duplicate_site_resolution = _empty_duplicate_site_resolution()
            return DuplicateSiteResolutionResult(
                phospho=direct_phospho,
                site_metadata=direct_site_metadata,
                dropped_row_count=0,
                duplicate_site_resolution=duplicate_site_resolution,
                metadata_conflicts=_empty_metadata_conflicts(),
                duplicate_aggregation_diagnostics=(
                    self._build_duplicate_aggregation_diagnostics(
                        input_phospho=phospho,
                        output_phospho=direct_phospho,
                        duplicate_site_resolution=duplicate_site_resolution,
                        duplicate_site_policy=resolved_policy,
                    )
                ),
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
        validate_no_conflicting_identity_collisions(
            site_metadata=site_metadata.loc[duplicate_mask],
            display_ids=constructed_site_id.loc[duplicate_mask],
            field_name=(
                "dataset build request preprocessing site-matrix duplicate-site "
                "identity validation"
            ),
            error_type=PhosPyInputError,
        )
        metadata_conflicts = self._metadata_conflict_detector.detect(
            site_metadata=site_metadata.loc[duplicate_mask],
            constructed_site_id=constructed_site_id.loc[duplicate_mask],
        )
        conflict_site_ids = set(metadata_conflicts.loc[:, "site_id"].astype(str))

        if resolved_policy is SiteMatrixDuplicateSitePolicy.ERROR:
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

        if resolved_policy is SiteMatrixDuplicateSitePolicy.FIRST:
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
                duplicate_site_policy=resolved_policy,
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
                duplicate_aggregation_diagnostics=(
                    self._build_duplicate_aggregation_diagnostics(
                        input_phospho=phospho,
                        output_phospho=selected_phospho,
                        duplicate_site_resolution=duplicate_site_resolution,
                        duplicate_site_policy=resolved_policy,
                    )
                ),
            )

        if resolved_policy is SiteMatrixDuplicateSitePolicy.MAX_MEAN_SIGNAL:
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
                duplicate_site_policy=resolved_policy,
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
                duplicate_aggregation_diagnostics=(
                    self._build_duplicate_aggregation_diagnostics(
                        input_phospho=phospho,
                        output_phospho=selected_phospho,
                        duplicate_site_resolution=duplicate_site_resolution,
                        duplicate_site_policy=resolved_policy,
                    )
                ),
            )

        if resolved_policy in _AGGREGATE_DUPLICATE_SITE_POLICIES:
            grouped_metadata = resolve_aggregate_site_metadata(
                site_metadata=site_metadata,
                constructed_site_id=constructed_site_id,
                metadata_conflicts=metadata_conflicts,
            )
            grouped_values = phospho.groupby(constructed_site_id, sort=False)
            if resolved_policy is SiteMatrixDuplicateSitePolicy.AGGREGATE_MEAN:
                grouped_phospho = cast(pd.DataFrame, grouped_values.mean())
            else:
                grouped_phospho = cast(pd.DataFrame, grouped_values.median())
            grouped_phospho.index = pd.Index(
                grouped_phospho.index.astype(str), name=_SITE_ID_COLUMN
            )
            duplicate_site_resolution = self._build_duplicate_site_resolution(
                duplicate_work=duplicate_work,
                site_metadata=site_metadata,
                selected_rows=duplicate_work.index,
                duplicate_site_policy=resolved_policy,
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
                duplicate_aggregation_diagnostics=(
                    self._build_duplicate_aggregation_diagnostics(
                        input_phospho=phospho,
                        output_phospho=grouped_phospho,
                        duplicate_site_resolution=duplicate_site_resolution,
                        duplicate_site_policy=resolved_policy,
                    )
                ),
            )

        raise RuntimeError("site-matrix duplicate policy dispatch fell through")

    @staticmethod
    def _build_duplicate_site_resolution(
        *,
        duplicate_work: pd.DataFrame,
        site_metadata: pd.DataFrame,
        selected_rows: pd.Index,
        duplicate_site_policy: SiteMatrixDuplicateSitePolicy,
        retained_reason: str,
        dropped_reason: str | None,
        conflict_site_ids: set[str],
        aggregated: bool = False,
    ) -> pd.DataFrame:
        if duplicate_work.empty:
            return _empty_duplicate_site_resolution()

        selected_row_ids = set(selected_rows.astype(str).tolist())
        source_metadata = site_metadata.loc[duplicate_work.index]
        is_aggregate_policy = (
            duplicate_site_policy in _AGGREGATE_DUPLICATE_SITE_POLICIES
        )
        missing_value_policy: object = (
            DuplicateMissingValuePolicy.SKIP_MISSING_VALUES.value
            if is_aggregate_policy
            else pd.NA
        )
        metadata_resolution_policy: object = _metadata_resolution_policy_text(
            duplicate_site_policy=duplicate_site_policy,
            is_aggregate_policy=is_aggregate_policy,
            for_row_resolution=True,
        )
        resolution = pd.DataFrame(
            {
                "site_id": duplicate_work.loc[:, _SITE_ID_COLUMN].astype(str).tolist(),
                "source_row_id": duplicate_work.loc[:, "source_row_id"]
                .astype(str)
                .tolist(),
                "retained": duplicate_work.index.astype(str).isin(selected_row_ids),
                "resolution_policy": duplicate_site_policy.value,
                "aggregation_method": duplicate_site_policy.value,
                "missing_value_policy": missing_value_policy,
                "metadata_resolution_policy": metadata_resolution_policy,
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

    @staticmethod
    def _build_duplicate_aggregation_diagnostics(
        *,
        input_phospho: pd.DataFrame,
        output_phospho: pd.DataFrame,
        duplicate_site_resolution: pd.DataFrame,
        duplicate_site_policy: SiteMatrixDuplicateSitePolicy,
    ) -> dict[str, object]:
        missing_cells_before = _missing_cell_count(input_phospho)
        missing_cells_after = _missing_cell_count(output_phospho)
        is_aggregate_policy = (
            duplicate_site_policy in _AGGREGATE_DUPLICATE_SITE_POLICIES
        )
        missing_value_policy = (
            DuplicateMissingValuePolicy.SKIP_MISSING_VALUES.value
            if is_aggregate_policy
            else "not_applicable_row_selection"
        )
        metadata_resolution_policy = _metadata_resolution_policy_text(
            duplicate_site_policy=duplicate_site_policy,
            is_aggregate_policy=is_aggregate_policy,
            for_row_resolution=False,
        )
        return {
            "aggregation_method": duplicate_site_policy.value,
            "missing_value_policy": missing_value_policy,
            "duplicate_group_count": (
                0
                if duplicate_site_resolution.empty
                else int(duplicate_site_resolution.loc[:, "site_id"].nunique())
            ),
            "rows_collapsed_count": max(
                int(len(input_phospho.index) - len(output_phospho.index)), 0
            ),
            "missing_cells_before_aggregation": missing_cells_before,
            "missing_cells_after_aggregation": missing_cells_after,
            "aggregation_reduced_missingness": (
                missing_cells_after < missing_cells_before
            ),
            "metadata_resolution_policy": metadata_resolution_policy,
        }


def _empty_duplicate_site_resolution() -> pd.DataFrame:
    return dataframe_from_duplicate_site_resolution_rows(())


def _missing_cell_count(frame: pd.DataFrame) -> int:
    return int(frame.isna().to_numpy().sum())


def _metadata_resolution_policy_text(
    *,
    duplicate_site_policy: SiteMatrixDuplicateSitePolicy,
    is_aggregate_policy: bool,
    for_row_resolution: bool,
) -> str | object:
    if is_aggregate_policy:
        return "first_non_missing_value_per_site_then_set_conflicting_fields_to_missing"
    if duplicate_site_policy is SiteMatrixDuplicateSitePolicy.FIRST:
        return "retain_earliest_input_row_per_site"
    if duplicate_site_policy is SiteMatrixDuplicateSitePolicy.MAX_MEAN_SIGNAL:
        return "retain_row_ranked_by_observed_values_then_mean_signal_then_input_order"
    if duplicate_site_policy is SiteMatrixDuplicateSitePolicy.ERROR:
        return "error_on_duplicate_sites"
    return pd.NA if for_row_resolution else "not_applicable"
