from __future__ import annotations

from enum import Enum
from typing import cast

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.models import DuplicateSiteResolutionResult
from phospy.science.datasets.preprocessing.policy_models import (
    SiteMatrixDuplicateSitePolicy,
)
from phospy.science.datasets.preprocessing.report_schema import (
    DUPLICATE_SITE_RESOLUTION_COLUMNS,
    dataframe_from_duplicate_site_resolution_rows,
)
from phospy.science.datasets.preprocessing.stages.site_matrix_components.metadata import (
    MetadataConflictDetector,
    _empty_metadata_conflicts,
    _resolve_source_metadata_column,
    resolve_aggregate_site_metadata,
)
from phospy.science.sites.identifiers import canonicalize_site_series
from phospy.science.sites.validation import (
    require_no_mixed_site_key_isoform_scope,
    require_site_key_series,
)

_SITE_KEY_COLUMN = "site_key"
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
        scientific_row_key: pd.Series | None = None,
        constructed_display_id: pd.Series | None = None,
        duplicate_site_policy: SiteMatrixDuplicateSitePolicy | str,
    ) -> DuplicateSiteResolutionResult:
        scientific_row_key = _require_valid_site_key_series(
            scientific_row_key,
            expected_index=phospho.index,
        )
        if constructed_display_id is None:
            raise PhosPyInputError(
                "duplicate-site resolution requires display_id metadata for diagnostics"
            )
        constructed_display_id = _require_display_id_series(
            constructed_display_id,
            expected_index=phospho.index,
        )
        output_index_name = _SITE_KEY_COLUMN
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
            empty_site_index = pd.Index([], name=output_index_name)
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

        duplicate_mask = scientific_row_key.duplicated(keep=False)
        if not bool(duplicate_mask.any()):
            final_site_index = pd.Index(
                scientific_row_key.astype(str).tolist(),
                name=output_index_name,
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
                _SITE_KEY_COLUMN: scientific_row_key.astype(str),
                "display_id": constructed_display_id.astype(str),
                "source_row_id": phospho.index.astype(str),
                "observed_values": phospho.notna().sum(axis=1),
                "mean_signal": phospho.mean(axis=1, skipna=True),
                "row_order": range(len(phospho.index)),
            },
            index=phospho.index,
        )
        duplicate_work = dedupe_work.loc[duplicate_mask].copy()
        duplicate_work.loc[:, "n_source_rows"] = (
            duplicate_work.groupby(_SITE_KEY_COLUMN, sort=False)
            .size()
            .reindex(duplicate_work.loc[:, _SITE_KEY_COLUMN])
            .to_numpy()
        )
        _raise_on_conflicting_duplicate_site_key_context(
            site_metadata=site_metadata,
            scientific_row_key=scientific_row_key,
            display_id=constructed_display_id,
        )
        metadata_conflicts = self._metadata_conflict_detector.detect(
            site_metadata=site_metadata.loc[duplicate_mask],
            scientific_row_key=scientific_row_key.loc[duplicate_mask],
            display_id=constructed_display_id.loc[duplicate_mask],
        )
        conflict_site_keys = set(metadata_conflicts.loc[:, "site_key"].astype(str))

        if resolved_policy is SiteMatrixDuplicateSitePolicy.ERROR:
            duplicate_site_keys = (
                scientific_row_key.loc[duplicate_mask]
                .astype(str)
                .drop_duplicates()
                .head(3)
            )
            preview = ", ".join(duplicate_site_keys.tolist())
            raise PhosPyInputError(
                "dataset build request preprocessing site-matrix construction found "
                "multiple input rows that resolved to the same analysis-ready "
                "site_key values (duplicate site_key values): "
                f"{preview}. Duplicate site_key rows are a scientific ambiguity "
                "because row retention, aggregation, or collapse changes the "
                "analysis-ready phosphosite evidence model. The default "
                "site_matrix.duplicate_site_policy='error' does not choose one "
                "source row silently. To resolve duplicates intentionally, configure "
                "DatasetSiteMatrixConfig(policy='build_from_metadata', "
                "duplicate_site_policy='<policy>') with one of the explicit "
                "non-error policies: 'max_mean_signal', 'first', "
                "'aggregate_mean', or 'aggregate_median'. Inspect "
                "dataset.preprocessing_report.duplicate_site_resolution and "
                "metadata_conflicts after using a non-error policy."
            )

        if resolved_policy is SiteMatrixDuplicateSitePolicy.FIRST:
            selected_rows = (
                pd.DataFrame(
                    {_SITE_KEY_COLUMN: scientific_row_key}, index=phospho.index
                )
                .drop_duplicates(_SITE_KEY_COLUMN, keep="first")
                .index
            )
            selected_phospho = phospho.loc[selected_rows].copy()
            selected_site_metadata = site_metadata.loc[selected_rows].copy()
            selected_site_keys = scientific_row_key.loc[selected_rows]
            final_site_index = pd.Index(
                selected_site_keys.astype(str).tolist(), name=output_index_name
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
                conflict_site_keys=conflict_site_keys,
                affected_sample_columns=_affected_sample_columns(phospho),
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
                    _SITE_KEY_COLUMN: scientific_row_key,
                    "display_id": constructed_display_id.astype(str),
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
                    [_SITE_KEY_COLUMN, "observed_values", "mean_signal", "row_order"],
                    ascending=[True, False, False, True],
                    kind="stable",
                    na_position="last",
                )
                .drop_duplicates(_SITE_KEY_COLUMN, keep="first")
                .index
            )
            selected_phospho = phospho.loc[selected_rows].copy()
            selected_site_metadata = site_metadata.loc[selected_rows].copy()
            selected_site_keys = scientific_row_key.loc[selected_rows]
            final_site_index = pd.Index(
                selected_site_keys.astype(str).tolist(), name=output_index_name
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
                conflict_site_keys=conflict_site_keys,
                affected_sample_columns=_affected_sample_columns(phospho),
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
                scientific_row_key=scientific_row_key,
                display_id=constructed_display_id,
                metadata_conflicts=metadata_conflicts,
            )
            grouped_metadata.index = pd.Index(
                grouped_metadata.index.astype(str),
                name=output_index_name,
            )
            grouped_values = phospho.groupby(scientific_row_key, sort=False)
            if resolved_policy is SiteMatrixDuplicateSitePolicy.AGGREGATE_MEAN:
                grouped_phospho = cast(pd.DataFrame, grouped_values.mean())
            else:
                grouped_phospho = cast(pd.DataFrame, grouped_values.median())
            grouped_phospho.index = pd.Index(
                grouped_phospho.index.astype(str), name=output_index_name
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
                conflict_site_keys=conflict_site_keys,
                affected_sample_columns=_affected_sample_columns(phospho),
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
        conflict_site_keys: set[str],
        affected_sample_columns: tuple[str, ...],
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
                "site_key": duplicate_work.loc[:, _SITE_KEY_COLUMN]
                .astype(str)
                .tolist(),
                "display_id": duplicate_work.loc[:, "display_id"].astype(str).tolist(),
                "site_id": duplicate_work.loc[:, "display_id"].astype(str).tolist(),
                "source_row_id": duplicate_work.loc[:, "source_row_id"]
                .astype(str)
                .tolist(),
                "retained": duplicate_work.index.astype(str).isin(selected_row_ids),
                "resolution_policy": duplicate_site_policy.value,
                "aggregation_method": duplicate_site_policy.value,
                "missing_value_policy": missing_value_policy,
                "metadata_resolution_policy": metadata_resolution_policy,
                "affected_sample_columns": [affected_sample_columns]
                * int(len(duplicate_work.index)),
                "observed_values": duplicate_work.loc[:, "observed_values"].to_numpy(),
                "mean_signal": duplicate_work.loc[:, "mean_signal"].to_numpy(),
                "n_source_rows": duplicate_work.loc[:, "n_source_rows"].to_numpy(),
                "source_organism": _resolve_source_metadata_column(
                    source_metadata, "organism"
                ),
                "source_protein_namespace": _resolve_source_metadata_column(
                    source_metadata, "protein_namespace"
                ),
                "source_protein_identifier": _resolve_source_metadata_column(
                    source_metadata, "protein_identifier"
                ),
                "source_protein_id": _resolve_source_metadata_column(
                    source_metadata, "protein_id"
                ),
                "source_protein_accession": _resolve_source_metadata_column(
                    source_metadata, "protein_accession"
                ),
                "source_isoform_id": _resolve_source_metadata_column(
                    source_metadata, "isoform_id"
                ),
                "source_gene_symbol": _resolve_source_metadata_column(
                    source_metadata, "gene_symbol"
                ),
                "source_site": _resolve_source_metadata_column(source_metadata, "site"),
                "source_residue": _resolve_source_metadata_column(
                    source_metadata, "residue"
                ),
                "source_position": _resolve_source_metadata_column(
                    source_metadata, "position"
                ),
                "source_site_position": _resolve_source_metadata_column(
                    source_metadata, "site_position"
                ),
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
            resolution.loc[:, "site_key"].astype(str).isin(conflict_site_keys)
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
                else int(duplicate_site_resolution.loc[:, "site_key"].nunique())
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
            "affected_sample_columns": _affected_sample_columns(input_phospho),
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


def _require_valid_site_key_series(
    scientific_row_key: pd.Series | None,
    *,
    expected_index: pd.Index,
) -> pd.Series:
    if scientific_row_key is None:
        raise PhosPyInputError(
            "duplicate-site resolution requires site_key row identity; "
            "display_id and constructed GENE;SITE; labels are diagnostic metadata "
            "only"
        )
    if scientific_row_key.name != _SITE_KEY_COLUMN:
        raise PhosPyInputError(
            "duplicate-site resolution requires scientific_row_key.name='site_key'"
        )
    if not scientific_row_key.index.equals(expected_index):
        raise PhosPyInputError(
            "duplicate-site resolution requires site_key values aligned to phospho rows"
        )
    site_keys = require_site_key_series(
        scientific_row_key.astype("object"),
        field_name="duplicate-site resolution site_key",
        error_type=PhosPyInputError,
    )
    require_no_mixed_site_key_isoform_scope(
        site_keys=site_keys,
        field_name="duplicate-site resolution site_key",
        error_type=PhosPyInputError,
    )
    return pd.Series(
        site_keys.astype(str).tolist(),
        index=expected_index.copy(),
        name=_SITE_KEY_COLUMN,
        dtype="object",
    )


def _require_display_id_series(
    display_id: pd.Series,
    *,
    expected_index: pd.Index,
) -> pd.Series:
    if not display_id.index.equals(expected_index):
        raise PhosPyInputError(
            "duplicate-site resolution requires display_id metadata aligned to "
            "phospho rows"
        )
    display_ids = canonicalize_site_series(
        display_id.astype("object"),
        field_name="duplicate-site resolution display_id",
        error_type=PhosPyInputError,
    )
    return pd.Series(
        display_ids.astype(str).tolist(),
        index=expected_index.copy(),
        name="display_id",
        dtype="object",
    )


def _affected_sample_columns(phospho: pd.DataFrame) -> tuple[str, ...]:
    return tuple(str(column) for column in phospho.columns.tolist())


def _raise_on_conflicting_duplicate_site_key_context(
    *,
    site_metadata: pd.DataFrame,
    scientific_row_key: pd.Series,
    display_id: pd.Series,
) -> None:
    duplicate_mask = scientific_row_key.duplicated(keep=False)
    if not bool(duplicate_mask.any()):
        return
    duplicate_metadata = site_metadata.loc[duplicate_mask]
    duplicate_site_key = scientific_row_key.loc[duplicate_mask].astype(str)
    duplicate_display = display_id.loc[duplicate_mask].astype(str)
    context_columns = (
        "organism",
        "protein_namespace",
        "protein_identifier",
        "protein_id",
        "protein_accession",
        "isoform_id",
        "residue",
        "position",
        "site_position",
    )
    for column_name in context_columns:
        if column_name not in duplicate_metadata.columns:
            continue
        column_values = duplicate_metadata.loc[:, column_name]
        normalized = column_values.astype("string").str.strip().fillna("").astype(str)
        grouped = normalized.groupby(duplicate_site_key, sort=False)
        for site_key, values in grouped:
            non_empty = [value for value in values.tolist() if value != ""]
            if len(set(non_empty)) > 1:
                display_values = tuple(
                    dict.fromkeys(
                        duplicate_display.loc[duplicate_site_key == site_key]
                        .astype(str)
                        .tolist()
                    )
                )
                raise PhosPyInputError(
                    "dataset build request preprocessing site-matrix construction "
                    "found conflicting protein context for duplicate site_key rows; "
                    f"site_key={site_key!r}, display_id_values={display_values!r}, "
                    f"column={column_name!r}"
                )
