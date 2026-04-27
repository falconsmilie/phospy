"""Focused components for site-matrix duplicate/conflict/audit handling."""

from __future__ import annotations

import pandas as pd

from phospy.api.configs import (
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_ERROR,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL,
)
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DuplicateSiteResolutionResult,
)
from phospy.datasets.preprocessing.report_schema import (
    DUPLICATE_SITE_RESOLUTION_COLUMNS,
    METADATA_CONFLICT_COLUMNS,
    PreprocessingRowAuditRow,
    dataframe_from_duplicate_site_resolution_rows,
    dataframe_from_metadata_conflict_rows,
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
_DEFAULT_METADATA_CONFLICT_FIELDS = (
    "protein_id",
    "gene_symbol",
    "site",
    "site_sequence",
    "residue",
    "position",
)


class MetadataConflictDetector:
    """Detect metadata conflicts for rows mapping to the same constructed site id."""

    def __init__(
        self,
        *,
        conflict_fields: tuple[str, ...] = _DEFAULT_METADATA_CONFLICT_FIELDS,
    ) -> None:
        self._conflict_fields = conflict_fields

    def detect(
        self,
        *,
        site_metadata: pd.DataFrame,
        constructed_site_id: pd.Series,
    ) -> pd.DataFrame:
        if site_metadata.empty:
            return _empty_metadata_conflicts()

        duplicate_groups = site_metadata.assign(
            **{
                _SITE_ID_COLUMN: constructed_site_id.astype(str),
                "source_row_id": site_metadata.index.astype(str),
            }
        )
        records: list[dict[str, object]] = []
        for site_id, group in duplicate_groups.groupby(_SITE_ID_COLUMN, sort=False):
            source_row_ids = tuple(group.loc[:, "source_row_id"].astype(str).tolist())
            for field in self._conflict_fields:
                if field not in group.columns:
                    continue
                observed_values = group.loc[:, field]
                distinct_values = []
                for value in observed_values.tolist():
                    normalized = _normalize_metadata_value(value)
                    if normalized not in distinct_values:
                        distinct_values.append(normalized)
                if len(distinct_values) <= 1:
                    continue
                records.append(
                    {
                        "site_id": str(site_id),
                        "field": field,
                        "values": tuple(str(value) for value in distinct_values),
                        "n_distinct_values": len(distinct_values),
                        "source_row_ids": source_row_ids,
                    }
                )
        if not records:
            return _empty_metadata_conflicts()
        return pd.DataFrame.from_records(records, columns=METADATA_CONFLICT_COLUMNS)


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
            grouped_values = phospho.groupby(constructed_site_id, sort=False)
            if (
                duplicate_site_policy
                == DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN
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


class SiteMatrixRowAuditBuilder:
    """Build row-audit records for site-matrix sequence/missing/duplicate decisions."""

    def build(
        self,
        *,
        dropped_missing_sequence_rows: tuple[tuple[str, str], ...],
        dropped_incomplete_rows: tuple[tuple[str, str, int], ...],
        duplicate_site_resolution: pd.DataFrame,
        site_matrix_policy: str,
        site_matrix_missing_data_policy: str,
        site_matrix_duplicate_site_policy: str,
        required_observed_count: int,
    ) -> list[PreprocessingRowAuditRow]:
        records: list[PreprocessingRowAuditRow] = []
        for source_row_id, site_id in dropped_missing_sequence_rows:
            records.append(
                PreprocessingRowAuditRow(
                    stage=DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
                    action="dropped",
                    reason=(
                        "dropped because site_metadata.site_sequence is missing or blank"
                    ),
                    source_row_id=source_row_id,
                    site_id=site_id,
                    retained=False,
                    retained_row_id=pd.NA,
                    source_rows=(source_row_id,),
                    retained_row=pd.NA,
                    parameter_snapshot={
                        "site_matrix_policy": site_matrix_policy,
                        "site_matrix_missing_data_policy": site_matrix_missing_data_policy,
                        "site_matrix_duplicate_site_policy": (
                            site_matrix_duplicate_site_policy
                        ),
                    },
                )
            )
        for source_row_id, site_id, observed_value_count in dropped_incomplete_rows:
            records.append(
                PreprocessingRowAuditRow(
                    stage=DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
                    action="dropped",
                    reason="dropped by site_matrix missing-data policy",
                    source_row_id=source_row_id,
                    site_id=site_id,
                    retained=False,
                    retained_row_id=pd.NA,
                    source_rows=(source_row_id,),
                    retained_row=pd.NA,
                    parameter_snapshot={
                        "site_matrix_policy": site_matrix_policy,
                        "site_matrix_missing_data_policy": site_matrix_missing_data_policy,
                        "site_matrix_duplicate_site_policy": (
                            site_matrix_duplicate_site_policy
                        ),
                        "observed_values": int(observed_value_count),
                        "required_observed_count": int(required_observed_count),
                    },
                )
            )
        records.extend(
            self._build_duplicate_site_row_audit_records(
                duplicate_site_resolution=duplicate_site_resolution,
                site_matrix_policy=site_matrix_policy,
                site_matrix_duplicate_site_policy=site_matrix_duplicate_site_policy,
            )
        )
        return records

    @staticmethod
    def _build_duplicate_site_row_audit_records(
        *,
        duplicate_site_resolution: pd.DataFrame,
        site_matrix_policy: str,
        site_matrix_duplicate_site_policy: str,
    ) -> list[PreprocessingRowAuditRow]:
        if duplicate_site_resolution.empty:
            return []

        source_rows_by_site = (
            duplicate_site_resolution.groupby("site_id", sort=False)["source_row_id"]
            .apply(lambda values: tuple(str(value) for value in values.tolist()))
            .to_dict()
        )
        retained_row_by_site = (
            duplicate_site_resolution.loc[duplicate_site_resolution.loc[:, "retained"]]
            .groupby("site_id", sort=False)["source_row_id"]
            .first()
            .astype(str)
            .to_dict()
        )
        aggregated = site_matrix_duplicate_site_policy in {
            DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
            DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN,
        }
        records: list[PreprocessingRowAuditRow] = []
        for row in duplicate_site_resolution.itertuples(index=False):
            site_id = str(row.site_id)
            source_row_id = str(row.source_row_id)
            source_rows = source_rows_by_site.get(site_id, (source_row_id,))
            if aggregated:
                action = "aggregated"
                retained = True
                retained_row_id: object = site_id
                retained_row: object = site_id
                reason = _optional_text(
                    row.retained_reason,
                    fallback=(
                        "contributed to site-level aggregate from duplicate source rows"
                    ),
                )
            else:
                retained = bool(row.retained)
                action = "retained" if retained else "collapsed"
                retained_row_id = retained_row_by_site.get(site_id, pd.NA)
                retained_row = (
                    retained_row_id if retained_row_id is not pd.NA else pd.NA
                )
                reason = _optional_text(
                    row.retained_reason if retained else row.dropped_reason,
                    fallback="duplicate-site resolution decision",
                )
            records.append(
                PreprocessingRowAuditRow(
                    stage=DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
                    action=action,
                    reason=reason,
                    source_row_id=source_row_id,
                    site_id=site_id,
                    retained=retained,
                    retained_row_id=retained_row_id,
                    source_rows=source_rows,
                    retained_row=retained_row,
                    parameter_snapshot={
                        "site_matrix_policy": site_matrix_policy,
                        "duplicate_site_policy": site_matrix_duplicate_site_policy,
                        "site_matrix_duplicate_site_policy": (
                            site_matrix_duplicate_site_policy
                        ),
                        "observed_values": _optional_int(row.observed_values),
                        "mean_signal": _optional_float(row.mean_signal),
                        "metadata_conflict_detected": bool(
                            _optional_bool(row.metadata_conflict_detected)
                        ),
                    },
                )
            )
        return records


def _empty_duplicate_site_resolution() -> pd.DataFrame:
    return dataframe_from_duplicate_site_resolution_rows(())


def _empty_metadata_conflicts() -> pd.DataFrame:
    return dataframe_from_metadata_conflict_rows(())


def _resolve_source_metadata_column(
    site_metadata: pd.DataFrame, column: str
) -> pd.Series:
    if column not in site_metadata.columns:
        return pd.Series(pd.NA, index=site_metadata.index)
    return site_metadata.loc[:, column].copy()


def _normalize_metadata_value(value: object) -> str:
    if pd.isna(value):
        return "<NA>"
    if isinstance(value, str):
        trimmed = value.strip()
        return trimmed if trimmed else "<EMPTY>"
    return str(value)


def _optional_text(value: object, *, fallback: str) -> str:
    if _is_missing_scalar(value):
        return fallback
    text = str(value).strip()
    return text if text else fallback


def _optional_int(value: object) -> int | None:
    if _is_missing_scalar(value):
        return None
    return int(value)


def _optional_float(value: object) -> float | None:
    if _is_missing_scalar(value):
        return None
    return float(value)


def _optional_bool(value: object) -> bool | None:
    if _is_missing_scalar(value):
        return None
    return bool(value)


def _is_missing_scalar(value: object) -> bool:
    if value is pd.NA:
        return True
    if isinstance(value, (list, tuple, dict)):
        return False
    try:
        return bool(pd.isna(value))
    except TypeError:
        return False


__all__ = [
    "DuplicateSiteResolver",
    "MetadataConflictDetector",
    "SiteMatrixRowAuditBuilder",
]
