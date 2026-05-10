from __future__ import annotations

from dataclasses import dataclass
from typing import cast

import pandas as pd

from phospy.datasets.preprocessing.report_schema import (
    METADATA_CONFLICT_COLUMNS,
    dataframe_from_metadata_conflict_rows,
)
from phospy.policy_models import (
    SiteMatrixDuplicateSitePolicy,
    SiteMatrixMissingDataPolicy,
    SiteMatrixPolicy,
)

_DEFAULT_METADATA_CONFLICT_FIELDS = (
    "protein_id",
    "gene_symbol",
    "site",
    "site_sequence",
    "residue",
    "position",
)
_DEFAULT_ROW_DROP_STATS_ATTR = "site_matrix_row_drop_stats"
_DEFAULT_SITE_MATRIX_POLICY_ATTR = "site_matrix_policy"
_DEFAULT_SITE_MATRIX_PROVENANCE_ATTR = "site_matrix_provenance"
_SITE_ID_COLUMN = "site_id"


@dataclass(frozen=True, slots=True)
class SiteMatrixProvenanceResult:
    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    row_drop_stats: dict[str, int | str]
    site_matrix_provenance: dict[str, object]
    diagnostics: dict[str, object]


class SiteMatrixProvenanceBuilder:
    """Attach final site-matrix attrs and build diagnostics/provenance payloads."""

    def __init__(
        self,
        *,
        row_drop_stats_attr: str = _DEFAULT_ROW_DROP_STATS_ATTR,
        site_matrix_policy_attr: str = _DEFAULT_SITE_MATRIX_POLICY_ATTR,
        site_matrix_provenance_attr: str = _DEFAULT_SITE_MATRIX_PROVENANCE_ATTR,
    ) -> None:
        self._row_drop_stats_attr = row_drop_stats_attr
        self._site_matrix_policy_attr = site_matrix_policy_attr
        self._site_matrix_provenance_attr = site_matrix_provenance_attr

    def build(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        input_rows: int,
        dropped_missing_sequence: int,
        dropped_incomplete_values: int,
        missing_data_policy: SiteMatrixMissingDataPolicy | str,
        required_observed_count: int,
        deduplicated_site_rows: int,
        duplicate_site_policy: SiteMatrixDuplicateSitePolicy | str,
        site_matrix_policy: SiteMatrixPolicy | str,
        dropped_missing_sequence_row_ids: tuple[str, ...],
        dropped_incomplete_row_ids: tuple[str, ...],
        dropped_row_ids: tuple[str, ...],
        duplicate_site_resolution: pd.DataFrame | None,
        duplicate_aggregation_diagnostics: dict[str, object] | None,
    ) -> SiteMatrixProvenanceResult:
        resolved_missing_data_policy = SiteMatrixMissingDataPolicy.parse(
            missing_data_policy,
            field_name="site_matrix.missing_data_policy",
        )
        resolved_duplicate_site_policy = SiteMatrixDuplicateSitePolicy.parse(
            duplicate_site_policy,
            field_name="site_matrix.duplicate_site_policy",
        )
        resolved_site_matrix_policy = SiteMatrixPolicy.parse(
            site_matrix_policy,
            field_name="site_matrix.policy",
        )
        row_drop_stats = {
            "input_rows": int(input_rows),
            "dropped_missing_sequence": dropped_missing_sequence,
            "dropped_incomplete_values": dropped_incomplete_values,
            "missing_data_policy": resolved_missing_data_policy.value,
            "required_observed_count": required_observed_count,
            "deduplicated_site_rows": deduplicated_site_rows,
            "duplicate_site_policy": resolved_duplicate_site_policy.value,
            "retained_rows": int(len(phospho.index)),
        }
        final_phospho = phospho.copy()
        final_site_metadata = site_metadata.copy()
        final_phospho.attrs[self._row_drop_stats_attr] = row_drop_stats.copy()
        final_site_metadata.attrs[self._row_drop_stats_attr] = row_drop_stats.copy()
        final_phospho.attrs[self._site_matrix_policy_attr] = (
            resolved_site_matrix_policy.value
        )
        final_site_metadata.attrs[self._site_matrix_policy_attr] = (
            resolved_site_matrix_policy.value
        )

        site_matrix_provenance = {
            "dropped_missing_sequence_row_ids": dropped_missing_sequence_row_ids,
            "dropped_incomplete_row_ids": dropped_incomplete_row_ids,
            "dropped_row_ids": dropped_row_ids,
            "duplicate_site_policy": resolved_duplicate_site_policy.value,
            "missing_data_policy": resolved_missing_data_policy.value,
            "required_observed_count": required_observed_count,
            "final_constructed_site_ids": tuple(
                str(site_id) for site_id in final_phospho.index.tolist()
            ),
        }
        if duplicate_aggregation_diagnostics is not None:
            site_matrix_provenance["duplicate_aggregation"] = dict(
                duplicate_aggregation_diagnostics
            )
        final_phospho.attrs[self._site_matrix_provenance_attr] = site_matrix_provenance
        final_site_metadata.attrs[self._site_matrix_provenance_attr] = (
            site_matrix_provenance.copy()
        )

        diagnostics = dict(site_matrix_provenance)
        diagnostics["final_constructed_site_ids"] = [
            str(site_id) for site_id in final_phospho.index.tolist()
        ]
        if duplicate_aggregation_diagnostics is not None:
            diagnostics["duplicate_aggregation"] = dict(
                duplicate_aggregation_diagnostics
            )
        if duplicate_site_resolution is not None:
            diagnostics["duplicate_site_decisions"] = _records_from_frame(
                duplicate_site_resolution
            )
        return SiteMatrixProvenanceResult(
            phospho=final_phospho,
            site_metadata=final_site_metadata,
            row_drop_stats=row_drop_stats,
            site_matrix_provenance=site_matrix_provenance,
            diagnostics=diagnostics,
        )


class MetadataConflictDetector:
    """Detect metadata conflicts for rows mapping to the same constructed site id."""

    def __init__(
        self,
        *,
        conflict_fields: tuple[str, ...] | None = None,
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
        conflict_fields = (
            tuple(site_metadata.columns)
            if self._conflict_fields is None
            else self._conflict_fields
        )
        records: list[dict[str, object]] = []
        for site_id, group in duplicate_groups.groupby(_SITE_ID_COLUMN, sort=False):
            source_row_ids = tuple(group.loc[:, "source_row_id"].astype(str).tolist())
            for field in conflict_fields:
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


def resolve_aggregate_site_metadata(
    *,
    site_metadata: pd.DataFrame,
    constructed_site_id: pd.Series,
    metadata_conflicts: pd.DataFrame,
) -> pd.DataFrame:
    if site_metadata.empty:
        return site_metadata.copy()

    metadata_columns = list(site_metadata.columns)
    grouped_metadata = cast(
        pd.DataFrame,
        site_metadata.assign(**{_SITE_ID_COLUMN: constructed_site_id.to_numpy()})
        .groupby(_SITE_ID_COLUMN, sort=False)[metadata_columns]
        .first(),
    )
    grouped_metadata.index = pd.Index(
        grouped_metadata.index.astype(str), name=_SITE_ID_COLUMN
    )

    if metadata_conflicts.empty:
        return cast(pd.DataFrame, grouped_metadata)

    conflict_records = metadata_conflicts.loc[:, ["site_id", "field"]].drop_duplicates()
    for conflict in conflict_records.to_dict(orient="records"):
        site_id = str(conflict["site_id"])
        field = str(conflict["field"])
        if site_id in grouped_metadata.index and field in grouped_metadata.columns:
            grouped_metadata.at[site_id, field] = pd.NA

    return cast(pd.DataFrame, grouped_metadata)


def _empty_metadata_conflicts() -> pd.DataFrame:
    return dataframe_from_metadata_conflict_rows(())


def _resolve_source_metadata_column(
    site_metadata: pd.DataFrame, column: str
) -> pd.Series:
    if column not in site_metadata.columns:
        return pd.Series(pd.NA, index=site_metadata.index)
    return site_metadata.loc[:, column].copy()


def _normalize_metadata_value(value: object) -> str:
    if bool(pd.Series((value,), dtype="object").isna().iat[0]):
        return "<NA>"
    if isinstance(value, str):
        trimmed = value.strip()
        return trimmed if trimmed else "<EMPTY>"
    return str(value)


def _is_missing_scalar(value: object) -> bool:
    if value is pd.NA:
        return True
    if isinstance(value, (list, tuple, dict)):
        return False
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _records_from_frame(frame: pd.DataFrame) -> list[dict[str, object]]:
    if frame.empty:
        return []
    records: list[dict[str, object]] = []
    for raw_record in frame.to_dict(orient="records"):
        record: dict[str, object] = {}
        for key, value in raw_record.items():
            if isinstance(value, tuple):
                record[str(key)] = [item for item in value]
            elif _is_missing_scalar(value):
                record[str(key)] = None
            else:
                record[str(key)] = value
        records.append(record)
    return records
