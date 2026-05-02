from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.datasets.preprocessing.report_schema import (
    METADATA_CONFLICT_COLUMNS,
    dataframe_from_metadata_conflict_rows,
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
        missing_data_policy: str,
        required_observed_count: int,
        deduplicated_site_rows: int,
        duplicate_site_policy: str,
        site_matrix_policy: str,
        dropped_missing_sequence_row_ids: tuple[str, ...],
        dropped_incomplete_row_ids: tuple[str, ...],
        dropped_row_ids: tuple[str, ...],
        duplicate_site_resolution: pd.DataFrame | None,
    ) -> SiteMatrixProvenanceResult:
        row_drop_stats = {
            "input_rows": int(input_rows),
            "dropped_missing_sequence": dropped_missing_sequence,
            "dropped_incomplete_values": dropped_incomplete_values,
            "missing_data_policy": missing_data_policy,
            "required_observed_count": required_observed_count,
            "deduplicated_site_rows": deduplicated_site_rows,
            "duplicate_site_policy": duplicate_site_policy,
            "retained_rows": int(len(phospho.index)),
        }
        final_phospho = phospho.copy()
        final_site_metadata = site_metadata.copy()
        final_phospho.attrs[self._row_drop_stats_attr] = row_drop_stats.copy()
        final_site_metadata.attrs[self._row_drop_stats_attr] = row_drop_stats.copy()
        final_phospho.attrs[self._site_matrix_policy_attr] = site_matrix_policy
        final_site_metadata.attrs[self._site_matrix_policy_attr] = site_matrix_policy

        site_matrix_provenance = {
            "dropped_missing_sequence_row_ids": dropped_missing_sequence_row_ids,
            "dropped_incomplete_row_ids": dropped_incomplete_row_ids,
            "dropped_row_ids": dropped_row_ids,
            "duplicate_site_policy": duplicate_site_policy,
            "missing_data_policy": missing_data_policy,
            "required_observed_count": required_observed_count,
            "final_constructed_site_ids": tuple(
                str(site_id) for site_id in final_phospho.index.tolist()
            ),
        }
        final_phospho.attrs[self._site_matrix_provenance_attr] = site_matrix_provenance
        final_site_metadata.attrs[self._site_matrix_provenance_attr] = (
            site_matrix_provenance.copy()
        )

        diagnostics = dict(site_matrix_provenance)
        diagnostics["final_constructed_site_ids"] = [
            str(site_id) for site_id in final_phospho.index.tolist()
        ]
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


def _empty_metadata_conflicts() -> pd.DataFrame:
    return dataframe_from_metadata_conflict_rows(())


def _resolve_source_metadata_column(
    site_metadata: pd.DataFrame, column: str
) -> pd.Series:
    if column not in site_metadata.columns:
        return pd.Series(pd.NA, index=site_metadata.index)
    return site_metadata.loc[:, column].copy()


def _normalize_metadata_value(value: object) -> str:
    try:
        if bool(pd.isna(value)):
            return "<NA>"
    except TypeError:
        pass
    if isinstance(value, str):
        trimmed = value.strip()
        return trimmed if trimmed else "<EMPTY>"
    return str(value)


def _is_missing_scalar(value: object) -> bool:
    if value is pd.NA:
        return True
    if isinstance(value, (list, tuple, dict)):
        return False
    try:
        return bool(pd.isna(value))
    except TypeError:
        return False


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
