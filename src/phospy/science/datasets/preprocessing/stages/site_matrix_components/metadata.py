from __future__ import annotations

from dataclasses import dataclass
from typing import cast

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.policy_models import (
    SiteMatrixDuplicateSitePolicy,
    SiteMatrixMissingDataPolicy,
    SiteMatrixPolicy,
)
from phospy.science.datasets.preprocessing.report_schema import (
    METADATA_CONFLICT_COLUMNS,
    dataframe_from_metadata_conflict_rows,
)
from phospy.science.sites.identifiers import canonicalize_site_series
from phospy.science.sites.validation import require_site_key_series

_DEFAULT_METADATA_CONFLICT_FIELDS = (
    "protein_id",
    "gene_symbol",
    "site",
    "site_sequence",
    "localisation_probability",
    "residue",
    "site_position",
    "position",
)
_DEFAULT_ROW_DROP_STATS_ATTR = "site_matrix_row_drop_stats"
_DEFAULT_SITE_MATRIX_POLICY_ATTR = "site_matrix_policy"
_DEFAULT_SITE_MATRIX_PROVENANCE_ATTR = "site_matrix_provenance"
_SITE_KEY_COLUMN = "site_key"


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
        }
        site_matrix_provenance["final_constructed_site_ids"] = tuple(
            str(site_id) for site_id in final_phospho.index.tolist()
        )
        if duplicate_aggregation_diagnostics is not None:
            site_matrix_provenance["duplicate_aggregation"] = dict(
                duplicate_aggregation_diagnostics
            )
        final_phospho.attrs[self._site_matrix_provenance_attr] = site_matrix_provenance
        final_site_metadata.attrs[self._site_matrix_provenance_attr] = (
            site_matrix_provenance.copy()
        )

        diagnostics = dict(site_matrix_provenance)
        diagnostics["final_site_keys"] = list(
            str(site_id) for site_id in final_phospho.index.tolist()
        )
        diagnostics["final_constructed_site_ids"] = list(
            diagnostics["final_constructed_site_ids"]
        )
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
    """Detect metadata conflicts for rows mapping to the same scientific row key."""

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
        scientific_row_key: pd.Series | None = None,
        display_id: pd.Series | None = None,
    ) -> pd.DataFrame:
        if site_metadata.empty:
            return _empty_metadata_conflicts()
        if scientific_row_key is None:
            raise PhosPyInputError(
                "metadata conflict detection requires site_key row identity"
            )
        if display_id is None:
            raise PhosPyInputError(
                "metadata conflict detection requires display_id metadata"
            )
        if scientific_row_key.name != _SITE_KEY_COLUMN:
            raise PhosPyInputError(
                "metadata conflict detection requires scientific_row_key.name='site_key'"
            )
        if not scientific_row_key.index.equals(site_metadata.index):
            raise PhosPyInputError(
                "metadata conflict detection requires site_key values aligned to "
                "site_metadata rows"
            )
        if not display_id.index.equals(site_metadata.index):
            raise PhosPyInputError(
                "metadata conflict detection requires display_id values aligned to "
                "site_metadata rows"
            )
        site_keys = require_site_key_series(
            scientific_row_key.astype("object"),
            field_name="metadata conflict detection site_key",
            error_type=PhosPyInputError,
        )
        display_ids = canonicalize_site_series(
            display_id.astype("object"),
            field_name="metadata conflict detection display_id",
            error_type=PhosPyInputError,
        )

        duplicate_groups = site_metadata.assign(
            **{
                _SITE_KEY_COLUMN: site_keys.astype(str),
                "display_id": display_ids.astype(str),
                "source_row_id": site_metadata.index.astype(str),
            }
        )
        conflict_fields = (
            tuple(site_metadata.columns)
            if self._conflict_fields is None
            else self._conflict_fields
        )
        records: list[dict[str, object]] = []
        for site_key, group in duplicate_groups.groupby(_SITE_KEY_COLUMN, sort=False):
            source_row_ids = tuple(group.loc[:, "source_row_id"].astype(str).tolist())
            display_label = (
                group.loc[:, "display_id"].astype(str).iloc[0]
                if "display_id" in group.columns
                else ""
            )
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
                        "site_key": str(site_key),
                        "display_id": str(display_label),
                        "site_id": str(display_label),
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
    scientific_row_key: pd.Series,
    display_id: pd.Series,
    metadata_conflicts: pd.DataFrame,
) -> pd.DataFrame:
    if site_metadata.empty:
        return site_metadata.copy()

    metadata_columns = list(site_metadata.columns)
    grouped_metadata = cast(
        pd.DataFrame,
        site_metadata.assign(
            **{
                _SITE_KEY_COLUMN: scientific_row_key.to_numpy(),
                "display_id": display_id.to_numpy(),
            }
        )
        .groupby(_SITE_KEY_COLUMN, sort=False)[metadata_columns]
        .first(),
    )
    grouped_metadata.index = pd.Index(
        grouped_metadata.index.astype(str), name=_SITE_KEY_COLUMN
    )

    if metadata_conflicts.empty:
        return cast(pd.DataFrame, grouped_metadata)

    conflict_records = metadata_conflicts.loc[
        :, ["site_key", "field"]
    ].drop_duplicates()
    for conflict in conflict_records.to_dict(orient="records"):
        site_key = str(conflict["site_key"])
        field = str(conflict["field"])
        if site_key in grouped_metadata.index and field in grouped_metadata.columns:
            grouped_metadata.at[site_key, field] = pd.NA

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
