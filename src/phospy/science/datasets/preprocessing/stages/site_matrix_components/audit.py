from __future__ import annotations

import pandas as pd

from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
)
from phospy.science.datasets.preprocessing.policy_models import (
    SiteMatrixDuplicateSitePolicy,
    SiteMatrixMissingDataPolicy,
    SiteMatrixPolicy,
)
from phospy.science.datasets.preprocessing.report_schema import PreprocessingRowAuditRow
from phospy.science.datasets.preprocessing.stages.site_matrix_components.metadata import (
    _is_missing_scalar,
)


class SiteMatrixRowAuditBuilder:
    """Build row-audit records for site-matrix sequence/missing/duplicate decisions."""

    def build(
        self,
        *,
        dropped_missing_sequence_rows: tuple[tuple[str, str], ...],
        dropped_incomplete_rows: tuple[tuple[str, str, int], ...],
        duplicate_site_resolution: pd.DataFrame,
        site_matrix_policy: SiteMatrixPolicy | str,
        site_matrix_missing_data_policy: SiteMatrixMissingDataPolicy | str,
        site_matrix_duplicate_site_policy: SiteMatrixDuplicateSitePolicy | str,
        required_observed_count: int,
    ) -> list[PreprocessingRowAuditRow]:
        resolved_site_matrix_policy = SiteMatrixPolicy.parse(
            site_matrix_policy,
            field_name="site_matrix.policy",
        )
        resolved_site_matrix_missing_data_policy = SiteMatrixMissingDataPolicy.parse(
            site_matrix_missing_data_policy,
            field_name="site_matrix.missing_data_policy",
        )
        resolved_site_matrix_duplicate_site_policy = (
            SiteMatrixDuplicateSitePolicy.parse(
                site_matrix_duplicate_site_policy,
                field_name="site_matrix.duplicate_site_policy",
            )
        )
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
                        "site_matrix_policy": resolved_site_matrix_policy.value,
                        "site_matrix_missing_data_policy": (
                            resolved_site_matrix_missing_data_policy.value
                        ),
                        "site_matrix_duplicate_site_policy": (
                            resolved_site_matrix_duplicate_site_policy.value
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
                        "site_matrix_policy": resolved_site_matrix_policy.value,
                        "site_matrix_missing_data_policy": (
                            resolved_site_matrix_missing_data_policy.value
                        ),
                        "site_matrix_duplicate_site_policy": (
                            resolved_site_matrix_duplicate_site_policy.value
                        ),
                        "observed_values": int(observed_value_count),
                        "required_observed_count": int(required_observed_count),
                    },
                )
            )
        records.extend(
            self._build_duplicate_site_row_audit_records(
                duplicate_site_resolution=duplicate_site_resolution,
                site_matrix_policy=resolved_site_matrix_policy,
                site_matrix_duplicate_site_policy=resolved_site_matrix_duplicate_site_policy,
            )
        )
        return records

    @staticmethod
    def _build_duplicate_site_row_audit_records(
        *,
        duplicate_site_resolution: pd.DataFrame,
        site_matrix_policy: SiteMatrixPolicy,
        site_matrix_duplicate_site_policy: SiteMatrixDuplicateSitePolicy,
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
            SiteMatrixDuplicateSitePolicy.AGGREGATE_MEAN,
            SiteMatrixDuplicateSitePolicy.AGGREGATE_MEDIAN,
        }
        records: list[PreprocessingRowAuditRow] = []
        for row in duplicate_site_resolution.to_dict(orient="records"):
            site_id = str(row.get("site_id", ""))
            source_row_id = str(row.get("source_row_id", ""))
            source_rows = source_rows_by_site.get(site_id, (source_row_id,))
            if aggregated:
                action = "aggregated"
                retained = True
                retained_row_id: object = site_id
                retained_row: object = site_id
                reason = _optional_text(
                    row.get("retained_reason"),
                    fallback=(
                        "contributed to site-level aggregate from duplicate source rows"
                    ),
                )
            else:
                retained = bool(row.get("retained"))
                action = "retained" if retained else "collapsed"
                retained_row_id = retained_row_by_site.get(site_id, pd.NA)
                retained_row = (
                    retained_row_id if retained_row_id is not pd.NA else pd.NA
                )
                reason = _optional_text(
                    row.get("retained_reason")
                    if retained
                    else row.get("dropped_reason"),
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
                        "site_matrix_policy": site_matrix_policy.value,
                        "duplicate_site_policy": site_matrix_duplicate_site_policy.value,
                        "site_matrix_duplicate_site_policy": (
                            site_matrix_duplicate_site_policy.value
                        ),
                        "observed_values": _optional_int(row.get("observed_values")),
                        "mean_signal": _optional_float(row.get("mean_signal")),
                        "metadata_conflict_detected": bool(
                            _optional_bool(row.get("metadata_conflict_detected"))
                        ),
                    },
                )
            )
        return records


def _optional_text(value: object, *, fallback: str) -> str:
    if _is_missing_scalar(value):
        return fallback
    text = str(value).strip()
    return text if text else fallback


def _optional_int(value: object | None) -> int | None:
    if _is_missing_scalar(value):
        return None
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if value.is_integer():
            return int(value)
        return int(value)
    if isinstance(value, str):
        return int(value.strip())
    raise TypeError(f"cannot coerce value of type {type(value).__name__} to int")


def _optional_float(value: object | None) -> float | None:
    if _is_missing_scalar(value):
        return None
    if isinstance(value, (int, float, bool)):
        return float(value)
    if isinstance(value, str):
        return float(value.strip())
    raise TypeError(f"cannot coerce value of type {type(value).__name__} to float")


def _optional_bool(value: object | None) -> bool | None:
    if _is_missing_scalar(value):
        return None
    return bool(value)
