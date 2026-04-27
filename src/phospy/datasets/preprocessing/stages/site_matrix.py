"""Site-matrix construction stage for dataset preprocessing."""

from __future__ import annotations

import re
from dataclasses import replace

import pandas as pd

from phospy.api.configs import (
    DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
    DATASET_SITE_MATRIX_POLICY_AS_INPUT,
    DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA,
)
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DuplicateSiteResolutionResult,
    PreprocessingStageResult,
    PreprocessingState,
    append_row_audit_records,
)
from phospy.datasets.preprocessing.report_rows import (
    report_rows_from_duplicate_site_resolution_dataframe,
    report_rows_from_metadata_conflicts_dataframe,
    report_rows_from_row_audit_rows,
)
from phospy.datasets.preprocessing.report_schema import (
    PreprocessingRowAuditRow,
)
from phospy.datasets.preprocessing.stages.site_matrix_components import (
    DuplicateSiteResolver,
    MetadataConflictDetector,
    SiteMatrixRowAuditBuilder,
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
_SITE_MATRIX_PROVENANCE_ATTR = "site_matrix_provenance"
_INTERNAL_SITE_MATRIX_MISSING_DATA_POLICY_RETAIN_MISSING = "retain_missing"
_INTERNAL_SITE_MATRIX_MISSING_DATA_POLICY_REQUIRE_MIN_OBSERVED_VALUES = (
    "require_min_observed_values"
)


class SiteMatrixStage:
    """Build site-matrix-ready phospho rows from site metadata when requested.

    This stage ports the historical-baseline site-matrix policy surface behind
    `site_matrix.policy='build_from_metadata'`.
    """

    stage_key = DATASET_PREPROCESSING_STAGE_SITE_MATRIX

    def __init__(
        self,
        *,
        duplicate_site_resolver: DuplicateSiteResolver | None = None,
        row_audit_builder: SiteMatrixRowAuditBuilder | None = None,
    ) -> None:
        self._duplicate_site_resolver = (
            DuplicateSiteResolver()
            if duplicate_site_resolver is None
            else duplicate_site_resolver
        )
        self._row_audit_builder = (
            SiteMatrixRowAuditBuilder()
            if row_audit_builder is None
            else row_audit_builder
        )

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        policy = state.plan.site_matrix_policy
        if policy == DATASET_SITE_MATRIX_POLICY_AS_INPUT:
            return PreprocessingStageResult(
                state=state,
                diagnostics={
                    "dropped_row_ids": (),
                    "dropped_row_count": 0,
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "stage executed",
                    "diagnostics": {},
                },
            )
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
        constructed_site_id = _build_site_identifier(gene_symbol=gene_symbol, site=site)
        (
            with_sequence_phospho,
            with_sequence_site_metadata,
            with_sequence_site_id,
            dropped_missing_sequence,
            dropped_missing_sequence_rows,
        ) = _select_rows_with_usable_sequence_support(
            phospho=state.phospho,
            site_metadata=state.site_metadata,
            constructed_site_id=constructed_site_id,
        )

        (
            policy_filtered_phospho,
            dropped_incomplete_values,
            required_observed_count,
            dropped_incomplete_rows,
        ) = _apply_missing_data_policy(
            phospho=with_sequence_phospho,
            constructed_site_id=with_sequence_site_id,
            missing_data_policy=state.plan.site_matrix_missing_data_policy,
            minimum_observed_values=state.plan.site_matrix_minimum_observed_values,
        )
        policy_filtered_site_metadata = with_sequence_site_metadata.loc[
            policy_filtered_phospho.index
        ]
        policy_filtered_site_id = with_sequence_site_id.loc[
            policy_filtered_phospho.index
        ]

        duplicate_site_result = self._duplicate_site_resolver.resolve(
            phospho=policy_filtered_phospho,
            site_metadata=policy_filtered_site_metadata,
            constructed_site_id=policy_filtered_site_id,
            duplicate_site_policy=state.plan.site_matrix_duplicate_site_policy,
        )

        final_phospho = duplicate_site_result.phospho.sort_index(kind="stable")
        final_site_metadata = duplicate_site_result.site_metadata.reindex(
            final_phospho.index
        )
        final_site_index = pd.Index(
            final_phospho.index.tolist(), name=state.phospho.index.name
        )
        final_phospho.index = final_site_index
        final_site_metadata.index = final_site_index.copy()
        duplicate_dropped_row_ids = tuple(
            str(row_id)
            for row_id in duplicate_site_result.duplicate_site_resolution.loc[
                ~duplicate_site_result.duplicate_site_resolution.loc[:, "retained"],
                "source_row_id",
            ]
            .astype(str)
            .tolist()
        )
        dropped_missing_sequence_row_ids = tuple(
            row_id for row_id, _ in dropped_missing_sequence_rows
        )
        dropped_incomplete_row_ids = tuple(
            row_id for row_id, _, _ in dropped_incomplete_rows
        )
        dropped_row_ids = _unique_strings_preserve_order(
            (
                *dropped_missing_sequence_row_ids,
                *dropped_incomplete_row_ids,
                *duplicate_dropped_row_ids,
            )
        )
        row_audit_records = self._row_audit_builder.build(
            dropped_missing_sequence_rows=dropped_missing_sequence_rows,
            dropped_incomplete_rows=dropped_incomplete_rows,
            duplicate_site_resolution=duplicate_site_result.duplicate_site_resolution,
            site_matrix_policy=state.plan.site_matrix_policy,
            site_matrix_missing_data_policy=state.plan.site_matrix_missing_data_policy,
            site_matrix_duplicate_site_policy=state.plan.site_matrix_duplicate_site_policy,
            required_observed_count=required_observed_count,
        )
        state_with_row_audit = append_row_audit_records(state, row_audit_records)

        row_drop_stats = {
            "input_rows": int(len(state.phospho.index)),
            "dropped_missing_sequence": dropped_missing_sequence,
            "dropped_incomplete_values": dropped_incomplete_values,
            "missing_data_policy": state.plan.site_matrix_missing_data_policy,
            "required_observed_count": required_observed_count,
            "deduplicated_site_rows": duplicate_site_result.dropped_row_count,
            "duplicate_site_policy": state.plan.site_matrix_duplicate_site_policy,
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
        site_matrix_provenance = {
            "dropped_missing_sequence_row_ids": dropped_missing_sequence_row_ids,
            "dropped_incomplete_row_ids": dropped_incomplete_row_ids,
            "dropped_row_ids": dropped_row_ids,
            "duplicate_site_policy": state.plan.site_matrix_duplicate_site_policy,
            "missing_data_policy": state.plan.site_matrix_missing_data_policy,
            "required_observed_count": required_observed_count,
            "final_constructed_site_ids": tuple(
                str(site_id) for site_id in final_phospho.index.tolist()
            ),
        }
        final_phospho.attrs[_SITE_MATRIX_PROVENANCE_ATTR] = site_matrix_provenance
        final_site_metadata.attrs[_SITE_MATRIX_PROVENANCE_ATTR] = (
            site_matrix_provenance.copy()
        )

        next_state = replace(
            state_with_row_audit,
            phospho=final_phospho,
            site_metadata=final_site_metadata,
            duplicate_site_resolution=duplicate_site_result.duplicate_site_resolution,
            metadata_conflicts=duplicate_site_result.metadata_conflicts,
        )
        diagnostics = dict(site_matrix_provenance)
        diagnostics["final_constructed_site_ids"] = [
            str(site_id) for site_id in final_phospho.index.tolist()
        ]
        if duplicate_site_result.duplicate_site_resolution is not None:
            diagnostics["duplicate_site_decisions"] = _records_from_frame(
                duplicate_site_result.duplicate_site_resolution
            )
        stage_report_rows = (
            report_rows_from_row_audit_rows(row_audit_records)
            + report_rows_from_duplicate_site_resolution_dataframe(
                duplicate_site_result.duplicate_site_resolution
            )
            + report_rows_from_metadata_conflicts_dataframe(
                duplicate_site_result.metadata_conflicts
            )
        )
        return PreprocessingStageResult(
            state=next_state,
            report_rows=stage_report_rows,
            diagnostics={
                "dropped_row_ids": dropped_row_ids,
                "dropped_row_count": int(len(dropped_row_ids)),
                "imputed_cell_count": 0,
                "imputed_row_ids": (),
                "notes": "stage executed",
                "diagnostics": diagnostics,
            },
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


def _select_rows_with_usable_sequence_support(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    constructed_site_id: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, int, tuple[tuple[str, str], ...]]:
    if _SITE_SEQUENCE_COLUMN in site_metadata.columns:
        site_sequence = _resolve_optional_string_column(
            site_metadata,
            column_name=_SITE_SEQUENCE_COLUMN,
        )
    else:
        site_sequence = pd.Series(
            pd.NA,
            index=site_metadata.index.copy(),
            dtype="string",
            name=_SITE_SEQUENCE_COLUMN,
        )
    has_sequence = site_sequence.notna()
    dropped_rows = tuple(
        (str(row_id), str(site_id))
        for row_id, site_id in zip(
            phospho.index[~has_sequence].tolist(),
            constructed_site_id.loc[~has_sequence].astype(str).tolist(),
            strict=True,
        )
    )
    return (
        phospho.loc[has_sequence],
        site_metadata.loc[has_sequence],
        constructed_site_id.loc[has_sequence],
        int((~has_sequence).sum()),
        dropped_rows,
    )


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
        preview = ", ".join(normalized_site.loc[invalid_site].astype(str).head(3))
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
    constructed_site_id: pd.Series,
    missing_data_policy: str,
    minimum_observed_values: int | None,
) -> tuple[pd.DataFrame, int, int, tuple[tuple[str, str, int], ...]]:
    observed_counts = phospho.notna().sum(axis=1)
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
    dropped_row_details = tuple(
        (str(row_id), str(site_id), int(observed_value_count))
        for row_id, site_id, observed_value_count in zip(
            phospho.index[~retained_mask].tolist(),
            constructed_site_id.loc[~retained_mask].astype(str).tolist(),
            observed_counts.loc[~retained_mask].tolist(),
            strict=True,
        )
    )
    return filtered, dropped_rows, required_observed_count, dropped_row_details


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


def _is_missing_scalar(value: object) -> bool:
    if value is pd.NA:
        return True
    if isinstance(value, (list, tuple, dict)):
        return False
    try:
        return bool(pd.isna(value))
    except TypeError:
        return False


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
        "duplicate_site_policy="
        f"{str(row_drop_stats.get('duplicate_site_policy', 'max_mean_signal'))}, "
        f"other_dropped_rows={other_dropped_rows}, "
        f"retained_rows={retained_rows}"
    )


def _unique_strings_preserve_order(values: tuple[str, ...]) -> tuple[str, ...]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return tuple(ordered)


def _apply_duplicate_site_policy(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    constructed_site_id: pd.Series,
    duplicate_site_policy: str,
) -> DuplicateSiteResolutionResult:
    """Compatibility wrapper for legacy direct tests of duplicate policy logic."""

    return DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        constructed_site_id=constructed_site_id,
        duplicate_site_policy=duplicate_site_policy,
    )


def _build_metadata_conflicts(
    *,
    site_metadata: pd.DataFrame,
    constructed_site_id: pd.Series,
) -> pd.DataFrame:
    """Compatibility wrapper for legacy direct metadata-conflict tests."""

    return MetadataConflictDetector().detect(
        site_metadata=site_metadata,
        constructed_site_id=constructed_site_id,
    )


def _build_site_matrix_row_audit_records(
    *,
    dropped_missing_sequence_rows: tuple[tuple[str, str], ...],
    dropped_incomplete_rows: tuple[tuple[str, str, int], ...],
    duplicate_site_resolution: pd.DataFrame,
    site_matrix_policy: str,
    site_matrix_missing_data_policy: str,
    site_matrix_duplicate_site_policy: str,
    required_observed_count: int,
) -> list[PreprocessingRowAuditRow]:
    """Compatibility wrapper for legacy row-audit construction tests."""

    return SiteMatrixRowAuditBuilder().build(
        dropped_missing_sequence_rows=dropped_missing_sequence_rows,
        dropped_incomplete_rows=dropped_incomplete_rows,
        duplicate_site_resolution=duplicate_site_resolution,
        site_matrix_policy=site_matrix_policy,
        site_matrix_missing_data_policy=site_matrix_missing_data_policy,
        site_matrix_duplicate_site_policy=site_matrix_duplicate_site_policy,
        required_observed_count=required_observed_count,
    )


__all__ = ["SiteMatrixStage"]
