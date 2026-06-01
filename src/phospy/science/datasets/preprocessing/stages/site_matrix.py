"""Site-matrix construction stage for dataset preprocessing."""

from __future__ import annotations

from dataclasses import replace

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DuplicateSiteResolutionResult,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
    PreprocessingStateTableKey,
    append_row_audit_records,
)
from phospy.science.datasets.preprocessing.policy_models import (
    SiteMatrixDuplicateSitePolicy,
    SiteMatrixMissingDataPolicy,
    SiteMatrixPolicy,
)
from phospy.science.datasets.preprocessing.report_rows import (
    report_rows_from_duplicate_site_resolution_dataframe,
    report_rows_from_metadata_conflicts_dataframe,
    report_rows_from_row_audit_rows,
)
from phospy.science.datasets.preprocessing.report_schema import (
    PreprocessingRowAuditRow,
)
from phospy.science.datasets.preprocessing.stage_contract import (
    PreprocessingStageContract,
)
from phospy.science.datasets.preprocessing.stages.site_matrix_components import (
    DuplicateSiteResolver,
    MetadataConflictDetector,
    MissingDataSiteFilter,
    SequenceSupportFilter,
    SiteMatrixAssembler,
    SiteMatrixProvenanceBuilder,
    SiteMatrixRowAuditBuilder,
)
from phospy.science.sites.identifiers import canonicalize_site_components_series
from phospy.validation.datasets.protein_scoped_site_identity import (
    enforce_site_key_column,
)

_GENE_SYMBOL_COLUMN = "gene_symbol"
_SITE_COLUMN = "site"
_SITE_ID_COLUMN = "site_id"
_SITE_KEY_COLUMN = "site_key"
_REQUIRED_SITE_METADATA_COLUMNS = (
    _GENE_SYMBOL_COLUMN,
    _SITE_COLUMN,
)
_ROW_DROP_STATS_ATTR = "site_matrix_row_drop_stats"
_SITE_MATRIX_POLICY_ATTR = "site_matrix_policy"
_SITE_MATRIX_PROVENANCE_ATTR = "site_matrix_provenance"


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
        sequence_support_filter: SequenceSupportFilter | None = None,
        missing_data_site_filter: MissingDataSiteFilter | None = None,
        site_matrix_assembler: SiteMatrixAssembler | None = None,
        site_matrix_provenance_builder: SiteMatrixProvenanceBuilder | None = None,
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
        self._sequence_support_filter = (
            SequenceSupportFilter()
            if sequence_support_filter is None
            else sequence_support_filter
        )
        self._missing_data_site_filter = (
            MissingDataSiteFilter()
            if missing_data_site_filter is None
            else missing_data_site_filter
        )
        self._site_matrix_assembler = (
            SiteMatrixAssembler()
            if site_matrix_assembler is None
            else site_matrix_assembler
        )
        self._site_matrix_provenance_builder = (
            SiteMatrixProvenanceBuilder(
                row_drop_stats_attr=_ROW_DROP_STATS_ATTR,
                site_matrix_policy_attr=_SITE_MATRIX_POLICY_ATTR,
                site_matrix_provenance_attr=_SITE_MATRIX_PROVENANCE_ATTR,
            )
            if site_matrix_provenance_builder is None
            else site_matrix_provenance_builder
        )

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        policy = state.plan.site_matrix_policy
        if policy is SiteMatrixPolicy.AS_INPUT:
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
        if policy is not SiteMatrixPolicy.BUILD_FROM_METADATA:
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
        constructed_display_id = _build_site_identifier(
            gene_symbol=gene_symbol,
            site=site,
        )
        scientific_row_key = _resolve_scientific_row_key(
            site_metadata=state.site_metadata,
            fallback_constructed_display_id=constructed_display_id,
        )
        sequence_filter_result = self._sequence_support_filter.filter(
            phospho=state.phospho,
            site_metadata=state.site_metadata,
            scientific_row_key=scientific_row_key,
        )

        missing_data_result = self._missing_data_site_filter.filter(
            phospho=sequence_filter_result.phospho,
            scientific_row_key=sequence_filter_result.scientific_row_key,
            missing_data_policy=state.plan.site_matrix_missing_data_policy,
            minimum_observed_values=state.plan.site_matrix_minimum_observed_values,
        )
        policy_filtered_site_metadata = sequence_filter_result.site_metadata.loc[
            missing_data_result.phospho.index
        ]
        policy_filtered_row_key = sequence_filter_result.scientific_row_key.loc[
            missing_data_result.phospho.index
        ]
        policy_filtered_display_id = constructed_display_id.loc[
            missing_data_result.phospho.index
        ]

        duplicate_site_result = self._duplicate_site_resolver.resolve(
            phospho=missing_data_result.phospho,
            site_metadata=policy_filtered_site_metadata,
            scientific_row_key=policy_filtered_row_key,
            constructed_display_id=policy_filtered_display_id,
            duplicate_site_policy=state.plan.site_matrix_duplicate_site_policy,
        )

        assembled = self._site_matrix_assembler.assemble(
            duplicate_site_result=duplicate_site_result,
            output_index_name=(
                None
                if state.phospho.index.name is None
                else str(state.phospho.index.name)
            ),
            dropped_missing_sequence_rows=sequence_filter_result.dropped_rows,
            dropped_incomplete_rows=missing_data_result.dropped_rows,
        )
        row_audit_records = self._row_audit_builder.build(
            dropped_missing_sequence_rows=sequence_filter_result.dropped_rows,
            dropped_incomplete_rows=missing_data_result.dropped_rows,
            duplicate_site_resolution=duplicate_site_result.duplicate_site_resolution,
            site_matrix_policy=state.plan.site_matrix_policy,
            site_matrix_missing_data_policy=state.plan.site_matrix_missing_data_policy,
            site_matrix_duplicate_site_policy=state.plan.site_matrix_duplicate_site_policy,
            required_observed_count=missing_data_result.required_observed_count,
        )
        state_with_row_audit = append_row_audit_records(state, row_audit_records)

        provenance = self._site_matrix_provenance_builder.build(
            phospho=assembled.phospho,
            site_metadata=assembled.site_metadata,
            input_rows=int(len(state.phospho.index)),
            dropped_missing_sequence=sequence_filter_result.dropped_row_count,
            dropped_incomplete_values=missing_data_result.dropped_row_count,
            missing_data_policy=state.plan.site_matrix_missing_data_policy,
            required_observed_count=missing_data_result.required_observed_count,
            deduplicated_site_rows=duplicate_site_result.dropped_row_count,
            duplicate_site_policy=state.plan.site_matrix_duplicate_site_policy,
            site_matrix_policy=policy,
            dropped_missing_sequence_row_ids=assembled.dropped_missing_sequence_row_ids,
            dropped_incomplete_row_ids=assembled.dropped_incomplete_row_ids,
            dropped_row_ids=assembled.dropped_row_ids,
            duplicate_site_resolution=duplicate_site_result.duplicate_site_resolution,
            duplicate_aggregation_diagnostics=(
                duplicate_site_result.duplicate_aggregation_diagnostics
            ),
        )
        if provenance.phospho.empty:
            diagnostics = _format_row_drop_diagnostics(provenance.row_drop_stats)
            raise PhosPyInputError(
                "dataset build request preprocessing site-matrix construction "
                f"produced no retained rows after filtering; {diagnostics}"
            )

        next_state = replace(
            state_with_row_audit,
            phospho=provenance.phospho,
            site_metadata=provenance.site_metadata,
            duplicate_site_resolution=duplicate_site_result.duplicate_site_resolution,
            metadata_conflicts=duplicate_site_result.metadata_conflicts,
        )
        diagnostics = provenance.diagnostics
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
                "dropped_row_ids": assembled.dropped_row_ids,
                "dropped_row_count": int(len(assembled.dropped_row_ids)),
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


def _select_rows_with_usable_sequence_support(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    scientific_row_key: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, int, tuple[tuple[str, str], ...]]:
    result = SequenceSupportFilter().filter(
        phospho=phospho,
        site_metadata=site_metadata,
        scientific_row_key=scientific_row_key,
    )
    return (
        result.phospho,
        result.site_metadata,
        result.scientific_row_key,
        result.dropped_row_count,
        result.dropped_rows,
    )


def _build_site_identifier(
    *,
    gene_symbol: pd.Series,
    site: pd.Series,
) -> pd.Series:
    return canonicalize_site_components_series(
        gene_symbol=gene_symbol,
        site=site,
        field_name=(
            "dataset build request preprocessing site-matrix construction "
            "site_metadata.gene_symbol/site"
        ),
        error_type=PhosPyInputError,
        output_name=_SITE_ID_COLUMN,
    )


def _apply_missing_data_policy(
    *,
    phospho: pd.DataFrame,
    scientific_row_key: pd.Series,
    missing_data_policy: SiteMatrixMissingDataPolicy,
    minimum_observed_values: int | None,
) -> tuple[pd.DataFrame, int, int, tuple[tuple[str, str, int], ...]]:
    result = MissingDataSiteFilter().filter(
        phospho=phospho,
        scientific_row_key=scientific_row_key,
        missing_data_policy=missing_data_policy,
        minimum_observed_values=minimum_observed_values,
    )
    return (
        result.phospho,
        result.dropped_row_count,
        result.required_observed_count,
        result.dropped_rows,
    )


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


def _apply_duplicate_site_policy(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    scientific_row_key: pd.Series,
    constructed_display_id: pd.Series,
    duplicate_site_policy: SiteMatrixDuplicateSitePolicy,
) -> DuplicateSiteResolutionResult:
    """Compatibility wrapper for legacy direct tests of duplicate policy logic."""

    return DuplicateSiteResolver().resolve(
        phospho=phospho,
        site_metadata=site_metadata,
        scientific_row_key=scientific_row_key,
        constructed_display_id=constructed_display_id,
        duplicate_site_policy=duplicate_site_policy,
    )


def _build_metadata_conflicts(
    *,
    site_metadata: pd.DataFrame,
    scientific_row_key: pd.Series,
    constructed_display_id: pd.Series,
) -> pd.DataFrame:
    """Compatibility wrapper for legacy direct metadata-conflict tests."""

    return MetadataConflictDetector().detect(
        site_metadata=site_metadata,
        scientific_row_key=scientific_row_key,
        display_id=constructed_display_id,
    )


def _resolve_scientific_row_key(
    *,
    site_metadata: pd.DataFrame,
    fallback_constructed_display_id: pd.Series,
    prefer_site_key: bool = True,
) -> pd.Series:
    if prefer_site_key and _SITE_KEY_COLUMN in site_metadata.columns:
        try:
            site_keys = enforce_site_key_column(
                site_metadata=site_metadata,
                field_name=(
                    "dataset build request preprocessing site-matrix construction "
                    "site_metadata"
                ),
                error_type=PhosPyInputError,
                column_name=_SITE_KEY_COLUMN,
            )
        except PhosPyInputError:
            site_keys = pd.Series(
                site_metadata.loc[:, _SITE_KEY_COLUMN].astype(str).tolist(),
                index=fallback_constructed_display_id.index.copy(),
                name=_SITE_KEY_COLUMN,
                dtype="object",
            )
        return pd.Series(
            site_keys.astype(str).tolist(),
            index=fallback_constructed_display_id.index.copy(),
            name=_SITE_KEY_COLUMN,
            dtype="object",
        )
    return pd.Series(
        fallback_constructed_display_id.astype(str).tolist(),
        index=fallback_constructed_display_id.index.copy(),
        name=_SITE_ID_COLUMN,
        dtype="object",
    )


def _build_site_matrix_row_audit_records(
    *,
    dropped_missing_sequence_rows: tuple[tuple[str, str], ...],
    dropped_incomplete_rows: tuple[tuple[str, str, int], ...],
    duplicate_site_resolution: pd.DataFrame,
    site_matrix_policy: SiteMatrixPolicy,
    site_matrix_missing_data_policy: SiteMatrixMissingDataPolicy,
    site_matrix_duplicate_site_policy: SiteMatrixDuplicateSitePolicy,
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


def _resolve_operation(plan: PreprocessingPlan) -> str:
    return plan.site_matrix_policy.value


def _resolve_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    return {
        "site_matrix_policy": plan.site_matrix_policy.value,
        "site_matrix_duplicate_site_policy": plan.site_matrix_duplicate_site_policy.value,
        "site_matrix_missing_data_policy": plan.site_matrix_missing_data_policy.value,
        "site_matrix_minimum_observed_values": plan.site_matrix_minimum_observed_values,
    }


SITE_MATRIX_STAGE_CONTRACT = PreprocessingStageContract(
    stage_key=DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    display_label=DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    provenance_stage=DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    operation_name=_resolve_operation,
    serialize_parameters=_resolve_parameters,
    consumed_input_tables=(
        PreprocessingStateTableKey.DATASET_PHOSPHO,
        PreprocessingStateTableKey.DATASET_SITE_METADATA,
    ),
    produced_output_tables=(
        PreprocessingStateTableKey.DATASET_PHOSPHO,
        PreprocessingStateTableKey.DATASET_SITE_METADATA,
        PreprocessingStateTableKey.REPORT_DUPLICATE_SITE_RESOLUTION,
        PreprocessingStateTableKey.REPORT_METADATA_CONFLICTS,
        PreprocessingStateTableKey.REPORT_ROW_AUDIT,
    ),
    stage_factory=SiteMatrixStage,
    backend="pandas",
    diagnostics_metadata={
        "known_diagnostics_fields": (
            "dropped_missing_sequence_row_ids",
            "dropped_incomplete_row_ids",
            "dropped_row_ids",
            "duplicate_site_policy",
            "missing_data_policy",
            "required_observed_count",
            "final_site_keys",
            "final_constructed_site_ids",
            "duplicate_aggregation",
            "duplicate_site_decisions",
        )
    },
)


__all__ = ["SITE_MATRIX_STAGE_CONTRACT", "SiteMatrixStage"]
