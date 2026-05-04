"""Internal dataset preprocessing pipeline orchestration."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace
from typing import TypedDict

import pandas as pd

from phospy.datasets.preprocessing.models import (
    PreprocessingReportRow,
    PreprocessingStage,
    PreprocessingStageExecution,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.datasets.preprocessing.report_rows import (
    validate_preprocessing_report_row,
)
from phospy.datasets.preprocessing.stage_registry import (
    PreprocessingStageMetadata,
    merge_preprocessing_stage_metadata,
)
from phospy.datasets.preprocessing.stages.comparisons import ComparisonsStage
from phospy.datasets.preprocessing.stages.intensity_transform import (
    IntensityTransformStage,
)
from phospy.datasets.preprocessing.stages.missing_data import MissingDataStage
from phospy.datasets.preprocessing.stages.normalisation import NormalisationStage
from phospy.datasets.preprocessing.stages.site_matrix import SiteMatrixStage
from phospy.datasets.preprocessing.stages.site_sequence_resolution import (
    SiteSequenceResolutionStage,
)
from phospy.datasets.preprocessing.stages.total_protein_correction import (
    TotalProteinCorrectionStage,
)
from phospy.errors.build import DatasetBuildError
from phospy.provenance.hashing import fingerprint_optional_table, hash_table
from phospy.provenance.models import (
    PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V2,
    TableFingerprint,
)

_RESERVED_DIAGNOSTIC_KEYS = frozenset(
    {
        "dropped_row_ids",
        "dropped_row_count",
        "imputed_cell_count",
        "imputed_row_ids",
        "notes",
        "diagnostics",
    }
)


class _NormalizedStageDiagnostics(TypedDict):
    dropped_row_ids: tuple[str, ...]
    dropped_row_count: int
    imputed_cell_count: int
    imputed_row_ids: tuple[str, ...]
    notes: str | None
    diagnostics: dict[str, object]


class PreprocessingPipeline:
    """Apply ordered preprocessing stages for interpreted dataset input."""

    def __init__(
        self,
        *,
        stage_registry: tuple[PreprocessingStage, ...] | None = None,
        stage_metadata_registry: tuple[PreprocessingStageMetadata, ...] | None = None,
    ) -> None:
        stages = stage_registry or (
            SiteSequenceResolutionStage(),
            MissingDataStage(),
            IntensityTransformStage(),
            TotalProteinCorrectionStage(),
            SiteMatrixStage(),
            NormalisationStage(),
            ComparisonsStage(),
        )
        self._stages_by_key = {stage.stage_key: stage for stage in stages}
        if len(self._stages_by_key) != len(stages):
            raise DatasetBuildError(
                "dataset preprocessing stage registry contains duplicate stage keys"
            )
        self._stage_metadata_by_key = merge_preprocessing_stage_metadata(
            stage_metadata_registry
        )

    def run(self, state: PreprocessingState) -> PreprocessingState:
        final_state, _ = self.run_with_trace(state)
        return final_state

    def run_with_trace(
        self,
        state: PreprocessingState,
    ) -> tuple[PreprocessingState, tuple[PreprocessingStageExecution, ...]]:
        current = state
        trace: list[PreprocessingStageExecution] = []
        report_rows: list[PreprocessingReportRow] = list(current.report_rows)
        for stage_key in current.plan.stage_order:
            stage = self._stages_by_key.get(stage_key)
            if stage is None:
                raise DatasetBuildError(
                    "dataset preprocessing plan references an unsupported stage: "
                    f"{stage_key}"
                )
            previous = current
            stage_result = stage.run(current)
            if not isinstance(stage_result, PreprocessingStageResult):
                raise DatasetBuildError(
                    "dataset preprocessing stage returned an invalid result payload: "
                    f"{stage_key}"
                )
            current = stage_result.state
            report_rows.extend(_normalize_report_rows(stage_result.report_rows))
            input_hash = hash_table(
                previous.phospho,
                name=f"{stage_key}.input.phospho",
            )
            output_hash = hash_table(
                current.phospho,
                name=f"{stage_key}.output.phospho",
            )
            diagnostics = _normalize_stage_diagnostics(
                raw=stage_result.diagnostics,
                previous=previous,
                current=current,
            )
            stage_metadata = self._stage_metadata_by_key.get(stage_key)
            if stage_metadata is None:
                raise DatasetBuildError(
                    "dataset preprocessing stage metadata is not registered for "
                    f"stage {stage_key!r}"
                )
            consumed_input_tables = _collect_stage_table_fingerprints(
                state=previous,
                table_names=stage_metadata.consumed_input_tables,
            )
            produced_output_tables = _collect_stage_table_fingerprints(
                state=current,
                table_names=stage_metadata.produced_output_tables,
            )
            trace.append(
                PreprocessingStageExecution(
                    stage=stage_metadata.provenance_stage,
                    operation=stage_metadata.operation_name(previous.plan),
                    parameters=stage_metadata.serialize_parameters(previous.plan),
                    input_shape=(
                        int(previous.phospho.shape[0]),
                        int(previous.phospho.shape[1]),
                    ),
                    output_shape=(
                        int(current.phospho.shape[0]),
                        int(current.phospho.shape[1]),
                    ),
                    input_hash=input_hash,
                    output_hash=output_hash,
                    dropped_row_ids=tuple(diagnostics["dropped_row_ids"]),
                    dropped_row_count=int(diagnostics["dropped_row_count"]),
                    schema_version=PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V2,
                    consumed_input_tables=consumed_input_tables,
                    produced_output_tables=produced_output_tables,
                    backend=stage_metadata.backend,
                    random_seed=_resolve_random_seed(
                        diagnostics=diagnostics["diagnostics"]
                    ),
                    is_deterministic=True,
                    imputed_cell_count=int(diagnostics["imputed_cell_count"]),
                    imputed_row_ids=tuple(diagnostics["imputed_row_ids"]),
                    notes=diagnostics["notes"],
                    diagnostics=dict(diagnostics["diagnostics"]),
                )
            )
        if report_rows:
            current = replace(current, report_rows=tuple(report_rows))
        return current, tuple(trace)


def _normalize_stage_diagnostics(
    *,
    raw: Mapping[str, object],
    previous: PreprocessingState,
    current: PreprocessingState,
) -> _NormalizedStageDiagnostics:
    default_dropped_row_ids = _resolve_dropped_row_ids(
        before=previous.phospho.index,
        after=current.phospho.index,
    )
    dropped_row_ids = _coerce_string_tuple(
        raw.get("dropped_row_ids", default_dropped_row_ids)
    )
    dropped_row_count = _coerce_int(
        raw.get("dropped_row_count"),
        default=len(dropped_row_ids),
    )

    default_imputed_row_ids, default_imputed_cell_count = _resolve_imputation_summary(
        before=previous.phospho,
        after=current.phospho,
    )
    imputed_row_ids = _coerce_string_tuple(
        raw.get("imputed_row_ids", default_imputed_row_ids)
    )
    imputed_cell_count = _coerce_int(
        raw.get("imputed_cell_count"),
        default=default_imputed_cell_count,
    )

    notes_raw = raw.get("notes", "stage executed")
    notes = None if notes_raw is None else str(notes_raw)

    nested_diagnostics = raw.get("diagnostics")
    if isinstance(nested_diagnostics, Mapping):
        diagnostics = {str(key): value for key, value in nested_diagnostics.items()}
    else:
        diagnostics = {
            str(key): value
            for key, value in raw.items()
            if key not in _RESERVED_DIAGNOSTIC_KEYS
        }

    return {
        "dropped_row_ids": dropped_row_ids,
        "dropped_row_count": dropped_row_count,
        "imputed_cell_count": imputed_cell_count,
        "imputed_row_ids": imputed_row_ids,
        "notes": notes,
        "diagnostics": diagnostics,
    }


def _coerce_string_tuple(value: object) -> tuple[str, ...]:
    if isinstance(value, tuple):
        return tuple(str(item) for item in value)
    if isinstance(value, list):
        return tuple(str(item) for item in value)
    return ()


def _coerce_int(value: object, *, default: int) -> int:
    if value is None:
        return int(default)
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value)
    if isinstance(value, str):
        return int(value.strip())
    return int(default)


def _normalize_report_rows(
    rows: Sequence[PreprocessingReportRow],
) -> tuple[PreprocessingReportRow, ...]:
    normalized: list[PreprocessingReportRow] = []
    for row in rows:
        normalized.append(validate_preprocessing_report_row(row))
    return tuple(normalized)


def _resolve_dropped_row_ids(*, before: pd.Index, after: pd.Index) -> tuple[str, ...]:
    after_values = {str(value) for value in after.tolist()}
    return tuple(
        str(value) for value in before.tolist() if str(value) not in after_values
    )


def _resolve_imputation_summary(
    *,
    before: pd.DataFrame,
    after: pd.DataFrame,
) -> tuple[tuple[str, ...], int]:
    if before.empty or after.empty:
        return (), 0
    aligned_before = before.reindex(after.index)
    imputed_mask = aligned_before.isna() & after.notna()
    imputed_cell_count = int(imputed_mask.to_numpy().sum())
    if imputed_cell_count == 0:
        return (), 0
    row_flags = imputed_mask.any(axis=1).to_numpy(dtype=bool, copy=False)
    row_ids = tuple(
        str(site_id)
        for site_id, flagged in zip(after.index.tolist(), row_flags, strict=True)
        if bool(flagged)
    )
    return row_ids, imputed_cell_count


def _resolve_random_seed(*, diagnostics: Mapping[str, object]) -> int | None:
    value = diagnostics.get("random_seed")
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    return None


def _collect_stage_table_fingerprints(
    *,
    state: PreprocessingState,
    table_names: tuple[str, ...],
) -> tuple[TableFingerprint, ...]:
    fingerprints: list[TableFingerprint] = []
    for table_name in table_names:
        table = _resolve_state_table(state=state, table_name=table_name)
        fingerprint = fingerprint_optional_table(table, name=table_name)
        if fingerprint is None:
            continue
        fingerprints.append(fingerprint)
    return tuple(fingerprints)


def _resolve_state_table(
    *,
    state: PreprocessingState,
    table_name: str,
) -> pd.DataFrame | None:
    if table_name == "dataset.phospho":
        return state.phospho
    if table_name == "dataset.site_metadata":
        return state.site_metadata
    if table_name == "dataset.sample_metadata":
        return state.sample_metadata
    if table_name == "dataset.total":
        return state.total
    if table_name == "dataset.comparisons":
        return state.comparisons
    if table_name == "report.comparison_group_stats":
        return state.comparison_group_stats
    if table_name == "report.comparison_pair_stats":
        return state.comparison_pair_stats
    if table_name == "report.duplicate_site_resolution":
        return state.duplicate_site_resolution
    if table_name == "report.metadata_conflicts":
        return state.metadata_conflicts
    if table_name == "report.row_audit":
        return state.row_audit
    return None


__all__ = ["PreprocessingPipeline"]
