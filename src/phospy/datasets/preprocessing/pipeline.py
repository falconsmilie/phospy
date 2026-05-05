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
    build_registered_preprocessing_stage_instances,
    resolve_registered_preprocessing_stages,
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
        resolved_metadata_registry = resolve_registered_preprocessing_stages(
            stage_metadata_registry
        )
        stages = stage_registry or build_registered_preprocessing_stage_instances(
            resolved_metadata_registry
        )
        self._stages_by_key = {stage.stage_key: stage for stage in stages}
        if len(self._stages_by_key) != len(stages):
            raise DatasetBuildError(
                "dataset preprocessing stage registry contains duplicate stage keys"
            )
        self._stage_metadata_by_key = {
            metadata.stage_key: metadata for metadata in resolved_metadata_registry
        }

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
                stage_key=stage_key,
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
                    stage=stage_metadata.provenance_stage_key,
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
                        stage_key=stage_key, diagnostics=diagnostics["diagnostics"]
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
    stage_key: str,
    raw: Mapping[str, object],
    previous: PreprocessingState,
    current: PreprocessingState,
) -> _NormalizedStageDiagnostics:
    default_dropped_row_ids = _resolve_dropped_row_ids(
        before=previous.phospho.index,
        after=current.phospho.index,
    )
    dropped_row_ids = _coerce_string_tuple(
        raw.get("dropped_row_ids", default_dropped_row_ids),
        stage_key=stage_key,
        field_name="dropped_row_ids",
    )
    dropped_row_count = _coerce_int(
        raw.get("dropped_row_count"),
        stage_key=stage_key,
        field_name="dropped_row_count",
        default=len(dropped_row_ids),
    )

    default_imputed_row_ids, default_imputed_cell_count = _resolve_imputation_summary(
        before=previous.phospho,
        after=current.phospho,
    )
    imputed_row_ids = _coerce_string_tuple(
        raw.get("imputed_row_ids", default_imputed_row_ids),
        stage_key=stage_key,
        field_name="imputed_row_ids",
    )
    imputed_cell_count = _coerce_int(
        raw.get("imputed_cell_count"),
        stage_key=stage_key,
        field_name="imputed_cell_count",
        default=default_imputed_cell_count,
    )

    notes_raw = raw.get("notes", "stage executed")
    if notes_raw is None:
        notes = None
    elif isinstance(notes_raw, str):
        notes = notes_raw
    else:
        raise DatasetBuildError(
            "dataset preprocessing stage diagnostics parse error: "
            f"stage={stage_key!r}, field='notes', expected string or null, got "
            f"{notes_raw!r} ({type(notes_raw).__name__})"
        )

    nested_diagnostics = raw.get("diagnostics")
    if isinstance(nested_diagnostics, Mapping):
        diagnostics: dict[str, object] = {}
        for key, value in nested_diagnostics.items():
            if not isinstance(key, str):
                raise DatasetBuildError(
                    "dataset preprocessing stage diagnostics parse error: "
                    f"stage={stage_key!r}, field='diagnostics.<key>', expected "
                    f"string, got {key!r} ({type(key).__name__})"
                )
            diagnostics[key] = value
    elif "diagnostics" in raw:
        raise DatasetBuildError(
            "dataset preprocessing stage diagnostics parse error: "
            f"stage={stage_key!r}, field='diagnostics', expected object, got "
            f"{nested_diagnostics!r} ({type(nested_diagnostics).__name__})"
        )
    else:
        diagnostics = {}
        for key, value in raw.items():
            if key in _RESERVED_DIAGNOSTIC_KEYS:
                continue
            if not isinstance(key, str):
                raise DatasetBuildError(
                    "dataset preprocessing stage diagnostics parse error: "
                    f"stage={stage_key!r}, field='diagnostics.<key>', expected "
                    f"string, got {key!r} ({type(key).__name__})"
                )
            diagnostics[key] = value

    return {
        "dropped_row_ids": dropped_row_ids,
        "dropped_row_count": dropped_row_count,
        "imputed_cell_count": imputed_cell_count,
        "imputed_row_ids": imputed_row_ids,
        "notes": notes,
        "diagnostics": diagnostics,
    }


def _coerce_string_tuple(
    value: object,
    *,
    stage_key: str,
    field_name: str,
) -> tuple[str, ...]:
    if isinstance(value, tuple):
        return _require_string_sequence(
            value,
            stage_key=stage_key,
            field_name=field_name,
        )
    if isinstance(value, list):
        return _require_string_sequence(
            value,
            stage_key=stage_key,
            field_name=field_name,
        )
    raise DatasetBuildError(
        "dataset preprocessing stage diagnostics parse error: "
        f"stage={stage_key!r}, field={field_name!r}, expected array of strings, "
        f"got {value!r} ({type(value).__name__})"
    )


def _require_string_sequence(
    value: tuple[object, ...] | list[object],
    *,
    stage_key: str,
    field_name: str,
) -> tuple[str, ...]:
    resolved: list[str] = []
    for index, item in enumerate(value):
        if not isinstance(item, str):
            raise DatasetBuildError(
                "dataset preprocessing stage diagnostics parse error: "
                f"stage={stage_key!r}, field={field_name!r}[{index}], expected "
                f"string, got {item!r} ({type(item).__name__})"
            )
        resolved.append(item)
    return tuple(resolved)


def _coerce_int(
    value: object,
    *,
    stage_key: str,
    field_name: str,
    default: int,
) -> int:
    if value is None:
        return default
    if isinstance(value, bool):
        raise DatasetBuildError(
            "dataset preprocessing stage diagnostics parse error: "
            f"stage={stage_key!r}, field={field_name!r}, expected int (bool is not "
            f"accepted), got {value!r} (bool)"
        )
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        raise DatasetBuildError(
            "dataset preprocessing stage diagnostics parse error: "
            f"stage={stage_key!r}, field={field_name!r}, expected int (floats are "
            f"not accepted), got {value!r} (float)"
        )
    if isinstance(value, str):
        raise DatasetBuildError(
            "dataset preprocessing stage diagnostics parse error: "
            f"stage={stage_key!r}, field={field_name!r}, expected int, got "
            f"{value!r} (str)"
        )
    raise DatasetBuildError(
        "dataset preprocessing stage diagnostics parse error: "
        f"stage={stage_key!r}, field={field_name!r}, expected int, got "
        f"{value!r} ({type(value).__name__})"
    )


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


def _resolve_random_seed(
    *,
    stage_key: str,
    diagnostics: Mapping[str, object],
) -> int | None:
    value = diagnostics.get("random_seed")
    if value is None:
        return None
    return _coerce_int(
        value,
        stage_key=stage_key,
        field_name="diagnostics.random_seed",
        default=0,
    )


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
