"""Private preprocessing stage diagnostics normalization services."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

import pandas as pd

from phospy.errors.build import DatasetBuildError
from phospy.science.datasets.preprocessing.models import PreprocessingState

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


@dataclass(frozen=True, slots=True)
class _StageDiagnosticsDefaults:
    dropped_row_ids: tuple[str, ...]
    imputed_row_ids: tuple[str, ...]
    imputed_cell_count: int


@dataclass(frozen=True, slots=True)
class _NormalizedStageDiagnostics:
    dropped_row_ids: tuple[str, ...]
    dropped_row_count: int
    imputed_cell_count: int
    imputed_row_ids: tuple[str, ...]
    notes: str | None
    diagnostics: Mapping[str, object]


class _StageDiagnosticsDefaultsResolver:
    """Resolve compatibility defaults that depend on stage state transitions."""

    def run(
        self,
        *,
        previous: PreprocessingState,
        current: PreprocessingState,
    ) -> _StageDiagnosticsDefaults:
        imputed_row_ids, imputed_cell_count = _resolve_imputation_summary(
            before=previous.phospho,
            after=current.phospho,
        )
        return _StageDiagnosticsDefaults(
            dropped_row_ids=_resolve_dropped_row_ids(
                before=previous.phospho.index,
                after=current.phospho.index,
            ),
            imputed_row_ids=imputed_row_ids,
            imputed_cell_count=imputed_cell_count,
        )


class _StageDiagnosticsNormalizer:
    """Normalize stage-owned diagnostics into a deterministic immutable record."""

    def run(
        self,
        *,
        stage_key: str,
        raw: Mapping[str, object],
        defaults: _StageDiagnosticsDefaults,
    ) -> _NormalizedStageDiagnostics:
        dropped_row_ids = _coerce_string_tuple(
            raw.get("dropped_row_ids", defaults.dropped_row_ids),
            stage_key=stage_key,
            field_name="dropped_row_ids",
        )
        dropped_row_count = _coerce_int(
            raw.get("dropped_row_count"),
            stage_key=stage_key,
            field_name="dropped_row_count",
            default=len(dropped_row_ids),
        )

        imputed_row_ids = _coerce_string_tuple(
            raw.get("imputed_row_ids", defaults.imputed_row_ids),
            stage_key=stage_key,
            field_name="imputed_row_ids",
        )
        imputed_cell_count = _coerce_int(
            raw.get("imputed_cell_count"),
            stage_key=stage_key,
            field_name="imputed_cell_count",
            default=defaults.imputed_cell_count,
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

        diagnostics = _normalize_nested_diagnostics(stage_key=stage_key, raw=raw)

        return _NormalizedStageDiagnostics(
            dropped_row_ids=dropped_row_ids,
            dropped_row_count=dropped_row_count,
            imputed_cell_count=imputed_cell_count,
            imputed_row_ids=imputed_row_ids,
            notes=notes,
            diagnostics=MappingProxyType(diagnostics),
        )


def _normalize_nested_diagnostics(
    *,
    stage_key: str,
    raw: Mapping[str, object],
) -> dict[str, object]:
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
        return diagnostics
    if "diagnostics" in raw:
        raise DatasetBuildError(
            "dataset preprocessing stage diagnostics parse error: "
            f"stage={stage_key!r}, field='diagnostics', expected object, got "
            f"{nested_diagnostics!r} ({type(nested_diagnostics).__name__})"
        )

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
    return diagnostics


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
