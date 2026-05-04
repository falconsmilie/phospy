"""Strict diagnostics parsing for preprocessing execution trace payloads."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Final, NoReturn

from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingStageExecution,
)
from phospy.datasets.processing_state import SiteSequenceResolutionRowDiagnostic
from phospy.errors.build import DatasetBuildError

_MISSING: Final = object()
_SITE_SEQUENCE_ALLOWED_FIELDS = frozenset(
    {
        "configured",
        "mode",
        "flank_size",
        "fasta_source_path",
        "fasta_source_label",
        "fasta_sha256",
        "resolver_version",
        "resolved_site_count",
        "unresolved_site_count",
        "unresolved_counts_by_reason",
        "filled_missing_count",
        "replaced_existing_count",
        "preserved_existing_count",
        "existing_sequence_conflict_count",
        "conflict_policy",
        "accession_column",
        "site_column",
        "row_status",
        "row_diagnostics",
    }
)
_SITE_SEQUENCE_ROW_DIAGNOSTIC_ALLOWED_FIELDS = frozenset(
    {
        "row_index",
        "row_id",
        "site_id",
        "status",
        "existing_site_sequence",
        "fasta_site_sequence",
        "resolved_site_sequence",
        "action",
        "reason",
        "conflict_policy",
        "resolver_version",
        "fasta_source_path",
        "fasta_sha256",
    }
)


@dataclass(frozen=True, slots=True)
class ProcessingTraceDiagnostics:
    """Parsed stage diagnostics extracted from preprocessing execution trace."""

    total_protein_correction: dict[str, object] | None
    missing_data: dict[str, object] | None
    site_sequence_resolution: dict[str, object] | None

    @classmethod
    def from_trace(
        cls,
        preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    ) -> ProcessingTraceDiagnostics:
        return cls(
            total_protein_correction=_resolve_stage_diagnostics(
                preprocessing_trace=preprocessing_trace,
                stage=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
            ),
            missing_data=_resolve_stage_diagnostics(
                preprocessing_trace=preprocessing_trace,
                stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
            ),
            site_sequence_resolution=_resolve_stage_diagnostics(
                preprocessing_trace=preprocessing_trace,
                stage=DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
            ),
        )

    @staticmethod
    def resolve_optional_string(
        diagnostics: Mapping[str, object] | None,
        *,
        stage: str,
        key: str,
        default: str | None,
    ) -> str | None:
        if diagnostics is None:
            return default
        value = diagnostics.get(key, _MISSING)
        if value is _MISSING:
            return default
        if value is None:
            return None
        if not isinstance(value, str):
            _raise_diagnostics_error(
                stage=stage,
                field=key,
                value=value,
                expected="string",
            )
        return value

    @staticmethod
    def resolve_optional_bool(
        diagnostics: Mapping[str, object] | None,
        *,
        stage: str,
        key: str,
        default: bool | None,
    ) -> bool | None:
        if diagnostics is None or key not in diagnostics:
            return default
        value = diagnostics.get(key)
        if value is None:
            return None
        if isinstance(value, bool):
            return value
        _raise_diagnostics_error(
            stage=stage,
            field=key,
            value=value,
            expected="bool",
        )
        raise AssertionError("unreachable")

    @staticmethod
    def resolve_optional_int(
        diagnostics: Mapping[str, object] | None,
        *,
        stage: str,
        key: str,
        default: int,
    ) -> int:
        if diagnostics is None:
            return default
        value = diagnostics.get(key, _MISSING)
        if value is _MISSING:
            return default
        if value is None:
            _raise_diagnostics_error(
                stage=stage,
                field=key,
                value=value,
                expected="int",
            )
        return _require_int(stage=stage, field=key, value=value)

    @staticmethod
    def resolve_optional_nullable_int(
        diagnostics: Mapping[str, object] | None,
        *,
        stage: str,
        key: str,
        default: int | None,
    ) -> int | None:
        if diagnostics is None:
            return default
        value = diagnostics.get(key, _MISSING)
        if value is _MISSING:
            return default
        if value is None:
            return None
        return _require_int(stage=stage, field=key, value=value)

    @staticmethod
    def resolve_optional_mapping_int(
        diagnostics: Mapping[str, object] | None,
        *,
        stage: str,
        key: str,
    ) -> dict[str, int]:
        if diagnostics is None:
            return {}
        value = diagnostics.get(key, _MISSING)
        if value is _MISSING:
            return {}
        if value is None:
            _raise_diagnostics_error(
                stage=stage,
                field=key,
                value=value,
                expected="object mapping string keys to int values",
            )
        if not isinstance(value, Mapping):
            _raise_diagnostics_error(
                stage=stage,
                field=key,
                value=value,
                expected="object mapping string keys to int values",
            )
        resolved: dict[str, int] = {}
        for raw_key, raw_value in value.items():
            if not isinstance(raw_key, str):
                _raise_diagnostics_error(
                    stage=stage,
                    field=f"{key}.<key>",
                    value=raw_key,
                    expected="string",
                )
            resolved[raw_key] = _require_int(
                stage=stage,
                field=f"{key}.{raw_key}",
                value=raw_value,
            )
        return resolved

    @staticmethod
    def validate_site_sequence_resolution_payload(
        diagnostics: Mapping[str, object] | None,
    ) -> None:
        if diagnostics is None:
            return
        stage = DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION
        unsupported = sorted(
            key for key in diagnostics if key not in _SITE_SEQUENCE_ALLOWED_FIELDS
        )
        if unsupported:
            raise DatasetBuildError(
                "dataset preprocessing diagnostics include unsupported fields for "
                f"stage {stage!r}: {', '.join(unsupported)}"
            )

    @staticmethod
    def resolve_site_sequence_row_diagnostics(
        diagnostics: Mapping[str, object] | None,
    ) -> tuple[SiteSequenceResolutionRowDiagnostic, ...]:
        if diagnostics is None:
            return ()
        stage = DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION
        raw_value = diagnostics.get("row_diagnostics")
        if raw_value is None:
            return ()
        if not isinstance(raw_value, list):
            _raise_diagnostics_error(
                stage=stage,
                field="row_diagnostics",
                value=raw_value,
                expected="array of row-diagnostic objects",
            )
            return ()
        resolved: list[SiteSequenceResolutionRowDiagnostic] = []
        for row_position, raw_row in enumerate(raw_value):
            if not isinstance(raw_row, Mapping):
                _raise_diagnostics_error(
                    stage=stage,
                    field=f"row_diagnostics[{row_position}]",
                    value=raw_row,
                    expected="object",
                )
            unsupported_row_fields = sorted(
                key
                for key in raw_row
                if key not in _SITE_SEQUENCE_ROW_DIAGNOSTIC_ALLOWED_FIELDS
            )
            if unsupported_row_fields:
                raise DatasetBuildError(
                    "dataset preprocessing diagnostics include unsupported fields for "
                    f"stage {stage!r}, row_diagnostics[{row_position}]: "
                    f"{', '.join(unsupported_row_fields)}"
                )
            row_index_raw = raw_row.get("row_index", _MISSING)
            if row_index_raw is _MISSING:
                _raise_diagnostics_error(
                    stage=stage,
                    field=f"row_diagnostics[{row_position}].row_index",
                    value=row_index_raw,
                    expected="required int >= 0",
                )
            row_index = _require_int(
                stage=stage,
                field=f"row_diagnostics[{row_position}].row_index",
                value=row_index_raw,
            )
            if row_index < 0:
                _raise_diagnostics_error(
                    stage=stage,
                    field=f"row_diagnostics[{row_position}].row_index",
                    value=row_index,
                    expected="int >= 0",
                )
            row_id_raw = raw_row.get("row_id", _MISSING)
            if row_id_raw is _MISSING:
                _raise_diagnostics_error(
                    stage=stage,
                    field=f"row_diagnostics[{row_position}].row_id",
                    value=row_id_raw,
                    expected="required string",
                )
            if not isinstance(row_id_raw, str):
                _raise_diagnostics_error(
                    stage=stage,
                    field=f"row_diagnostics[{row_position}].row_id",
                    value=row_id_raw,
                    expected="string",
                )
            resolved.append(
                SiteSequenceResolutionRowDiagnostic(
                    row_index=row_index,
                    row_id=row_id_raw,
                    site_id=ProcessingTraceDiagnostics.resolve_optional_string(
                        raw_row,
                        stage=stage,
                        key="site_id",
                        default=None,
                    ),
                    status=(
                        ProcessingTraceDiagnostics.resolve_optional_string(
                            raw_row,
                            stage=stage,
                            key="status",
                            default="unknown",
                        )
                        or "unknown"
                    ),
                    existing_site_sequence=(
                        ProcessingTraceDiagnostics.resolve_optional_string(
                            raw_row,
                            stage=stage,
                            key="existing_site_sequence",
                            default=None,
                        )
                    ),
                    fasta_site_sequence=(
                        ProcessingTraceDiagnostics.resolve_optional_string(
                            raw_row,
                            stage=stage,
                            key="fasta_site_sequence",
                            default=None,
                        )
                    ),
                    resolved_site_sequence=(
                        ProcessingTraceDiagnostics.resolve_optional_string(
                            raw_row,
                            stage=stage,
                            key="resolved_site_sequence",
                            default=None,
                        )
                    ),
                    action=(
                        ProcessingTraceDiagnostics.resolve_optional_string(
                            raw_row,
                            stage=stage,
                            key="action",
                            default="unknown",
                        )
                        or "unknown"
                    ),
                    reason=ProcessingTraceDiagnostics.resolve_optional_string(
                        raw_row,
                        stage=stage,
                        key="reason",
                        default=None,
                    ),
                    conflict_policy=ProcessingTraceDiagnostics.resolve_optional_string(
                        raw_row,
                        stage=stage,
                        key="conflict_policy",
                        default=None,
                    ),
                    resolver_version=ProcessingTraceDiagnostics.resolve_optional_string(
                        raw_row,
                        stage=stage,
                        key="resolver_version",
                        default=None,
                    ),
                    fasta_source_path=ProcessingTraceDiagnostics.resolve_optional_string(
                        raw_row,
                        stage=stage,
                        key="fasta_source_path",
                        default=None,
                    ),
                    fasta_sha256=ProcessingTraceDiagnostics.resolve_optional_string(
                        raw_row,
                        stage=stage,
                        key="fasta_sha256",
                        default=None,
                    ),
                )
            )
        return tuple(resolved)


def _resolve_stage_diagnostics(
    *,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    stage: str,
) -> dict[str, object] | None:
    if preprocessing_trace is None:
        return None
    for item in preprocessing_trace:
        if item.stage == stage:
            diagnostics = item.diagnostics
            if not isinstance(diagnostics, Mapping):
                _raise_diagnostics_error(
                    stage=stage,
                    field="diagnostics",
                    value=diagnostics,
                    expected="object",
                )
            resolved: dict[str, object] = {}
            for raw_key, raw_value in diagnostics.items():
                if not isinstance(raw_key, str):
                    _raise_diagnostics_error(
                        stage=stage,
                        field="diagnostics.<key>",
                        value=raw_key,
                        expected="string",
                    )
                resolved[raw_key] = raw_value
            return resolved
    return None


def _require_int(*, stage: str, field: str, value: object) -> int:
    if isinstance(value, bool):
        _raise_diagnostics_error(
            stage=stage,
            field=field,
            value=value,
            expected="int (bool is not accepted)",
        )
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        _raise_diagnostics_error(
            stage=stage,
            field=field,
            value=value,
            expected="int (floats are not accepted)",
        )
    _raise_diagnostics_error(
        stage=stage,
        field=field,
        value=value,
        expected="int",
    )
    raise AssertionError("unreachable")


def _raise_diagnostics_error(
    *,
    stage: str,
    field: str,
    value: object,
    expected: str,
) -> NoReturn:
    rendered_value = "<missing>" if value is _MISSING else repr(value)
    value_type = "missing" if value is _MISSING else type(value).__name__
    raise DatasetBuildError(
        "dataset preprocessing diagnostics parse error: "
        f"stage={stage!r}, field={field!r}, expected {expected}, "
        f"got {rendered_value} ({value_type})"
    )
