"""Diagnostics parsing and coercion for preprocessing trace payloads."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingStageExecution,
)
from phospy.datasets.processing_state import SiteSequenceResolutionRowDiagnostic


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
        key: str,
        default: str | None,
    ) -> str | None:
        if diagnostics is None:
            return default
        value = diagnostics.get(key, default)
        if value is None:
            return None
        return str(value)

    @staticmethod
    def resolve_optional_bool(
        diagnostics: Mapping[str, object] | None,
        *,
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
        return default

    @staticmethod
    def resolve_optional_int(
        diagnostics: Mapping[str, object] | None,
        *,
        key: str,
        default: int,
    ) -> int:
        if diagnostics is None:
            return int(default)
        value = diagnostics.get(key)
        if value is None:
            return int(default)
        if isinstance(value, bool):
            return int(value)
        if isinstance(value, int):
            return int(value)
        if isinstance(value, float):
            return int(value)
        if isinstance(value, str):
            stripped = value.strip()
            return int(default) if stripped == "" else int(stripped)
        return int(default)

    @staticmethod
    def resolve_optional_nullable_int(
        diagnostics: Mapping[str, object] | None,
        *,
        key: str,
        default: int | None,
    ) -> int | None:
        if diagnostics is None:
            return default
        value = diagnostics.get(key, default)
        if value is None:
            return None
        if isinstance(value, bool):
            return int(value)
        if isinstance(value, int):
            return int(value)
        if isinstance(value, float):
            return int(value)
        if isinstance(value, str):
            stripped = value.strip()
            return default if stripped == "" else int(stripped)
        return default

    @staticmethod
    def resolve_optional_mapping_int(
        diagnostics: Mapping[str, object] | None,
        *,
        key: str,
    ) -> dict[str, int]:
        if diagnostics is None:
            return {}
        value = diagnostics.get(key)
        if not isinstance(value, Mapping):
            return {}
        resolved: dict[str, int] = {}
        for raw_key, raw_value in value.items():
            normalized_key = str(raw_key)
            if isinstance(raw_value, bool):
                resolved[normalized_key] = int(raw_value)
                continue
            if isinstance(raw_value, int):
                resolved[normalized_key] = int(raw_value)
                continue
            if isinstance(raw_value, float):
                resolved[normalized_key] = int(raw_value)
                continue
            if isinstance(raw_value, str):
                stripped = raw_value.strip()
                if stripped == "":
                    continue
                resolved[normalized_key] = int(stripped)
                continue
        return resolved

    @staticmethod
    def resolve_site_sequence_row_diagnostics(
        diagnostics: Mapping[str, object] | None,
    ) -> tuple[SiteSequenceResolutionRowDiagnostic, ...]:
        if diagnostics is None:
            return ()
        raw_value = diagnostics.get("row_diagnostics")
        if not isinstance(raw_value, list):
            return ()
        resolved: list[SiteSequenceResolutionRowDiagnostic] = []
        for raw_row in raw_value:
            if not isinstance(raw_row, Mapping):
                continue
            row_index = ProcessingTraceDiagnostics.resolve_optional_int(
                raw_row,
                key="row_index",
                default=-1,
            )
            if row_index < 0:
                continue
            row_id = ProcessingTraceDiagnostics.resolve_optional_string(
                raw_row,
                key="row_id",
                default=None,
            )
            if row_id is None:
                continue
            resolved.append(
                SiteSequenceResolutionRowDiagnostic(
                    row_index=row_index,
                    row_id=row_id,
                    site_id=ProcessingTraceDiagnostics.resolve_optional_string(
                        raw_row,
                        key="site_id",
                        default=None,
                    ),
                    status=(
                        ProcessingTraceDiagnostics.resolve_optional_string(
                            raw_row,
                            key="status",
                            default="unknown",
                        )
                        or "unknown"
                    ),
                    existing_site_sequence=(
                        ProcessingTraceDiagnostics.resolve_optional_string(
                            raw_row,
                            key="existing_site_sequence",
                            default=None,
                        )
                    ),
                    fasta_site_sequence=(
                        ProcessingTraceDiagnostics.resolve_optional_string(
                            raw_row,
                            key="fasta_site_sequence",
                            default=None,
                        )
                    ),
                    resolved_site_sequence=(
                        ProcessingTraceDiagnostics.resolve_optional_string(
                            raw_row,
                            key="resolved_site_sequence",
                            default=None,
                        )
                    ),
                    action=(
                        ProcessingTraceDiagnostics.resolve_optional_string(
                            raw_row,
                            key="action",
                            default="unknown",
                        )
                        or "unknown"
                    ),
                    reason=ProcessingTraceDiagnostics.resolve_optional_string(
                        raw_row,
                        key="reason",
                        default=None,
                    ),
                    conflict_policy=ProcessingTraceDiagnostics.resolve_optional_string(
                        raw_row,
                        key="conflict_policy",
                        default=None,
                    ),
                    resolver_version=ProcessingTraceDiagnostics.resolve_optional_string(
                        raw_row,
                        key="resolver_version",
                        default=None,
                    ),
                    fasta_source_path=ProcessingTraceDiagnostics.resolve_optional_string(
                        raw_row,
                        key="fasta_source_path",
                        default=None,
                    ),
                    fasta_sha256=ProcessingTraceDiagnostics.resolve_optional_string(
                        raw_row,
                        key="fasta_sha256",
                        default=None,
                    ),
                )
            )
        return tuple(resolved)

    @staticmethod
    def with_default_string(
        diagnostics: dict[str, object] | None,
        *,
        key: str,
        default: str | None,
    ) -> dict[str, object] | None:
        if diagnostics is None:
            if default is None:
                return None
            return {key: default}
        resolved = dict(diagnostics)
        value = resolved.get(key)
        if value is None and default is not None:
            resolved[key] = default
        return resolved

    @staticmethod
    def with_default_int(
        diagnostics: dict[str, object] | None,
        *,
        key: str,
        default: int,
    ) -> dict[str, object] | None:
        if diagnostics is None:
            return None
        resolved = dict(diagnostics)
        value = resolved.get(key)
        if value is None:
            resolved[key] = default
        return resolved


def _resolve_stage_diagnostics(
    *,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None,
    stage: str,
) -> dict[str, object] | None:
    if preprocessing_trace is None:
        return None
    for item in preprocessing_trace:
        if item.stage == stage:
            return dict(item.diagnostics)
    return None
