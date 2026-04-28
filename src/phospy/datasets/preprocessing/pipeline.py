"""Internal dataset preprocessing pipeline orchestration."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace

import pandas as pd

from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingReportRow,
    PreprocessingStage,
    PreprocessingStageExecution,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.datasets.preprocessing.report_rows import (
    validate_preprocessing_report_row,
)
from phospy.datasets.preprocessing.stages.comparisons import ComparisonsStage
from phospy.datasets.preprocessing.stages.intensity_transform import (
    IntensityTransformStage,
)
from phospy.datasets.preprocessing.stages.missing_data import MissingDataStage
from phospy.datasets.preprocessing.stages.normalisation import NormalisationStage
from phospy.datasets.preprocessing.stages.site_matrix import SiteMatrixStage
from phospy.datasets.preprocessing.stages.total_protein_correction import (
    TotalProteinCorrectionStage,
)
from phospy.errors.build import DatasetBuildError
from phospy.provenance.hashing import hash_table

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

_STAGE_LABEL_TO_PARAMETERS: dict[str, tuple[str, ...]] = {
    DATASET_PREPROCESSING_STAGE_NORMALISATION: (),
    DATASET_PREPROCESSING_STAGE_MISSING_DATA: (
        "missing_data_policy",
        "missing_data_min_observed_values",
    ),
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX: (
        "site_matrix_policy",
        "site_matrix_duplicate_site_policy",
        "site_matrix_missing_data_policy",
        "site_matrix_minimum_observed_values",
    ),
    DATASET_PREPROCESSING_STAGE_COMPARISONS: (
        "comparison_building_policy",
        "comparison_sample_group_column",
        "comparison_pairs",
    ),
}


class PreprocessingPipeline:
    """Apply ordered preprocessing stages for interpreted dataset input."""

    def __init__(
        self,
        *,
        stage_registry: tuple[PreprocessingStage, ...] | None = None,
    ) -> None:
        stages = stage_registry or (
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
            trace.append(
                PreprocessingStageExecution(
                    stage=stage_key,
                    operation=_resolve_stage_operation(
                        plan=previous.plan,
                        stage=stage_key,
                    ),
                    parameters=_resolve_stage_parameters(
                        plan=previous.plan,
                        stage=stage_key,
                    ),
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
                    imputed_cell_count=int(diagnostics["imputed_cell_count"]),
                    imputed_row_ids=tuple(diagnostics["imputed_row_ids"]),
                    notes=diagnostics["notes"],
                    diagnostics=dict(diagnostics["diagnostics"]),
                )
            )
        if report_rows:
            current = replace(current, report_rows=tuple(report_rows))
        return current, tuple(trace)


def _resolve_stage_parameters(
    *,
    plan: PreprocessingPlan,
    stage: str,
) -> dict[str, object]:
    if stage == DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION:
        identity = plan.total_protein_correction_identity_policy
        return {
            "total_protein_correction_policy": plan.total_protein_correction_policy,
            "identity_mode": identity.mode,
            "phosphosite_key": identity.phosphosite_key,
            "total_protein_key": identity.total_protein_key,
            "mapping_phosphosite_key": identity.mapping_phosphosite_key,
            "mapping_total_protein_key": identity.mapping_total_protein_key,
            "mapping_table_fingerprint": identity.mapping_table_fingerprint,
            "mapping_table_row_count": (
                None if identity.mapping_table is None else len(identity.mapping_table)
            ),
            "duplicate_policy": identity.duplicate_policy,
            "unmatched_policy": identity.unmatched_policy,
        }
    if stage == DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM:
        return {"pseudocount": float(plan.intensity_transform_pseudocount)}
    parameter_names = _STAGE_LABEL_TO_PARAMETERS.get(stage, ())
    return {name: getattr(plan, name) for name in parameter_names}


def _resolve_stage_operation(*, plan: PreprocessingPlan, stage: str) -> str:
    if stage == DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM:
        return plan.intensity_transform_policy
    if stage == DATASET_PREPROCESSING_STAGE_NORMALISATION:
        return plan.normalisation_policy
    if stage == DATASET_PREPROCESSING_STAGE_MISSING_DATA:
        return plan.missing_data_policy
    if stage == DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION:
        return str(plan.total_protein_correction_policy)
    if stage == DATASET_PREPROCESSING_STAGE_SITE_MATRIX:
        return plan.site_matrix_policy
    if stage == DATASET_PREPROCESSING_STAGE_COMPARISONS:
        return plan.comparison_building_policy
    return "unsupported_stage"


def _normalize_stage_diagnostics(
    *,
    raw: Mapping[str, object],
    previous: PreprocessingState,
    current: PreprocessingState,
) -> dict[str, object]:
    default_dropped_row_ids = _resolve_dropped_row_ids(
        before=previous.phospho.index,
        after=current.phospho.index,
    )
    dropped_row_ids = _coerce_string_tuple(
        raw.get("dropped_row_ids", default_dropped_row_ids)
    )
    dropped_row_count = int(raw.get("dropped_row_count", len(dropped_row_ids)))

    default_imputed_row_ids, default_imputed_cell_count = _resolve_imputation_summary(
        before=previous.phospho,
        after=current.phospho,
    )
    imputed_row_ids = _coerce_string_tuple(
        raw.get("imputed_row_ids", default_imputed_row_ids)
    )
    imputed_cell_count = int(raw.get("imputed_cell_count", default_imputed_cell_count))

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


__all__ = ["PreprocessingPipeline"]
