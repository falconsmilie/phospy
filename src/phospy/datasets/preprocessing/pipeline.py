"""Internal dataset preprocessing pipeline orchestration."""

from __future__ import annotations

import pandas as pd

from phospy.api.configs import (
    DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY,
    DATASET_NORMALISATION_POLICY_NONE,
)
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingStage,
    PreprocessingStageExecution,
    PreprocessingState,
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

_STAGE_LABEL_TO_PARAMETERS: dict[str, tuple[str, ...]] = {
    DATASET_PREPROCESSING_STAGE_NORMALISATION: (),
    DATASET_PREPROCESSING_STAGE_MISSING_DATA: (
        "missing_data_policy",
        "missing_data_min_observed_values",
    ),
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION: (
        "total_protein_correction_policy",
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
            TotalProteinCorrectionStage(),
            SiteMatrixStage(),
            IntensityTransformStage(),
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
        for stage_key in current.plan.stage_order:
            stage = self._stages_by_key.get(stage_key)
            if stage is None:
                raise DatasetBuildError(
                    "dataset preprocessing plan references an unsupported stage: "
                    f"{stage_key}"
                )
            previous = current
            current = stage.run(current)
            input_hash = hash_table(
                previous.phospho,
                name=f"{stage_key}.input.phospho",
            )
            output_hash = hash_table(
                current.phospho,
                name=f"{stage_key}.output.phospho",
            )
            diagnostics = _resolve_stage_diagnostics(
                stage_key=stage_key,
                plan=previous.plan,
                previous=previous,
                current=current,
                input_hash=input_hash,
                output_hash=output_hash,
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
        return current, tuple(trace)


def _resolve_stage_parameters(
    *,
    plan: PreprocessingPlan,
    stage: str,
) -> dict[str, object]:
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
        return plan.total_protein_correction_policy
    if stage == DATASET_PREPROCESSING_STAGE_SITE_MATRIX:
        return plan.site_matrix_policy
    if stage == DATASET_PREPROCESSING_STAGE_COMPARISONS:
        return plan.comparison_building_policy
    return "unsupported_stage"


def _resolve_stage_diagnostics(
    *,
    stage_key: str,
    plan: PreprocessingPlan,
    previous: PreprocessingState,
    current: PreprocessingState,
    input_hash: str,
    output_hash: str,
) -> dict[str, object]:
    dropped_row_ids = _resolve_dropped_row_ids(
        before=previous.phospho.index,
        after=current.phospho.index,
    )
    details: dict[str, object] = {
        "dropped_row_ids": dropped_row_ids,
        "dropped_row_count": int(len(dropped_row_ids)),
        "imputed_cell_count": 0,
        "imputed_row_ids": (),
        "notes": "stage executed",
        "diagnostics": {},
    }

    if stage_key == DATASET_PREPROCESSING_STAGE_MISSING_DATA:
        imputed_row_ids, imputed_cell_count = _resolve_imputation_summary(
            before=previous.phospho,
            after=current.phospho,
        )
        details["imputed_cell_count"] = int(imputed_cell_count)
        details["imputed_row_ids"] = imputed_row_ids
        details["diagnostics"] = {
            "min_observed_values": plan.missing_data_min_observed_values,
            "imputed_row_ids": list(imputed_row_ids),
        }
        return details

    if stage_key == DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION:
        details["diagnostics"] = {
            "policy": plan.total_protein_correction_policy,
            "matched_rows": int(current.phospho.shape[0]),
            "total_table_hash": (
                None
                if previous.total is None
                else hash_table(
                    previous.total,
                    name="total_protein_correction.total",
                )
            ),
            "input_phospho_hash": input_hash,
            "output_phospho_hash": output_hash,
        }
        return details

    if stage_key == DATASET_PREPROCESSING_STAGE_SITE_MATRIX:
        site_matrix_diagnostics = current.phospho.attrs.get("site_matrix_provenance")
        if isinstance(site_matrix_diagnostics, dict):
            details["diagnostics"] = dict(site_matrix_diagnostics)
            dropped_ids = site_matrix_diagnostics.get("dropped_row_ids")
            if isinstance(dropped_ids, tuple):
                details["dropped_row_ids"] = dropped_ids
                details["dropped_row_count"] = int(len(dropped_ids))
        if current.duplicate_site_resolution is not None:
            details["diagnostics"] = {
                **dict(details["diagnostics"]),
                "duplicate_site_decisions": _records_from_frame(
                    current.duplicate_site_resolution
                ),
            }
        details["diagnostics"] = {
            **dict(details["diagnostics"]),
            "final_constructed_site_ids": [
                str(site_id) for site_id in current.phospho.index.tolist()
            ],
        }
        return details

    if stage_key == DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM:
        affected_matrices = ["phospho"]
        diagnostics: dict[str, object] = {
            "policy": plan.intensity_transform_policy,
            "pseudocount": float(plan.intensity_transform_pseudocount),
            "affected_matrices": affected_matrices,
            "input_phospho_hash": input_hash,
            "output_phospho_hash": output_hash,
        }
        if previous.total is not None and current.total is not None:
            diagnostics["input_total_hash"] = hash_table(
                previous.total,
                name="intensity_transform.input.total",
            )
            diagnostics["output_total_hash"] = hash_table(
                current.total,
                name="intensity_transform.output.total",
            )
            affected_matrices.append("total")
        details["diagnostics"] = diagnostics
        return details

    if stage_key == DATASET_PREPROCESSING_STAGE_NORMALISATION:
        diagnostics = {
            "policy": plan.normalisation_policy,
            "affected_columns": [
                str(column) for column in current.phospho.columns.tolist()
            ],
            "input_phospho_hash": input_hash,
            "output_phospho_hash": output_hash,
        }
        if plan.normalisation_policy != DATASET_NORMALISATION_POLICY_NONE:
            diagnostics["note"] = (
                "quantile normalisation used"
                if plan.normalisation_policy == "quantile"
                else "median centering used"
            )
        details["diagnostics"] = diagnostics
        return details

    if stage_key == DATASET_PREPROCESSING_STAGE_COMPARISONS:
        details["diagnostics"] = {
            "policy": plan.comparison_building_policy,
            "sample_group_column": plan.comparison_sample_group_column,
            "resolved_comparison_pairs": _resolve_comparison_pairs(current),
            "group_labels": _resolve_group_labels(current),
            "output_comparison_hash": (
                None
                if current.comparisons is None
                else hash_table(
                    current.comparisons,
                    name="comparisons.output.table",
                )
            ),
        }
        return details

    if (
        stage_key == DATASET_PREPROCESSING_STAGE_MISSING_DATA
        and plan.missing_data_policy == "forbid"
    ):
        details["notes"] = "stage executed (forbid policy passthrough)"
    if (
        stage_key == DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM
        and plan.intensity_transform_policy
        == DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY
    ):
        details["notes"] = "stage executed (identity transform passthrough)"
    return details


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


def _resolve_comparison_pairs(state: PreprocessingState) -> list[tuple[str, str]]:
    pair_stats = state.comparison_pair_stats
    if pair_stats is None or pair_stats.empty:
        return []
    pairs = pair_stats.loc[:, ["left_group", "right_group"]].drop_duplicates()
    return [
        (str(left), str(right))
        for left, right in pairs.itertuples(index=False, name=None)
    ]


def _resolve_group_labels(state: PreprocessingState) -> list[str]:
    group_stats = state.comparison_group_stats
    if group_stats is None or group_stats.empty:
        return []
    labels = group_stats.loc[:, "group"].astype(str).drop_duplicates().tolist()
    return [str(label) for label in labels]


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


__all__ = ["PreprocessingPipeline"]
