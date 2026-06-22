"""Group-aware phosphosite coverage filtering stage."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import replace
from typing import TypedDict

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.provenance.hashing import hash_table_tolerance
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
    PreprocessingStateTableKey,
    append_row_audit_records,
)
from phospy.science.datasets.preprocessing.report_rows import (
    report_rows_from_row_audit_rows,
)
from phospy.science.datasets.preprocessing.report_schema import PreprocessingRowAuditRow
from phospy.science.datasets.preprocessing.stage_contract import (
    PreprocessingStageContract,
)
from phospy.validation.datasets.group_coverage_filter import (
    GroupCoverageFilterMetadataValidator,
    ResolvedGroupCoverageFilterMetadata,
    require_numeric_group_coverage_matrix,
)

_REMOVAL_REASON = "insufficient finite coverage within configured sample groups"


class _ShapePayload(TypedDict):
    rows: int
    columns: int


class _CoverageDiagnostics(TypedDict):
    method: str
    group_column: str
    threshold_type: str
    min_finite_observations_per_group: int | None
    min_finite_fraction_per_group: float | None
    min_groups_passing_threshold: int
    input_feature_count: int
    retained_feature_count: int
    removed_feature_count: int
    removal_reason: str
    group_labels: list[str]
    group_sample_counts: dict[str, int]
    input_matrix_shape: _ShapePayload
    output_matrix_shape: _ShapePayload
    rows_dropped: bool
    dropped_row_ids: list[str]
    input_phospho_hash: str
    output_phospho_hash: str


class GroupCoverageFilterStage:
    """Filter phosphosite rows using configured group-wise finite coverage."""

    stage_key = DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER

    def __init__(
        self,
        *,
        metadata_validator: GroupCoverageFilterMetadataValidator | None = None,
    ) -> None:
        self._metadata_validator = (
            metadata_validator or GroupCoverageFilterMetadataValidator()
        )

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        if not state.plan.group_coverage_filter_enabled:
            diagnostics = _build_diagnostics(
                plan=state.plan,
                metadata=None,
                before=state.phospho,
                after=state.phospho,
                dropped_row_ids=(),
            )
            return PreprocessingStageResult(
                state=state,
                diagnostics=_stage_diagnostics_payload(
                    dropped_row_ids=(),
                    notes="group-aware coverage filter disabled by configuration",
                    diagnostics=diagnostics,
                ),
            )

        require_numeric_group_coverage_matrix(state.phospho)
        metadata = self._metadata_validator.run(
            phospho=state.phospho,
            sample_metadata=state.sample_metadata,
            group_column=state.plan.group_coverage_filter_group_column,
            min_groups_passing_threshold=(
                state.plan.group_coverage_filter_min_groups_passing_threshold
            ),
        )
        outcome = _evaluate_group_coverage(
            phospho=state.phospho,
            metadata=metadata,
            plan=state.plan,
        )
        retained_mask = outcome["retained_mask"]
        retained_count = int(retained_mask.sum())
        input_count = int(state.phospho.shape[0])
        if retained_count == 0:
            raise PhosPyInputError(
                "dataset preprocessing group_coverage_filter removed all "
                "phosphosites/features; "
                f"input_features={input_count}; retained=0; "
                f"removed={input_count}; group_column={metadata.group_column!r}; "
                f"{_threshold_summary(plan=state.plan)}; "
                f"reason={_REMOVAL_REASON}"
            )

        filtered_phospho = state.phospho.loc[retained_mask].copy(deep=True)
        filtered_site_metadata = state.site_metadata.loc[filtered_phospho.index].copy(
            deep=True
        )
        dropped_row_ids = tuple(
            str(row_id) for row_id in state.phospho.index[~retained_mask].tolist()
        )
        row_audit_records = _build_row_audit_records(
            plan=state.plan,
            metadata=metadata,
            outcome=outcome,
        )
        next_state = append_row_audit_records(state, row_audit_records)
        next_state = replace(
            next_state,
            phospho=filtered_phospho,
            site_metadata=filtered_site_metadata,
        )
        diagnostics = _build_diagnostics(
            plan=state.plan,
            metadata=metadata,
            before=state.phospho,
            after=filtered_phospho,
            dropped_row_ids=dropped_row_ids,
        )
        return PreprocessingStageResult(
            state=next_state,
            report_rows=report_rows_from_row_audit_rows(row_audit_records),
            diagnostics=_stage_diagnostics_payload(
                dropped_row_ids=dropped_row_ids,
                notes=(
                    "group-aware coverage filter retained "
                    f"{retained_count} of {input_count} features; "
                    f"removed {len(dropped_row_ids)}"
                ),
                diagnostics=diagnostics,
            ),
        )


class _CoverageOutcome(TypedDict):
    retained_mask: pd.Series
    finite_counts_by_group: dict[str, pd.Series]
    finite_fractions_by_group: dict[str, pd.Series]
    passing_groups_by_row: pd.Series


def _evaluate_group_coverage(
    *,
    phospho: pd.DataFrame,
    metadata: ResolvedGroupCoverageFilterMetadata,
    plan: PreprocessingPlan,
) -> _CoverageOutcome:
    numeric = phospho.astype("float64")
    finite_mask = pd.DataFrame(
        np.isfinite(numeric.to_numpy(dtype="float64", copy=False)),
        index=phospho.index.copy(),
        columns=phospho.columns.copy(),
    )
    group_passes: dict[str, pd.Series] = {}
    finite_counts_by_group: dict[str, pd.Series] = {}
    finite_fractions_by_group: dict[str, pd.Series] = {}
    threshold_count = plan.group_coverage_filter_min_finite_observations_per_group
    threshold_fraction = plan.group_coverage_filter_min_finite_fraction_per_group
    for group, samples in metadata.sample_order_by_group.items():
        group_counts = finite_mask.loc[:, list(samples)].sum(axis=1).astype("int64")
        finite_counts_by_group[group] = group_counts
        group_fraction = group_counts.astype("float64") / float(len(samples))
        finite_fractions_by_group[group] = group_fraction
        if threshold_count is not None:
            group_passes[group] = group_counts >= int(threshold_count)
        elif threshold_fraction is not None:
            group_passes[group] = group_fraction >= float(threshold_fraction)
        else:  # pragma: no cover - plan validation owns this branch.
            raise PhosPyInputError(
                "dataset preprocessing group_coverage_filter requires exactly one "
                "finite-observation threshold"
            )

    passing_groups_by_row = pd.DataFrame(group_passes, index=phospho.index).sum(axis=1)
    retained_mask = passing_groups_by_row >= int(
        plan.group_coverage_filter_min_groups_passing_threshold
    )
    return {
        "retained_mask": retained_mask.astype(bool),
        "finite_counts_by_group": finite_counts_by_group,
        "finite_fractions_by_group": finite_fractions_by_group,
        "passing_groups_by_row": passing_groups_by_row.astype("int64"),
    }


def _build_row_audit_records(
    *,
    plan: PreprocessingPlan,
    metadata: ResolvedGroupCoverageFilterMetadata,
    outcome: _CoverageOutcome,
) -> tuple[PreprocessingRowAuditRow, ...]:
    records: list[PreprocessingRowAuditRow] = []
    retained_mask = outcome["retained_mask"]
    retained_values = retained_mask.to_numpy(dtype=bool, copy=False)
    passing_group_values = outcome["passing_groups_by_row"].to_numpy(copy=False)
    dropped_positions = [
        position
        for position, retained in enumerate(retained_values)
        if not bool(retained)
    ]
    for position in dropped_positions:
        row_id = retained_mask.index[position]
        row_id_text = str(row_id)
        records.append(
            PreprocessingRowAuditRow(
                stage=DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER,
                action="dropped",
                reason=_REMOVAL_REASON,
                source_row_id=row_id_text,
                site_id=row_id_text,
                retained=False,
                retained_row_id=None,
                source_rows=(row_id_text,),
                retained_row=None,
                parameter_snapshot={
                    **_threshold_parameters(plan),
                    "group_column": metadata.group_column,
                    "groups_passing_threshold": int(passing_group_values[position]),
                    "finite_counts_by_group": _row_group_values(
                        outcome["finite_counts_by_group"],
                        row_position=position,
                        value_type="int",
                    ),
                    "finite_fractions_by_group": _row_group_values(
                        outcome["finite_fractions_by_group"],
                        row_position=position,
                        value_type="float",
                    ),
                    "stage_order": [str(stage) for stage in plan.stage_order],
                },
            )
        )
    return tuple(records)


def _row_group_values(
    values_by_group: Mapping[str, pd.Series],
    *,
    row_position: int,
    value_type: str,
) -> dict[str, int | float]:
    if value_type == "int":
        return {
            group: int(values.to_numpy(copy=False)[row_position])
            for group, values in values_by_group.items()
        }
    return {
        group: float(values.to_numpy(copy=False)[row_position])
        for group, values in values_by_group.items()
    }


def _build_diagnostics(
    *,
    plan: PreprocessingPlan,
    metadata: ResolvedGroupCoverageFilterMetadata | None,
    before: pd.DataFrame,
    after: pd.DataFrame,
    dropped_row_ids: tuple[str, ...],
) -> _CoverageDiagnostics:
    group_labels = [] if metadata is None else list(metadata.group_labels)
    group_sample_counts = (
        {}
        if metadata is None
        else {
            group: int(len(samples))
            for group, samples in metadata.sample_order_by_group.items()
        }
    )
    return {
        "method": "group_coverage_filter",
        "group_column": "" if metadata is None else str(metadata.group_column),
        "threshold_type": _threshold_type(plan),
        "min_finite_observations_per_group": (
            plan.group_coverage_filter_min_finite_observations_per_group
        ),
        "min_finite_fraction_per_group": (
            plan.group_coverage_filter_min_finite_fraction_per_group
        ),
        "min_groups_passing_threshold": int(
            plan.group_coverage_filter_min_groups_passing_threshold
        ),
        "input_feature_count": int(before.shape[0]),
        "retained_feature_count": int(after.shape[0]),
        "removed_feature_count": int(len(dropped_row_ids)),
        "removal_reason": _REMOVAL_REASON,
        "group_labels": group_labels,
        "group_sample_counts": group_sample_counts,
        "input_matrix_shape": _shape_payload(before),
        "output_matrix_shape": _shape_payload(after),
        "rows_dropped": bool(dropped_row_ids),
        "dropped_row_ids": list(dropped_row_ids),
        "input_phospho_hash": hash_table_tolerance(
            before,
            name="group_coverage_filter.input.phospho",
        ),
        "output_phospho_hash": hash_table_tolerance(
            after,
            name="group_coverage_filter.output.phospho",
        ),
    }


def _shape_payload(matrix: pd.DataFrame) -> _ShapePayload:
    return {
        "rows": int(matrix.shape[0]),
        "columns": int(matrix.shape[1]),
    }


def _stage_diagnostics_payload(
    *,
    dropped_row_ids: tuple[str, ...],
    notes: str,
    diagnostics: Mapping[str, object],
) -> dict[str, object]:
    return {
        "dropped_row_ids": dropped_row_ids,
        "dropped_row_count": int(len(dropped_row_ids)),
        "imputed_cell_count": 0,
        "imputed_row_ids": (),
        "notes": notes,
        "diagnostics": dict(diagnostics),
    }


def _threshold_type(plan: PreprocessingPlan) -> str:
    if plan.group_coverage_filter_min_finite_observations_per_group is not None:
        return "count"
    if plan.group_coverage_filter_min_finite_fraction_per_group is not None:
        return "fraction"
    return "none"


def _threshold_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    return {
        "threshold_type": _threshold_type(plan),
        "min_finite_observations_per_group": (
            plan.group_coverage_filter_min_finite_observations_per_group
        ),
        "min_finite_fraction_per_group": (
            plan.group_coverage_filter_min_finite_fraction_per_group
        ),
        "min_groups_passing_threshold": int(
            plan.group_coverage_filter_min_groups_passing_threshold
        ),
    }


def _threshold_summary(*, plan: PreprocessingPlan) -> str:
    parameters = _threshold_parameters(plan)
    if parameters["threshold_type"] == "count":
        return (
            "min_finite_observations_per_group="
            f"{parameters['min_finite_observations_per_group']}; "
            "min_groups_passing_threshold="
            f"{parameters['min_groups_passing_threshold']}"
        )
    return (
        "min_finite_fraction_per_group="
        f"{parameters['min_finite_fraction_per_group']}; "
        "min_groups_passing_threshold="
        f"{parameters['min_groups_passing_threshold']}"
    )


def _resolve_operation(plan: PreprocessingPlan) -> str:
    return _threshold_type(plan)


def _resolve_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    return {
        "group_column": plan.group_coverage_filter_group_column,
        **_threshold_parameters(plan),
    }


def _include_when(plan: PreprocessingPlan) -> bool:
    return bool(plan.group_coverage_filter_enabled)


GROUP_COVERAGE_FILTER_STAGE_CONTRACT = PreprocessingStageContract(
    stage_key=DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER,
    display_label=DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER,
    provenance_stage=DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER,
    operation_name=_resolve_operation,
    serialize_parameters=_resolve_parameters,
    consumed_input_tables=(
        PreprocessingStateTableKey.DATASET_PHOSPHO,
        PreprocessingStateTableKey.DATASET_SAMPLE_METADATA,
    ),
    produced_output_tables=(
        PreprocessingStateTableKey.DATASET_PHOSPHO,
        PreprocessingStateTableKey.DATASET_SITE_METADATA,
        PreprocessingStateTableKey.REPORT_ROW_AUDIT,
    ),
    stage_factory=GroupCoverageFilterStage,
    backend="numpy",
    include_when=_include_when,
    diagnostics_metadata={
        "known_diagnostics_fields": (
            "method",
            "group_column",
            "threshold_type",
            "min_finite_observations_per_group",
            "min_finite_fraction_per_group",
            "min_groups_passing_threshold",
            "input_feature_count",
            "retained_feature_count",
            "removed_feature_count",
            "removal_reason",
            "group_labels",
            "group_sample_counts",
            "input_matrix_shape",
            "output_matrix_shape",
            "rows_dropped",
            "dropped_row_ids",
            "input_phospho_hash",
            "output_phospho_hash",
        )
    },
)


__all__ = [
    "GROUP_COVERAGE_FILTER_STAGE_CONTRACT",
    "GroupCoverageFilterStage",
]
