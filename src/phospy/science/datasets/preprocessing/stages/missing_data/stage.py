"""Stage coordinator for missing-data preprocessing."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets._processing_state.json_contracts import (
    MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1,
    V1_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS,
)
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
    PreprocessingStateTableKey,
    append_row_audit_records,
)
from phospy.science.datasets.preprocessing.policy_models import MissingDataPolicy
from phospy.science.datasets.preprocessing.report_rows import (
    report_rows_from_row_audit_rows,
)
from phospy.science.datasets.preprocessing.report_schema import PreprocessingRowAuditRow
from phospy.science.datasets.preprocessing.stage_contract import (
    DeterminismKind,
    PreprocessingStageContract,
    PreprocessingStageFactoryContext,
)
from phospy.science.datasets.processing_state import JsonValue, MissingDataDiagnosticsV1

from .audit import (
    build_knn_audit_records,
    build_minprob_audit_records,
    build_row_median_audit_records,
)
from .diagnostics import (
    build_input_profile,
    build_missing_data_diagnostics,
    hash_imputation_mask,
)
from .forbid import fail_if_forbid_policy_has_missing_values
from .knn import run_knn_policy
from .minprob import run_minprob_policy
from .models import (
    KnnPolicyOutcome,
    MinProbPolicyOutcome,
    MissingDataInputProfile,
    RowMedianPolicyOutcome,
)
from .row_median import run_row_median_policy


class MissingDataStage:
    """Apply the configured missing-data policy to phospho/site tables."""

    stage_key = DATASET_PREPROCESSING_STAGE_MISSING_DATA

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        input_profile = build_input_profile(state.phospho)
        policy = state.plan.missing_data_policy

        if policy is MissingDataPolicy.FORBID:
            return _run_forbid_policy(state=state, input_profile=input_profile)
        if policy is MissingDataPolicy.IMPUTE_ROW_MEDIAN:
            return _run_row_median_policy(state=state, input_profile=input_profile)
        if policy is MissingDataPolicy.IMPUTE_KNN:
            return _run_knn_policy(state=state, input_profile=input_profile)
        if policy is MissingDataPolicy.IMPUTE_MINPROB:
            return _run_minprob_policy(state=state, input_profile=input_profile)

        raise PhosPyInputError(
            "dataset build request preprocessing_config contains an unsupported "
            "missing_data.policy"
        )


def _run_forbid_policy(
    *,
    state: PreprocessingState,
    input_profile: MissingDataInputProfile,
) -> PreprocessingStageResult:
    policy = MissingDataPolicy.FORBID
    forbid_failure_diagnostics = MissingDataDiagnosticsV1.from_mapping(
        build_missing_data_diagnostics(
            missing_data_policy=policy.value,
            imputation_method_id="forbid",
            imputation_method_family="strict_rejection",
            input_missing_cell_count=input_profile.input_missing_cell_count,
            output_missing_cell_count=input_profile.input_missing_cell_count,
            imputed_cell_count=0,
            affected_row_ids=input_profile.affected_row_ids,
            affected_column_ids=input_profile.affected_column_ids,
            imputed_row_ids=(),
            imputed_column_ids=(),
            dropped_row_ids=(),
            random_seed=None,
            method_parameters={},
            matrix_scale_requirement=None,
            stage_order=state.plan.stage_order,
            missingness_mask_hash=input_profile.missingness_mask_hash,
            left_censored_assumption=False,
            rows_not_imputable=input_profile.affected_row_ids,
            row_medians_used={},
            per_column_distribution_parameters=None,
            dropped_rows_above_max_missing_fraction=(),
            neighbour_count=None,
            distance_metric=None,
            imputation_mask_hash=None,
        ),
        field_name="dataset preprocessing stage 'missing_data' diagnostics",
    )
    fail_if_forbid_policy_has_missing_values(
        state.phospho,
        diagnostics=forbid_failure_diagnostics,
    )
    diagnostics = build_missing_data_diagnostics(
        missing_data_policy=policy.value,
        imputation_method_id=None,
        imputation_method_family=None,
        input_missing_cell_count=input_profile.input_missing_cell_count,
        output_missing_cell_count=input_profile.input_missing_cell_count,
        imputed_cell_count=0,
        affected_row_ids=input_profile.affected_row_ids,
        affected_column_ids=input_profile.affected_column_ids,
        imputed_row_ids=(),
        imputed_column_ids=(),
        dropped_row_ids=(),
        random_seed=None,
        method_parameters={},
        matrix_scale_requirement=None,
        stage_order=state.plan.stage_order,
        missingness_mask_hash=input_profile.missingness_mask_hash,
        left_censored_assumption=None,
        rows_not_imputable=(),
        row_medians_used={},
        per_column_distribution_parameters=None,
        dropped_rows_above_max_missing_fraction=(),
        neighbour_count=None,
        distance_metric=None,
        imputation_mask_hash=None,
    )
    return PreprocessingStageResult(
        state=state,
        diagnostics=_stage_diagnostics_payload(
            dropped_row_ids=(),
            imputed_cell_count=0,
            imputed_row_ids=(),
            notes=(
                "missing_data policy='forbid'; complete matrix confirmed; "
                "no imputation applied"
            ),
            diagnostics=diagnostics,
        ),
    )


def _run_row_median_policy(
    *,
    state: PreprocessingState,
    input_profile: MissingDataInputProfile,
) -> PreprocessingStageResult:
    outcome = run_row_median_policy(state)
    imputation_mask_hash = hash_imputation_mask(outcome.imputed_mask)
    row_audit_records = build_row_median_audit_records(
        plan=state.plan,
        input_profile=input_profile,
        outcome=outcome,
    )
    diagnostics = build_missing_data_diagnostics(
        missing_data_policy=state.plan.missing_data_policy.value,
        imputation_method_id="row_median",
        imputation_method_family="deterministic_row_statistic",
        input_missing_cell_count=input_profile.input_missing_cell_count,
        output_missing_cell_count=outcome.output_missing_cell_count,
        imputed_cell_count=outcome.imputed_cell_count,
        affected_row_ids=input_profile.affected_row_ids,
        affected_column_ids=input_profile.affected_column_ids,
        imputed_row_ids=outcome.imputed_row_ids,
        imputed_column_ids=outcome.imputed_column_ids,
        dropped_row_ids=outcome.dropped_row_ids,
        random_seed=None,
        method_parameters={"min_observed_values": int(outcome.min_observed_values)},
        matrix_scale_requirement=None,
        stage_order=state.plan.stage_order,
        missingness_mask_hash=input_profile.missingness_mask_hash,
        left_censored_assumption=False,
        rows_not_imputable=outcome.rows_not_imputable,
        row_medians_used=outcome.row_medians_used,
        per_column_distribution_parameters=None,
        dropped_rows_above_max_missing_fraction=(),
        neighbour_count=None,
        distance_metric=None,
        imputation_mask_hash=imputation_mask_hash,
    )
    return _finalize_outcome(
        state=state,
        outcome=outcome,
        row_audit_records=row_audit_records,
        notes=_build_imputation_execution_note(
            policy=state.plan.missing_data_policy.value,
            imputed_cell_count=outcome.imputed_cell_count,
            imputed_row_ids=outcome.imputed_row_ids,
            dropped_row_ids=outcome.dropped_row_ids,
            output_missing_cell_count=outcome.output_missing_cell_count,
        ),
        diagnostics=diagnostics,
    )


def _run_knn_policy(
    *,
    state: PreprocessingState,
    input_profile: MissingDataInputProfile,
) -> PreprocessingStageResult:
    outcome = run_knn_policy(state)
    imputation_mask_hash = hash_imputation_mask(outcome.imputed_mask)
    row_audit_records = build_knn_audit_records(
        plan=state.plan,
        input_profile=input_profile,
        outcome=outcome,
    )
    diagnostics = build_missing_data_diagnostics(
        missing_data_policy=state.plan.missing_data_policy.value,
        imputation_method_id="knn",
        imputation_method_family="nearest_neighbour",
        input_missing_cell_count=input_profile.input_missing_cell_count,
        output_missing_cell_count=outcome.output_missing_cell_count,
        imputed_cell_count=outcome.imputed_cell_count,
        affected_row_ids=input_profile.affected_row_ids,
        affected_column_ids=input_profile.affected_column_ids,
        imputed_row_ids=outcome.imputed_row_ids,
        imputed_column_ids=outcome.imputed_column_ids,
        dropped_row_ids=outcome.dropped_row_ids,
        random_seed=None,
        method_parameters={
            "k": int(outcome.k),
            "distance": outcome.distance,
            "max_missing_fraction_per_row": float(outcome.max_missing_fraction_per_row),
        },
        matrix_scale_requirement=None,
        stage_order=state.plan.stage_order,
        missingness_mask_hash=input_profile.missingness_mask_hash,
        left_censored_assumption=False,
        rows_not_imputable=outcome.rows_not_imputable,
        row_medians_used={},
        per_column_distribution_parameters=None,
        dropped_rows_above_max_missing_fraction=outcome.dropped_row_ids,
        neighbour_count=int(outcome.k),
        distance_metric=outcome.distance,
        imputation_mask_hash=imputation_mask_hash,
    )
    return _finalize_outcome(
        state=state,
        outcome=outcome,
        row_audit_records=row_audit_records,
        notes=_build_imputation_execution_note(
            policy=state.plan.missing_data_policy.value,
            imputed_cell_count=outcome.imputed_cell_count,
            imputed_row_ids=outcome.imputed_row_ids,
            dropped_row_ids=outcome.dropped_row_ids,
            output_missing_cell_count=outcome.output_missing_cell_count,
        ),
        diagnostics=diagnostics,
    )


def _run_minprob_policy(
    *,
    state: PreprocessingState,
    input_profile: MissingDataInputProfile,
) -> PreprocessingStageResult:
    outcome = run_minprob_policy(state)
    imputation_mask_hash = hash_imputation_mask(outcome.imputed_mask)
    row_audit_records = build_minprob_audit_records(
        plan=state.plan,
        input_profile=input_profile,
        outcome=outcome,
    )
    diagnostics = build_missing_data_diagnostics(
        missing_data_policy=state.plan.missing_data_policy.value,
        imputation_method_id="minprob",
        imputation_method_family="left_censored_random",
        input_missing_cell_count=input_profile.input_missing_cell_count,
        output_missing_cell_count=outcome.output_missing_cell_count,
        imputed_cell_count=outcome.imputed_cell_count,
        affected_row_ids=input_profile.affected_row_ids,
        affected_column_ids=input_profile.affected_column_ids,
        imputed_row_ids=outcome.imputed_row_ids,
        imputed_column_ids=outcome.imputed_column_ids,
        dropped_row_ids=outcome.dropped_row_ids,
        random_seed=outcome.seed,
        method_parameters={
            "q": float(outcome.q),
            "width": float(outcome.width),
            "seed": int(outcome.seed),
            "max_missing_fraction_per_row": float(outcome.max_missing_fraction_per_row),
        },
        matrix_scale_requirement="log2",
        stage_order=state.plan.stage_order,
        missingness_mask_hash=input_profile.missingness_mask_hash,
        left_censored_assumption=True,
        rows_not_imputable=outcome.rows_not_imputable,
        row_medians_used={},
        per_column_distribution_parameters=outcome.per_column_distribution_parameters,
        dropped_rows_above_max_missing_fraction=outcome.dropped_row_ids,
        neighbour_count=None,
        distance_metric=None,
        imputation_mask_hash=imputation_mask_hash,
    )
    return _finalize_outcome(
        state=state,
        outcome=outcome,
        row_audit_records=row_audit_records,
        notes=_build_imputation_execution_note(
            policy=state.plan.missing_data_policy.value,
            imputed_cell_count=outcome.imputed_cell_count,
            imputed_row_ids=outcome.imputed_row_ids,
            dropped_row_ids=outcome.dropped_row_ids,
            output_missing_cell_count=outcome.output_missing_cell_count,
        ),
        diagnostics=diagnostics,
    )


def _finalize_outcome(
    *,
    state: PreprocessingState,
    outcome: RowMedianPolicyOutcome | KnnPolicyOutcome | MinProbPolicyOutcome,
    row_audit_records: Sequence[PreprocessingRowAuditRow],
    notes: str,
    diagnostics: Mapping[str, JsonValue],
) -> PreprocessingStageResult:
    next_state = append_row_audit_records(state, row_audit_records)
    observation_mask = _observation_mask_from_outcome(outcome)
    return PreprocessingStageResult(
        state=replace(
            next_state,
            phospho=outcome.phospho,
            site_metadata=outcome.site_metadata,
            imputation_observation_mask=observation_mask,
        ),
        report_rows=report_rows_from_row_audit_rows(row_audit_records),
        diagnostics=_stage_diagnostics_payload(
            dropped_row_ids=outcome.dropped_row_ids,
            imputed_cell_count=outcome.imputed_cell_count,
            imputed_row_ids=outcome.imputed_row_ids,
            notes=notes,
            diagnostics=diagnostics,
        ),
    )


def _observation_mask_from_outcome(
    outcome: RowMedianPolicyOutcome | KnnPolicyOutcome | MinProbPolicyOutcome,
) -> pd.DataFrame:
    imputed_mask = outcome.imputed_mask
    if not imputed_mask.index.equals(outcome.phospho.index):
        raise PhosPyInputError(
            "dataset preprocessing stage 'missing_data' produced an "
            "imputation mask with rows not aligned to phospho output"
        )
    if not imputed_mask.columns.equals(outcome.phospho.columns):
        raise PhosPyInputError(
            "dataset preprocessing stage 'missing_data' produced an "
            "imputation mask with columns not aligned to phospho output"
        )
    return (~imputed_mask.astype(bool)).copy(deep=True)


def _stage_diagnostics_payload(
    *,
    dropped_row_ids: tuple[str, ...],
    imputed_cell_count: int,
    imputed_row_ids: tuple[str, ...],
    notes: str,
    diagnostics: Mapping[str, JsonValue],
) -> dict[str, object]:
    return {
        "dropped_row_ids": dropped_row_ids,
        "dropped_row_count": int(len(dropped_row_ids)),
        "imputed_cell_count": int(imputed_cell_count),
        "imputed_row_ids": imputed_row_ids,
        "notes": notes,
        "diagnostics": diagnostics,
    }


def _build_imputation_execution_note(
    *,
    policy: str,
    imputed_cell_count: int,
    imputed_row_ids: tuple[str, ...],
    dropped_row_ids: tuple[str, ...],
    output_missing_cell_count: int,
) -> str:
    return (
        f"missing_data policy={policy!r}; "
        f"imputed_cells={int(imputed_cell_count)}; "
        f"imputed_rows={int(len(imputed_row_ids))}; "
        f"dropped_rows={int(len(dropped_row_ids))}; "
        f"output_missing_cells={int(output_missing_cell_count)}"
    )


def _resolve_operation(plan: PreprocessingPlan) -> str:
    return plan.missing_data_policy.value


def _resolve_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    return {
        "missing_data_policy": plan.missing_data_policy.value,
        "missing_data_min_observed_values": plan.missing_data_min_observed_values,
        "missing_data_q": plan.missing_data_q,
        "missing_data_width": plan.missing_data_width,
        "missing_data_seed": plan.missing_data_seed,
        "missing_data_k": plan.missing_data_k,
        "missing_data_distance": plan.missing_data_distance,
        "missing_data_max_missing_fraction_per_row": (
            plan.missing_data_max_missing_fraction_per_row
        ),
    }


def _resolve_determinism_kind(plan: PreprocessingPlan) -> DeterminismKind:
    if plan.missing_data_policy is MissingDataPolicy.IMPUTE_MINPROB:
        return DeterminismKind.SEEDED_STOCHASTIC
    return DeterminismKind.DETERMINISTIC


def _build_missing_data_stage(
    _context: PreprocessingStageFactoryContext,
) -> MissingDataStage:
    return MissingDataStage()


MISSING_DATA_STAGE_CONTRACT = PreprocessingStageContract(
    stage_key=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    display_label=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    provenance_stage=DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    operation_name=_resolve_operation,
    serialize_parameters=_resolve_parameters,
    consumed_input_tables=(
        PreprocessingStateTableKey.DATASET_PHOSPHO,
        PreprocessingStateTableKey.DATASET_SITE_METADATA,
    ),
    produced_output_tables=(
        PreprocessingStateTableKey.DATASET_PHOSPHO,
        PreprocessingStateTableKey.DATASET_SITE_METADATA,
        PreprocessingStateTableKey.DATASET_IMPUTATION_OBSERVATION_MASK,
        PreprocessingStateTableKey.REPORT_ROW_AUDIT,
    ),
    stage_factory=_build_missing_data_stage,
    backend="pandas",
    determinism_kind=_resolve_determinism_kind,
    diagnostics_metadata={
        "diagnostics_schema_version": MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1,
        "known_diagnostics_fields": tuple(
            sorted(V1_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS)
        ),
    },
)


__all__ = ["MISSING_DATA_STAGE_CONTRACT", "MissingDataStage"]
