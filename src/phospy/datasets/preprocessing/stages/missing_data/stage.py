"""Stage coordinator for missing-data preprocessing."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace

from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    PreprocessingStageResult,
    PreprocessingState,
    append_row_audit_records,
)
from phospy.datasets.preprocessing.report_rows import report_rows_from_row_audit_rows
from phospy.datasets.preprocessing.report_schema import PreprocessingRowAuditRow
from phospy.datasets.processing_state import JsonValue, MissingDataDiagnosticsV1
from phospy.errors.input import PhosPyInputError
from phospy.policy_models import MissingDataPolicy

from .audit import (
    build_knn_audit_records,
    build_minprob_audit_records,
    build_row_median_audit_records,
)
from .diagnostics import build_input_profile, build_missing_data_diagnostics
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
    )
    return PreprocessingStageResult(
        state=state,
        diagnostics=_stage_diagnostics_payload(
            dropped_row_ids=(),
            imputed_cell_count=0,
            imputed_row_ids=(),
            diagnostics=diagnostics,
        ),
    )


def _run_row_median_policy(
    *,
    state: PreprocessingState,
    input_profile: MissingDataInputProfile,
) -> PreprocessingStageResult:
    outcome = run_row_median_policy(state)
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
    )
    return _finalize_outcome(
        state=state,
        outcome=outcome,
        row_audit_records=row_audit_records,
        diagnostics=diagnostics,
    )


def _run_knn_policy(
    *,
    state: PreprocessingState,
    input_profile: MissingDataInputProfile,
) -> PreprocessingStageResult:
    outcome = run_knn_policy(state)
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
    )
    return _finalize_outcome(
        state=state,
        outcome=outcome,
        row_audit_records=row_audit_records,
        diagnostics=diagnostics,
    )


def _run_minprob_policy(
    *,
    state: PreprocessingState,
    input_profile: MissingDataInputProfile,
) -> PreprocessingStageResult:
    outcome = run_minprob_policy(state)
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
    )
    return _finalize_outcome(
        state=state,
        outcome=outcome,
        row_audit_records=row_audit_records,
        diagnostics=diagnostics,
    )


def _finalize_outcome(
    *,
    state: PreprocessingState,
    outcome: RowMedianPolicyOutcome | KnnPolicyOutcome | MinProbPolicyOutcome,
    row_audit_records: Sequence[PreprocessingRowAuditRow],
    diagnostics: Mapping[str, JsonValue],
) -> PreprocessingStageResult:
    next_state = append_row_audit_records(state, row_audit_records)
    return PreprocessingStageResult(
        state=replace(
            next_state,
            phospho=outcome.phospho,
            site_metadata=outcome.site_metadata,
        ),
        report_rows=report_rows_from_row_audit_rows(row_audit_records),
        diagnostics=_stage_diagnostics_payload(
            dropped_row_ids=outcome.dropped_row_ids,
            imputed_cell_count=outcome.imputed_cell_count,
            imputed_row_ids=outcome.imputed_row_ids,
            diagnostics=diagnostics,
        ),
    )


def _stage_diagnostics_payload(
    *,
    dropped_row_ids: tuple[str, ...],
    imputed_cell_count: int,
    imputed_row_ids: tuple[str, ...],
    diagnostics: Mapping[str, JsonValue],
) -> dict[str, object]:
    return {
        "dropped_row_ids": dropped_row_ids,
        "dropped_row_count": int(len(dropped_row_ids)),
        "imputed_cell_count": int(imputed_cell_count),
        "imputed_row_ids": imputed_row_ids,
        "notes": "stage executed",
        "diagnostics": diagnostics,
    }
