"""Stage coordinator for missing-data preprocessing."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.configs.preprocessing import (
    DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICY_ERROR,
    DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICY_VERSION,
)
from phospy.science.datasets._processing_state.json_contracts import (
    MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1,
    MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V2,
    V1_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS,
    V2_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS,
)
from phospy.science.datasets.preprocessing.imputation_scale_policy import (
    imputation_input_scale_kind,
)
from phospy.science.datasets.preprocessing.missing_data_mask_hashing import (
    hash_group_aware_imputation_mask,
    hash_group_aware_mechanism_mask,
    hash_imputation_mask,
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
from phospy.science.datasets.processing_state import (
    JsonValue,
    MissingDataDiagnosticsV1,
    MissingDataDiagnosticsV2,
)
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
)
from phospy.science.transformations.quantitative_contracts import (
    ALL_QUANTITATIVE_MEANINGS,
    NegativeDomainPolicy,
    QuantitativeEvidenceRequirement,
    QuantitativeInformationLossKind,
    QuantitativeOperationContract,
    QuantitativeReversibilityKind,
    preserve_meaning_transition,
    preserve_quantitative_contract,
    preserve_scale_transition,
)

from .audit import (
    build_group_aware_audit_records,
    build_knn_audit_records,
    build_minprob_audit_records,
    build_row_median_audit_records,
)
from .diagnostics import (
    build_input_profile,
    build_missing_data_diagnostics,
)
from .forbid import fail_if_forbid_policy_has_missing_values
from .group_aware_routing import route_group_aware_missingness
from .knn import impute_knn_targets, run_knn_policy
from .minprob import impute_minprob_targets, run_minprob_policy
from .models import (
    GroupAwarePolicyOutcome,
    GroupAwareRoutingOutcome,
    GroupMissingnessClassification,
    KnnPolicyOutcome,
    MinProbPolicyOutcome,
    MissingDataInputProfile,
    RowImputationRecord,
    RowMedianPolicyOutcome,
)
from .row_median import run_row_median_policy

KNN_COLUMN_MEAN_FALLBACK_CAVEAT_CODE = "knn_column_mean_fallback_used"


class MissingDataStage:
    """Apply the configured missing-data policy to phospho/site tables."""

    stage_key = DATASET_PREPROCESSING_STAGE_MISSING_DATA

    def validate_before_quantitative_contract(
        self,
        state: PreprocessingState,
    ) -> None:
        del state
        return None

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        policy = state.plan.missing_data_policy
        input_profile = build_input_profile(state.phospho)

        if policy is MissingDataPolicy.IMPUTE_GROUP_AWARE:
            return _run_group_aware_policy(
                state=state,
                input_profile=input_profile,
            )
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
            imputation_input_scale=None,
            imputation_input_scale_source=None,
            imputation_operation_order=None,
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
        imputation_input_scale=None,
        imputation_input_scale_source=None,
        imputation_operation_order=None,
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
    imputation_input_scale = _required_imputation_input_scale_value(state.plan)
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
        method_parameters={
            "min_observed_values": int(outcome.min_observed_values),
            "input_scale": imputation_input_scale,
            "imputation_operation_order": (
                state.plan.missing_data_imputation_operation_order
            ),
        },
        matrix_scale_requirement=None,
        imputation_input_scale=imputation_input_scale,
        imputation_input_scale_source=state.plan.missing_data_input_scale_source,
        imputation_operation_order=state.plan.missing_data_imputation_operation_order,
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
    imputation_input_scale = _required_imputation_input_scale_value(state.plan)
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
            "no_overlap_policy": outcome.no_overlap_policy,
            "no_overlap_policy_version": int(outcome.no_overlap_policy_version),
            "input_scale": imputation_input_scale,
            "imputation_operation_order": (
                state.plan.missing_data_imputation_operation_order
            ),
        },
        matrix_scale_requirement=None,
        imputation_input_scale=imputation_input_scale,
        imputation_input_scale_source=state.plan.missing_data_input_scale_source,
        imputation_operation_order=state.plan.missing_data_imputation_operation_order,
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
        knn_no_overlap_policy=outcome.no_overlap_policy,
        knn_no_overlap_policy_version=outcome.no_overlap_policy_version,
        knn_nearest_neighbour_imputed_cell_count=(
            outcome.nearest_neighbour_imputed_cell_count
        ),
        knn_nearest_neighbour_imputed_row_ids=(
            outcome.nearest_neighbour_imputed_row_ids
        ),
        knn_nearest_neighbour_imputed_column_ids=(
            outcome.nearest_neighbour_imputed_column_ids
        ),
        knn_column_mean_fallback_imputed_cell_count=(
            outcome.column_mean_fallback_imputed_cell_count
        ),
        knn_column_mean_fallback_row_ids=outcome.column_mean_fallback_row_ids,
        knn_column_mean_fallback_column_ids=outcome.column_mean_fallback_column_ids,
        knn_nearest_neighbour_imputation_mask_hash=(
            outcome.nearest_neighbour_imputation_mask_hash
        ),
        knn_column_mean_fallback_imputation_mask_hash=(
            outcome.column_mean_fallback_imputation_mask_hash
        ),
        knn_fully_column_mean_fallback_row_ids=(
            outcome.fully_column_mean_fallback_row_ids
        ),
        diagnostic_caveat_codes=(
            (KNN_COLUMN_MEAN_FALLBACK_CAVEAT_CODE,)
            if outcome.column_mean_fallback_imputed_cell_count > 0
            else ()
        ),
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
    imputation_input_scale = _required_imputation_input_scale_value(state.plan)
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
            "input_scale": imputation_input_scale,
            "imputation_operation_order": (
                state.plan.missing_data_imputation_operation_order
            ),
        },
        matrix_scale_requirement="log2",
        imputation_input_scale=imputation_input_scale,
        imputation_input_scale_source=state.plan.missing_data_input_scale_source,
        imputation_operation_order=state.plan.missing_data_imputation_operation_order,
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


def _run_group_aware_policy(
    *,
    state: PreprocessingState,
    input_profile: MissingDataInputProfile,
) -> PreprocessingStageResult:
    (
        group_column,
        min_partial_observed_fraction,
        min_reference_observed_fraction,
        q,
        width,
        seed,
        k,
        distance,
        no_overlap_policy,
    ) = _require_group_aware_parameters(state.plan)

    routing = route_group_aware_missingness(
        phospho=state.phospho,
        sample_metadata=state.sample_metadata,
        group_column=group_column,
        min_partial_observed_fraction=min_partial_observed_fraction,
        min_reference_observed_fraction=min_reference_observed_fraction,
    )
    outcome = execute_group_aware_imputation(
        state=state,
        routing=routing,
        q=q,
        width=width,
        seed=seed,
        k=k,
        distance=distance,
        no_overlap_policy=no_overlap_policy,
    )
    row_audit_records = build_group_aware_audit_records(
        plan=state.plan,
        input_profile=input_profile,
        outcome=outcome,
    )
    imputation_input_scale = _required_imputation_input_scale_value(state.plan)
    diagnostics = build_missing_data_diagnostics(
        missing_data_policy=state.plan.missing_data_policy.value,
        imputation_method_id="group_aware_knn_minprob",
        imputation_method_family="group_aware_mixed_mechanism",
        input_missing_cell_count=input_profile.input_missing_cell_count,
        output_missing_cell_count=0,
        imputed_cell_count=outcome.imputed_cell_count,
        affected_row_ids=input_profile.affected_row_ids,
        affected_column_ids=input_profile.affected_column_ids,
        imputed_row_ids=outcome.imputed_row_ids,
        imputed_column_ids=outcome.imputed_column_ids,
        dropped_row_ids=routing.dropped_row_ids,
        random_seed=seed,
        method_parameters={
            "group_column": group_column,
            "min_partial_observed_fraction": min_partial_observed_fraction,
            "min_reference_observed_fraction": min_reference_observed_fraction,
            "k": k,
            "distance": distance,
            "no_overlap_policy": no_overlap_policy,
            "no_overlap_policy_version": (
                DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICY_VERSION
            ),
            "q": q,
            "width": width,
            "seed": seed,
            "knn_target_cell_count": outcome.knn_target_cell_count,
            "minprob_target_cell_count": outcome.minprob_target_cell_count,
            "knn_target_mask_hash": outcome.knn_target_mask_hash,
            "minprob_target_mask_hash": outcome.minprob_target_mask_hash,
            "mechanism_input": "original_retained_matrix",
            "resolved_group_samples": {
                group: list(samples)
                for group, samples in routing.resolved_groups.sample_order_by_group.items()
            },
            "input_scale": imputation_input_scale,
            "imputation_operation_order": (
                state.plan.missing_data_imputation_operation_order
            ),
        },
        matrix_scale_requirement="log2",
        imputation_input_scale=imputation_input_scale,
        imputation_input_scale_source=state.plan.missing_data_input_scale_source,
        imputation_operation_order=state.plan.missing_data_imputation_operation_order,
        stage_order=state.plan.stage_order,
        missingness_mask_hash=input_profile.missingness_mask_hash,
        left_censored_assumption=True,
        rows_not_imputable=routing.dropped_row_ids,
        row_medians_used={},
        per_column_distribution_parameters=(outcome.per_column_distribution_parameters),
        dropped_rows_above_max_missing_fraction=(),
        neighbour_count=k,
        distance_metric=distance,
        imputation_mask_hash=outcome.imputation_mask_hash,
        knn_no_overlap_policy=no_overlap_policy,
        knn_no_overlap_policy_version=(
            DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICY_VERSION
        ),
    )
    unsupported_partial_row_ids = tuple(
        record.row_id
        for record in routing.dropped_row_reasons
        if any(
            classification is GroupMissingnessClassification.UNSUPPORTED_PARTIAL
            for classification in record.reasons_by_group.values()
        )
    )
    unsupported_absence_row_ids = tuple(
        record.row_id
        for record in routing.dropped_row_reasons
        if any(
            classification is GroupMissingnessClassification.UNSUPPORTED_FULLY_MISSING
            for classification in record.reasons_by_group.values()
        )
    )
    sample_group_by_column = {
        str(sample): group
        for group, samples in routing.resolved_groups.sample_order_by_group.items()
        for sample in samples
    }
    routing_facts_by_row = {
        row_id: tuple(fact for fact in routing.group_facts if fact.row_id == row_id)
        for row_id in (*routing.retained_row_ids, *routing.dropped_row_ids)
    }
    diagnostics["diagnostics_schema_version"] = (
        MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V2
    )
    diagnostics["group_aware"] = {
        "group_column": routing.resolved_groups.group_column,
        "observed_group_sizes": {
            group: len(samples)
            for group, samples in routing.resolved_groups.sample_order_by_group.items()
        },
        "min_partial_observed_fraction": routing.min_partial_observed_fraction,
        "min_reference_observed_fraction": routing.min_reference_observed_fraction,
        "retained_row_count": len(routing.retained_row_ids),
        "dropped_unsupported_row_count": len(routing.dropped_row_ids),
        "knn_routed_cell_count": outcome.knn_target_cell_count,
        "minprob_routed_cell_count": outcome.minprob_target_cell_count,
        "knn_imputed_cell_count": outcome.knn_imputed_cell_count,
        "minprob_imputed_cell_count": outcome.minprob_imputed_cell_count,
        "knn_target_mask_hash": outcome.knn_target_mask_hash,
        "minprob_target_mask_hash": outcome.minprob_target_mask_hash,
        "knn_imputation_mask_hash": outcome.knn_imputation_mask_hash,
        "minprob_imputation_mask_hash": outcome.minprob_imputation_mask_hash,
        "unsupported_partial_group_count": routing.unsupported_partial_group_count,
        "unsupported_absence_group_count": (
            routing.unsupported_fully_missing_group_count
        ),
        "unsupported_partial_row_ids": list(unsupported_partial_row_ids),
        "unsupported_absence_row_ids": list(unsupported_absence_row_ids),
        "routed_rows": [
            {
                "row_id": row.row_id,
                "knn_imputed_columns": list(row.nearest_neighbour_imputed_columns),
                "minprob_imputed_columns": [
                    column
                    for column in row.imputed_columns
                    if column not in set(row.nearest_neighbour_imputed_columns)
                ],
                "knn_group_labels": list(
                    dict.fromkeys(
                        sample_group_by_column[column]
                        for column in row.nearest_neighbour_imputed_columns
                    )
                ),
                "minprob_group_labels": list(
                    dict.fromkeys(
                        sample_group_by_column[column]
                        for column in row.imputed_columns
                        if column not in set(row.nearest_neighbour_imputed_columns)
                    )
                ),
            }
            for row in outcome.imputed_rows
        ],
        "rejected_rows": [
            {
                "row_id": record.row_id,
                "unsupported_partial_groups": [
                    group
                    for group, classification in record.reasons_by_group.items()
                    if classification
                    is GroupMissingnessClassification.UNSUPPORTED_PARTIAL
                ],
                "unsupported_absence_groups": [
                    group
                    for group, classification in record.reasons_by_group.items()
                    if classification
                    is GroupMissingnessClassification.UNSUPPORTED_FULLY_MISSING
                ],
                "observed_finite_count_by_group": {
                    fact.group_label: fact.observed_finite_count
                    for fact in routing_facts_by_row[record.row_id]
                },
                "observed_fraction_by_group": {
                    fact.group_label: fact.observed_fraction
                    for fact in routing_facts_by_row[record.row_id]
                },
            }
            for record in routing.dropped_row_reasons
        ],
        "minprob_left_censored_assumption": True,
        "route_categories": [
            "partial_observation_knn",
            "asymmetric_absence_minprob",
            "unsupported_partial",
            "unsupported_absence",
        ],
        "mechanism_input": "original_retained_matrix",
    }
    diagnostics = MissingDataDiagnosticsV2.from_mapping(
        diagnostics,
        field_name="dataset preprocessing stage 'missing_data' diagnostics",
    ).to_payload()
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


def execute_group_aware_imputation(
    *,
    state: PreprocessingState,
    routing: GroupAwareRoutingOutcome,
    q: float,
    width: float,
    seed: int,
    k: int,
    distance: str,
    no_overlap_policy: str,
) -> GroupAwarePolicyOutcome:
    """Run both routed mechanisms and return their dedicated typed outcome."""

    original_retained = routing.retain_rows(state.phospho)
    retained_site_metadata = state.site_metadata.loc[original_retained.index].copy(
        deep=True
    )
    knn_target_mask = routing.knn_target_mask.loc[original_retained.index].copy(
        deep=True
    )
    minprob_target_mask = routing.minprob_target_mask.loc[original_retained.index].copy(
        deep=True
    )
    knn_result = impute_knn_targets(
        original_retained,
        target_mask=knn_target_mask,
        k=k,
        distance=distance,
        no_overlap_policy=no_overlap_policy,
        policy_name=MissingDataPolicy.IMPUTE_GROUP_AWARE.value,
    )
    minprob_result = impute_minprob_targets(
        original_retained,
        target_mask=minprob_target_mask,
        q=q,
        width=width,
        seed=seed,
    )

    knn_imputed_mask = knn_result.imputed_mask
    minprob_imputed_mask = minprob_result.imputed_mask
    final_imputed_mask = knn_imputed_mask | minprob_imputed_mask
    merged_values = original_retained.to_numpy(dtype=float, copy=True, na_value=np.nan)
    knn_values = knn_result.phospho.to_numpy(dtype=float, copy=False, na_value=np.nan)
    minprob_values = minprob_result.phospho.to_numpy(
        dtype=float, copy=False, na_value=np.nan
    )
    knn_mask_values = knn_imputed_mask.to_numpy(dtype=bool, copy=False)
    minprob_mask_values = minprob_imputed_mask.to_numpy(dtype=bool, copy=False)
    merged_values[knn_mask_values] = knn_values[knn_mask_values]
    merged_values[minprob_mask_values] = minprob_values[minprob_mask_values]
    merged = pd.DataFrame(
        merged_values,
        index=original_retained.index.copy(),
        columns=original_retained.columns.copy(),
    )

    _validate_group_aware_merge(
        original_retained=original_retained,
        final_matrix=merged,
        knn_target_mask=knn_target_mask,
        minprob_target_mask=minprob_target_mask,
        knn_imputed_mask=knn_imputed_mask,
        minprob_imputed_mask=minprob_imputed_mask,
        final_imputed_mask=final_imputed_mask,
        knn_column_mean_fallback_mask=(knn_result.column_mean_fallback_imputed_mask),
    )

    imputed_row_ids = _mask_row_ids(final_imputed_mask)
    imputed_column_ids = _mask_column_ids(final_imputed_mask)
    imputed_rows = tuple(
        RowImputationRecord(
            row_id=str(row_id),
            imputed_columns=tuple(
                str(column)
                for column in merged.columns[final_imputed_mask.loc[row_id]].tolist()
            ),
            imputed_cell_count=int(final_imputed_mask.loc[row_id].sum()),
            nearest_neighbour_imputed_columns=tuple(
                str(column)
                for column in merged.columns[knn_imputed_mask.loc[row_id]].tolist()
            ),
        )
        for row_id in merged.index[
            final_imputed_mask.any(axis=1).to_numpy(dtype=bool, copy=False)
        ]
    )
    return GroupAwarePolicyOutcome(
        phospho=merged,
        site_metadata=retained_site_metadata,
        imputed_mask=final_imputed_mask,
        knn_target_mask=knn_target_mask,
        minprob_target_mask=minprob_target_mask,
        knn_imputed_mask=knn_imputed_mask,
        minprob_imputed_mask=minprob_imputed_mask,
        routing=routing,
        q=q,
        width=width,
        seed=seed,
        k=k,
        distance=distance,
        no_overlap_policy=no_overlap_policy,
        per_column_distribution_parameters=(
            minprob_result.per_column_distribution_parameters
        ),
        dropped_row_ids=routing.dropped_row_ids,
        imputed_cell_count=int(final_imputed_mask.to_numpy().sum()),
        imputed_row_ids=imputed_row_ids,
        imputed_column_ids=imputed_column_ids,
        output_missing_cell_count=0,
        rows_not_imputable=routing.dropped_row_ids,
        imputed_rows=imputed_rows,
        knn_target_cell_count=int(knn_target_mask.to_numpy().sum()),
        minprob_target_cell_count=int(minprob_target_mask.to_numpy().sum()),
        knn_imputed_cell_count=int(knn_imputed_mask.to_numpy().sum()),
        minprob_imputed_cell_count=int(minprob_imputed_mask.to_numpy().sum()),
        knn_target_mask_hash=hash_group_aware_mechanism_mask(
            knn_target_mask, mechanism="knn", mask_kind="target"
        ),
        minprob_target_mask_hash=hash_group_aware_mechanism_mask(
            minprob_target_mask, mechanism="minprob", mask_kind="target"
        ),
        knn_imputation_mask_hash=hash_group_aware_mechanism_mask(
            knn_imputed_mask, mechanism="knn", mask_kind="imputation"
        ),
        minprob_imputation_mask_hash=hash_group_aware_mechanism_mask(
            minprob_imputed_mask, mechanism="minprob", mask_kind="imputation"
        ),
        imputation_mask_hash=hash_group_aware_imputation_mask(final_imputed_mask),
    )


def _require_group_aware_parameters(
    plan: PreprocessingPlan,
) -> tuple[str, float, float, float, float, int, int, str, str]:
    group_column_value = plan.missing_data_group_column
    min_partial_value = plan.missing_data_min_partial_observed_fraction
    min_reference_value = plan.missing_data_min_reference_observed_fraction
    q_value = plan.missing_data_q
    width_value = plan.missing_data_width
    seed_value = plan.missing_data_seed
    k_value = plan.missing_data_k
    distance_value = plan.missing_data_distance
    no_overlap_policy_value = plan.missing_data_no_overlap_policy
    values = (
        group_column_value,
        min_partial_value,
        min_reference_value,
        q_value,
        width_value,
        seed_value,
        k_value,
        distance_value,
        no_overlap_policy_value,
    )
    if any(value is None for value in values):
        raise PhosPyInputError(
            "dataset build request preprocessing_config.missing_data.policy="
            "'impute_group_aware' requires group routing, KNN, and MinProb parameters"
        )
    assert group_column_value is not None
    assert min_partial_value is not None
    assert min_reference_value is not None
    assert q_value is not None
    assert width_value is not None
    assert seed_value is not None
    assert k_value is not None
    assert distance_value is not None
    assert no_overlap_policy_value is not None
    group_column = str(group_column_value).strip()
    no_overlap_policy = str(no_overlap_policy_value).strip()
    if no_overlap_policy != DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICY_ERROR:
        raise PhosPyInputError(
            "missing_data.policy='impute_group_aware' requires "
            "missing_data.no_overlap_policy='error'"
        )
    return (
        group_column,
        float(min_partial_value),
        float(min_reference_value),
        float(q_value),
        float(width_value),
        int(seed_value),
        int(k_value),
        str(distance_value).strip(),
        no_overlap_policy,
    )


def _validate_group_aware_merge(
    *,
    original_retained: pd.DataFrame,
    final_matrix: pd.DataFrame,
    knn_target_mask: pd.DataFrame,
    minprob_target_mask: pd.DataFrame,
    knn_imputed_mask: pd.DataFrame,
    minprob_imputed_mask: pd.DataFrame,
    final_imputed_mask: pd.DataFrame,
    knn_column_mean_fallback_mask: pd.DataFrame,
) -> None:
    aligned_frames = (
        final_matrix,
        knn_target_mask,
        minprob_target_mask,
        knn_imputed_mask,
        minprob_imputed_mask,
        final_imputed_mask,
        knn_column_mean_fallback_mask,
    )
    if any(
        not frame.index.equals(original_retained.index)
        or not frame.columns.equals(original_retained.columns)
        for frame in aligned_frames
    ):
        raise RuntimeError(
            "group-aware imputation produced misaligned matrices or masks"
        )
    if bool((knn_target_mask & minprob_target_mask).to_numpy().any()):
        raise RuntimeError("group-aware KNN and MinProb target masks overlap")
    if bool((knn_imputed_mask & ~knn_target_mask).to_numpy().any()):
        raise RuntimeError("group-aware KNN imputed cells outside its target mask")
    if bool((minprob_imputed_mask & ~minprob_target_mask).to_numpy().any()):
        raise RuntimeError("group-aware MinProb imputed cells outside its target mask")
    expected_final_mask = knn_imputed_mask | minprob_imputed_mask
    if not final_imputed_mask.equals(expected_final_mask):
        raise RuntimeError(
            "group-aware final imputation mask is not the union of mechanism masks"
        )
    if bool(knn_column_mean_fallback_mask.to_numpy().any()):
        raise RuntimeError("group-aware KNN unexpectedly used column-mean fallback")
    if bool(final_matrix.isna().to_numpy().any()):
        raise PhosPyInputError(
            "dataset preprocessing stage 'missing_data' could not complete "
            "missing_data.policy='impute_group_aware' because missing values "
            "remain after merging the independently produced KNN and MinProb results"
        )
    original_values = original_retained.to_numpy(
        dtype=float, copy=False, na_value=np.nan
    )
    final_values = final_matrix.to_numpy(dtype=float, copy=False, na_value=np.nan)
    originally_observed = ~np.isnan(original_values)
    if not np.array_equal(
        original_values[originally_observed],
        final_values[originally_observed],
    ):
        raise RuntimeError(
            "group-aware imputation changed originally observed retained values"
        )


def _mask_row_ids(mask: pd.DataFrame) -> tuple[str, ...]:
    if mask.empty:
        return ()
    return tuple(
        str(row_id)
        for row_id in mask.index[
            mask.any(axis=1).to_numpy(dtype=bool, copy=False)
        ].tolist()
    )


def _mask_column_ids(mask: pd.DataFrame) -> tuple[str, ...]:
    if mask.empty:
        return ()
    return tuple(
        str(column)
        for column in mask.columns[
            mask.any(axis=0).to_numpy(dtype=bool, copy=False)
        ].tolist()
    )


def _finalize_outcome(
    *,
    state: PreprocessingState,
    outcome: (
        RowMedianPolicyOutcome
        | KnnPolicyOutcome
        | MinProbPolicyOutcome
        | GroupAwarePolicyOutcome
    ),
    row_audit_records: Sequence[PreprocessingRowAuditRow],
    notes: str,
    diagnostics: Mapping[str, JsonValue],
) -> PreprocessingStageResult:
    observation_mask = _observation_mask_from_outcome(outcome)
    next_state = append_row_audit_records(state, row_audit_records)
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
    outcome: (
        RowMedianPolicyOutcome
        | KnnPolicyOutcome
        | MinProbPolicyOutcome
        | GroupAwarePolicyOutcome
    ),
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


def _required_imputation_input_scale_value(plan: PreprocessingPlan) -> str:
    if plan.missing_data_input_scale is None:
        raise PhosPyInputError(
            "dataset preprocessing stage 'missing_data' requires "
            "missing_data_input_scale for imputation policies"
        )
    return plan.missing_data_input_scale.value


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
        "missing_data_no_overlap_policy": plan.missing_data_no_overlap_policy,
        "missing_data_group_column": plan.missing_data_group_column,
        "missing_data_min_partial_observed_fraction": (
            plan.missing_data_min_partial_observed_fraction
        ),
        "missing_data_min_reference_observed_fraction": (
            plan.missing_data_min_reference_observed_fraction
        ),
        "missing_data_input_scale": (
            None
            if plan.missing_data_input_scale is None
            else plan.missing_data_input_scale.value
        ),
        "missing_data_input_scale_source": plan.missing_data_input_scale_source,
        "missing_data_imputation_operation_order": (
            plan.missing_data_imputation_operation_order
        ),
    }


def _resolve_quantitative_contract(
    plan: PreprocessingPlan,
) -> QuantitativeOperationContract:
    policy = plan.missing_data_policy
    if policy is MissingDataPolicy.FORBID:
        return preserve_quantitative_contract(
            required_evidence=frozenset({QuantitativeEvidenceRequirement.NONE}),
            negative_domain_policy=NegativeDomainPolicy.PRESERVES_INPUT_DOMAIN,
            reversibility=QuantitativeReversibilityKind.REVERSIBLE,
            information_loss=QuantitativeInformationLossKind.NONE,
        )
    if policy in {
        MissingDataPolicy.IMPUTE_MINPROB,
        MissingDataPolicy.IMPUTE_GROUP_AWARE,
    }:
        return QuantitativeOperationContract(
            accepted_input_scale_kinds=frozenset({IntensityScaleKind.LOG2}),
            accepted_quantitative_meanings=frozenset(
                {
                    QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
                    QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
                    QuantitativeMeaning.MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE,
                    QuantitativeMeaning.UNKNOWN,
                }
            ),
            output_scale_transition=preserve_scale_transition(
                frozenset({IntensityScaleKind.LOG2}),
                output_scale_label="log2",
            ),
            output_meaning_transition=preserve_meaning_transition(
                frozenset(
                    {
                        QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
                        QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
                        QuantitativeMeaning.MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE,
                        QuantitativeMeaning.UNKNOWN,
                    }
                )
            ),
            preserves_abundance=False,
            negative_domain_policy=NegativeDomainPolicy.REQUIRES_LOG2_DOMAIN,
            required_evidence=frozenset(
                {
                    QuantitativeEvidenceRequirement.MISSINGNESS_MASK,
                    QuantitativeEvidenceRequirement.RANDOM_SEED,
                }
            ),
            reversibility=QuantitativeReversibilityKind.IRREVERSIBLE,
            information_loss=QuantitativeInformationLossKind.IMPUTATION,
        )
    if policy in {MissingDataPolicy.IMPUTE_ROW_MEDIAN, MissingDataPolicy.IMPUTE_KNN}:
        input_scale_kind = imputation_input_scale_kind(plan.missing_data_input_scale)
        return QuantitativeOperationContract(
            accepted_input_scale_kinds=frozenset({input_scale_kind}),
            accepted_quantitative_meanings=ALL_QUANTITATIVE_MEANINGS,
            output_scale_transition=preserve_scale_transition(
                frozenset({input_scale_kind}),
                output_scale_label=input_scale_kind.value,
            ),
            output_meaning_transition=preserve_meaning_transition(),
            preserves_abundance=False,
            required_evidence=frozenset(
                {QuantitativeEvidenceRequirement.MISSINGNESS_MASK}
            ),
            negative_domain_policy=NegativeDomainPolicy.PRESERVES_INPUT_DOMAIN,
            reversibility=QuantitativeReversibilityKind.IRREVERSIBLE,
            information_loss=QuantitativeInformationLossKind.IMPUTATION,
        )
    raise PhosPyInputError(
        "dataset build request preprocessing_config contains an unsupported "
        "missing_data.policy"
    )


def _resolve_determinism_kind(plan: PreprocessingPlan) -> DeterminismKind:
    if plan.missing_data_policy in {
        MissingDataPolicy.IMPUTE_MINPROB,
        MissingDataPolicy.IMPUTE_GROUP_AWARE,
    }:
        return DeterminismKind.SEEDED_STOCHASTIC
    return DeterminismKind.DETERMINISTIC


def _resolve_consumed_input_tables(
    plan: PreprocessingPlan,
) -> tuple[PreprocessingStateTableKey, ...]:
    tables = (
        PreprocessingStateTableKey.DATASET_PHOSPHO,
        PreprocessingStateTableKey.DATASET_SITE_METADATA,
    )
    if plan.missing_data_policy is MissingDataPolicy.IMPUTE_GROUP_AWARE:
        return (*tables, PreprocessingStateTableKey.DATASET_SAMPLE_METADATA)
    return tables


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
    quantitative_contract=_resolve_quantitative_contract,
    stage_factory=_build_missing_data_stage,
    backend="pandas",
    determinism_kind=_resolve_determinism_kind,
    consumed_input_tables_resolver=_resolve_consumed_input_tables,
    diagnostics_metadata={
        "diagnostics_schema_version": MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V2,
        "supported_diagnostics_schema_versions": (
            MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1,
            MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V2,
        ),
        "known_diagnostics_fields": tuple(
            sorted(V2_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS)
        ),
        "known_diagnostics_fields_by_version": {
            str(MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1): tuple(
                sorted(V1_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS)
            ),
            str(MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V2): tuple(
                sorted(V2_KNOWN_MISSING_DATA_DIAGNOSTICS_FIELDS)
            ),
        },
    },
)


__all__ = ["MISSING_DATA_STAGE_CONTRACT", "MissingDataStage"]
