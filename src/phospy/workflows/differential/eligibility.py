"""Differential feature eligibility resolution components."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, cast

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from phospy.contracts.configs.differential import (
    IMPUTED_VALUE_POLICY_REJECT,
    IMPUTED_VALUE_POLICY_WITHHOLD_IMPUTED_FEATURES,
    PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
)
from phospy.errors.validation import DatasetValidationError
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.design.models import Contrast, ExperimentalDesign
from phospy.science.differential.linear_model import (
    DifferentialDesignDecomposition,
    DifferentialDesignDecompositionError,
)
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_INVALID_NUMERIC_VALUES,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_OTHER,
)
from phospy.science.differential.models import (
    DifferentialAnalysisRequest as DifferentialComputationRequest,
)
from phospy.workflows._pandas_typing import (
    dataframe_column,
    dataframe_copy,
    dataframe_loc,
    index_as_strings,
    index_snapshot,
    series_copy,
)
from phospy.workflows.differential.imputation_inference import (
    imputation_inference_columns,
)
from phospy.workflows.differential.models import (
    DifferentialFeatureEligibilityInputs,
    DifferentialImputationPolicyInputs,
)

if TYPE_CHECKING:
    from phospy.workflows.differential.models import (
        InterpretedDifferentialAnalysisRequest,
    )


@dataclass(frozen=True, slots=True)
class DifferentialPreFitEligibilityResolution:
    """Pre-fit feature eligibility resolved by the interpreter."""

    imputation_policy_inputs: DifferentialImputationPolicyInputs | None
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs


@dataclass(frozen=True, slots=True)
class DifferentialPostFitEligibilityResolution:
    """Post-fit numerical eligibility resolved before public result expansion."""

    model_fit_feature_ids: tuple[str, ...]
    failed_model_fit_feature_ids: tuple[str, ...]
    valid_model_fit_feature_ids: tuple[str, ...]
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs | None


@dataclass(frozen=True, slots=True)
class DifferentialExecutionEligibilityResolution:
    """Execution matrix and row attrition facts after eligibility filtering."""

    computation_request: DifferentialComputationRequest
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs | None
    input_feature_ids: tuple[str, ...]
    model_fit_feature_ids: tuple[str, ...]
    failed_model_fit_feature_ids: tuple[str, ...]
    multiple_testing_feature_ids: tuple[str, ...]


class DifferentialPreFitEligibilityResolver:
    """Resolve finite-value, constant-row, and imputation-based eligibility."""

    def run(
        self,
        *,
        dataset_view: DatasetInternalView,
        matrix: pd.DataFrame,
        analysis_sample_ids: tuple[str, ...],
        design: ExperimentalDesign,
        contrasts: tuple[Contrast, ...],
        policy: str,
        max_fraction: float,
        minimum_condition_replicates: int,
    ) -> DifferentialPreFitEligibilityResolution:
        imputation_policy_inputs = _build_imputation_policy_inputs(
            dataset_view=dataset_view,
            matrix_index=matrix.index,
            analysis_sample_ids=analysis_sample_ids,
            design=design,
            contrasts=contrasts,
            policy=policy,
            max_fraction=max_fraction,
            minimum_condition_replicates=minimum_condition_replicates,
        )
        feature_eligibility_inputs = _build_feature_eligibility_inputs(
            matrix=matrix,
            imputation_policy_inputs=imputation_policy_inputs,
        )
        return DifferentialPreFitEligibilityResolution(
            imputation_policy_inputs=imputation_policy_inputs,
            feature_eligibility_inputs=feature_eligibility_inputs,
        )


class DifferentialPostFitEligibilityResolver:
    """Classify rows that fail model-fit numerical eligibility checks."""

    def run(
        self,
        *,
        computation_request: DifferentialComputationRequest,
        feature_eligibility_inputs: DifferentialFeatureEligibilityInputs | None,
    ) -> DifferentialPostFitEligibilityResolution:
        model_fit_feature_ids = feature_ids(computation_request.matrix.index)
        failed_model_fit_feature_ids = _failed_model_fit_feature_ids(
            computation_request
        )
        failed_model_fit_feature_id_set = set(failed_model_fit_feature_ids)
        valid_model_fit_feature_ids = tuple(
            feature_id
            for feature_id in model_fit_feature_ids
            if feature_id not in failed_model_fit_feature_id_set
        )
        resolved_feature_eligibility_inputs = feature_eligibility_inputs
        if (
            failed_model_fit_feature_ids
            and resolved_feature_eligibility_inputs is not None
        ):
            resolved_feature_eligibility_inputs = _with_failed_model_fit_status(
                feature_eligibility_inputs=resolved_feature_eligibility_inputs,
                failed_feature_ids=failed_model_fit_feature_ids,
            )
        return DifferentialPostFitEligibilityResolution(
            model_fit_feature_ids=model_fit_feature_ids,
            failed_model_fit_feature_ids=failed_model_fit_feature_ids,
            valid_model_fit_feature_ids=valid_model_fit_feature_ids,
            feature_eligibility_inputs=resolved_feature_eligibility_inputs,
        )


class DifferentialComputationEligibilityResolver:
    """Resolve execution-time feature eligibility before model fitting."""

    def __init__(
        self,
        *,
        post_fit_eligibility_resolver: DifferentialPostFitEligibilityResolver
        | None = None,
    ) -> None:
        self._post_fit_eligibility_resolver = (
            post_fit_eligibility_resolver or DifferentialPostFitEligibilityResolver()
        )

    def run(
        self,
        request: InterpretedDifferentialAnalysisRequest,
    ) -> DifferentialExecutionEligibilityResolution:
        computation_request = request.computation_request
        imputation_policy_inputs = request.imputation_policy_inputs
        feature_eligibility_inputs = request.feature_eligibility_inputs
        if feature_eligibility_inputs is not None:
            computation_request = _filter_computation_request_for_feature_eligibility(
                computation_request=request.computation_request,
                feature_eligibility_inputs=feature_eligibility_inputs,
            )
        elif imputation_policy_inputs is not None:
            computation_request = _filter_computation_request_for_imputation_policy(
                computation_request=request.computation_request,
                imputation_policy_inputs=imputation_policy_inputs,
            )
        if (
            request.execution_config.paired_design_policy
            == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
        ):
            model_fit_feature_ids = feature_ids(computation_request.matrix.index)
            return DifferentialExecutionEligibilityResolution(
                computation_request=computation_request,
                feature_eligibility_inputs=feature_eligibility_inputs,
                input_feature_ids=feature_ids(request.result_identity_metadata.index),
                model_fit_feature_ids=model_fit_feature_ids,
                failed_model_fit_feature_ids=(),
                multiple_testing_feature_ids=model_fit_feature_ids,
            )
        post_fit_eligibility = self._post_fit_eligibility_resolver.run(
            computation_request=computation_request,
            feature_eligibility_inputs=feature_eligibility_inputs,
        )
        model_fit_feature_ids = post_fit_eligibility.model_fit_feature_ids
        failed_model_fit_feature_ids = post_fit_eligibility.failed_model_fit_feature_ids
        if failed_model_fit_feature_ids:
            computation_request = _filter_computation_request_for_feature_ids(
                computation_request=computation_request,
                feature_ids=post_fit_eligibility.valid_model_fit_feature_ids,
                seam="differential.executor.model_fit_valid_features",
                next_action=(
                    "provide at least one feature whose fitted residual variance is "
                    "finite and positive under the resolved differential design"
                ),
                details={
                    "failed_model_fit_feature_ids": list(
                        failed_model_fit_feature_ids[:5]
                    ),
                    "failed_model_fit_count": int(len(failed_model_fit_feature_ids)),
                },
            )
            feature_eligibility_inputs = post_fit_eligibility.feature_eligibility_inputs
        return DifferentialExecutionEligibilityResolution(
            computation_request=computation_request,
            feature_eligibility_inputs=feature_eligibility_inputs,
            input_feature_ids=feature_ids(request.result_identity_metadata.index),
            model_fit_feature_ids=model_fit_feature_ids,
            failed_model_fit_feature_ids=failed_model_fit_feature_ids,
            multiple_testing_feature_ids=feature_ids(computation_request.matrix.index),
        )


def filter_matrix_for_feature_ids(
    *,
    matrix: pd.DataFrame,
    feature_ids: tuple[str, ...],
) -> pd.DataFrame:
    """Return a workflow-local matrix restricted to resolved feature IDs."""

    current_feature_ids = tuple(str(feature_id) for feature_id in matrix.index)
    if current_feature_ids == tuple(feature_ids):
        return dataframe_copy(matrix, deep=False)
    row_positions_by_feature_id = {
        str(feature_id): position for position, feature_id in enumerate(matrix.index)
    }
    missing_feature_ids = [
        feature_id
        for feature_id in feature_ids
        if feature_id not in row_positions_by_feature_id
    ]
    if missing_feature_ids:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.feature_eligibility_alignment",
            next_action=(
                "ensure the interpreted differential matrix contains every "
                "feature selected by pre-fit eligibility"
            ),
            details={"missing_feature_ids": missing_feature_ids[:5]},
            message_prefix="differential workflow boundary validation failed",
        )
    row_positions = [
        row_positions_by_feature_id[feature_id] for feature_id in feature_ids
    ]
    return pd.DataFrame(
        matrix.to_numpy(dtype=float)[row_positions, :],
        index=pd.Index(feature_ids, name=matrix.index.name),
        columns=index_snapshot(matrix.columns),
    )


def feature_ids(index: pd.Index) -> tuple[str, ...]:
    """Return stable string feature identifiers from a pandas index."""

    return tuple(str(feature_id) for feature_id in index.tolist())


def differential_status_counts(result_status: pd.Series) -> dict[str, int]:
    """Count differential result statuses without changing encounter order."""

    counts = result_status.astype(str).value_counts(sort=False)
    return {str(status): int(count) for status, count in counts.items()}


def _filter_computation_request_for_imputation_policy(
    *,
    computation_request: DifferentialComputationRequest,
    imputation_policy_inputs: DifferentialImputationPolicyInputs,
) -> DifferentialComputationRequest:
    return _filter_computation_request_for_feature_ids(
        computation_request=computation_request,
        feature_ids=imputation_policy_inputs.testable_feature_ids,
        seam="differential.executor.imputation_policy_testable_features",
        next_action=(
            "raise differential.imputed_value_max_fraction, require more "
            "observed values per condition, or use a non-imputed dataset"
        ),
        details={"policy": imputation_policy_inputs.policy},
    )


def _filter_computation_request_for_feature_eligibility(
    *,
    computation_request: DifferentialComputationRequest,
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs,
) -> DifferentialComputationRequest:
    return _filter_computation_request_for_feature_ids(
        computation_request=computation_request,
        feature_ids=feature_eligibility_inputs.testable_feature_ids,
        seam="differential.executor.feature_eligibility_testable_features",
        next_action=(
            "ensure differential feature eligibility is resolved before "
            "statistical execution and at least one feature remains testable"
        ),
        details={
            "status_counts": differential_status_counts(
                feature_eligibility_inputs.result_status
            )
        },
    )


def _filter_computation_request_for_feature_ids(
    *,
    computation_request: DifferentialComputationRequest,
    feature_ids: tuple[str, ...],
    seam: str,
    next_action: str,
    details: dict[str, object],
) -> DifferentialComputationRequest:
    testable_feature_ids = list(feature_ids)
    if not testable_feature_ids:
        raise WorkflowBoundaryError(
            seam=seam,
            next_action=next_action,
            details=details,
            message_prefix="differential workflow boundary validation failed",
        )
    current_feature_ids = tuple(
        str(feature_id) for feature_id in computation_request.matrix.index.tolist()
    )
    if current_feature_ids == tuple(testable_feature_ids):
        return computation_request
    row_positions_by_feature_id = {
        str(feature_id): position
        for position, feature_id in enumerate(computation_request.matrix.index)
    }
    missing_feature_ids = [
        feature_id
        for feature_id in testable_feature_ids
        if feature_id not in row_positions_by_feature_id
    ]
    if missing_feature_ids:
        raise WorkflowBoundaryError(
            seam=seam,
            next_action=(
                "ensure the interpreted computation matrix contains every "
                "testable feature selected by differential feature eligibility"
            ),
            details={**details, "missing_feature_ids": missing_feature_ids[:5]},
            message_prefix="differential workflow boundary validation failed",
        )
    row_positions = [
        row_positions_by_feature_id[feature_id] for feature_id in testable_feature_ids
    ]
    filtered_matrix = pd.DataFrame(
        computation_request.matrix.to_numpy(dtype=float)[row_positions, :],
        index=pd.Index(
            testable_feature_ids,
            name=computation_request.matrix.index.name,
        ),
        columns=index_snapshot(computation_request.matrix.columns),
    )
    return DifferentialComputationRequest(
        matrix=filtered_matrix,
        design=computation_request.design,
        contrasts=computation_request.contrasts,
        design_decomposition=computation_request.design_decomposition,
        empirical_bayes=computation_request.empirical_bayes,
        multiple_testing_method=computation_request.multiple_testing_method,
    )


def _build_feature_eligibility_inputs(
    *,
    matrix: pd.DataFrame,
    imputation_policy_inputs: DifferentialImputationPolicyInputs | None,
) -> DifferentialFeatureEligibilityInputs:
    feature_metadata = _base_feature_eligibility_metadata(matrix)
    statuses = feature_metadata[DIFFERENTIAL_RESULT_STATUS_COLUMN].astype(str)
    reasons = feature_metadata[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN].astype(str)

    if imputation_policy_inputs is not None:
        if (
            not imputation_policy_inputs.feature_metadata.index.equals(matrix.index)
            or not imputation_policy_inputs.result_status.index.equals(matrix.index)
            or not imputation_policy_inputs.result_status_reason.index.equals(
                matrix.index
            )
        ):
            raise WorkflowBoundaryError(
                seam="differential.interpreter.feature_eligibility_alignment",
                next_action=(
                    "ensure feature eligibility and imputation policy metadata use "
                    "the same feature index"
                ),
                message_prefix="differential workflow boundary validation failed",
            )
        imputation_metadata = imputation_policy_inputs.feature_metadata
        for column_name in (
            "imputed_cell_count",
            "observed_cell_count",
            "imputed_fraction",
        ):
            feature_metadata[column_name] = dataframe_column(
                imputation_metadata,
                column_name,
            ).to_numpy()
        feature_metadata["imputation_policy"] = imputation_policy_inputs.policy
        feature_metadata["imputation_fraction_threshold"] = (
            imputation_policy_inputs.max_fraction
        )

        imputation_statuses = imputation_policy_inputs.result_status.astype(str)
        imputation_reasons = imputation_policy_inputs.result_status_reason.astype(str)
        merged_statuses = statuses.copy(deep=True)
        merged_reasons = reasons.copy(deep=True)
        base_tested = statuses == DIFFERENTIAL_RESULT_STATUS_TESTED
        imputation_withheld = imputation_statuses != DIFFERENTIAL_RESULT_STATUS_TESTED
        imputation_applies = base_tested & imputation_withheld
        merged_statuses.loc[imputation_applies] = imputation_statuses.loc[
            imputation_applies
        ]
        merged_reasons.loc[imputation_applies] = imputation_reasons.loc[
            imputation_applies
        ]
        statuses = merged_statuses
        reasons = merged_reasons
        for column_name, values in imputation_inference_columns(
            feature_metadata=feature_metadata,
            result_status=statuses,
        ).items():
            feature_metadata[column_name] = values

    feature_metadata[DIFFERENTIAL_RESULT_STATUS_COLUMN] = statuses.to_numpy(dtype=str)
    feature_metadata[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] = reasons.to_numpy(
        dtype=str
    )
    result_status = pd.Series(
        statuses.to_numpy(dtype=str),
        index=index_snapshot(matrix.index),
        name=DIFFERENTIAL_RESULT_STATUS_COLUMN,
    )
    testable_feature_ids = tuple(
        str(feature_id)
        for feature_id, status in zip(
            result_status.index,
            result_status.to_numpy(dtype=str),
            strict=True,
        )
        if status == DIFFERENTIAL_RESULT_STATUS_TESTED
    )
    attach_to_result_tables = bool(
        imputation_policy_inputs is not None
        or (result_status != DIFFERENTIAL_RESULT_STATUS_TESTED).any()
    )
    return DifferentialFeatureEligibilityInputs(
        feature_metadata=feature_metadata,
        result_status=result_status,
        testable_feature_ids=testable_feature_ids,
        attach_to_result_tables=attach_to_result_tables,
    )


def _base_feature_eligibility_metadata(matrix: pd.DataFrame) -> pd.DataFrame:
    values = np.asarray(matrix.to_numpy(dtype=float), dtype=np.float64)
    finite_mask = np.isfinite(values)
    analysed_value_count = int(values.shape[1])
    observed_value_counts = finite_mask.sum(axis=1).astype(np.int64)
    invalid_numeric_counts = (~finite_mask).sum(axis=1).astype(np.int64)
    unique_observed_counts: list[int] = []
    statuses: list[str] = []
    reasons: list[str] = []

    for row_position in range(int(values.shape[0])):
        finite_values = values[row_position, finite_mask[row_position, :]]
        unique_count = int(np.unique(finite_values).size)
        unique_observed_counts.append(unique_count)
        invalid_count = int(invalid_numeric_counts[row_position])
        if invalid_count > 0:
            statuses.append(DIFFERENTIAL_RESULT_STATUS_WITHHELD_INVALID_NUMERIC_VALUES)
            reasons.append(
                "Feature contains non-finite numeric values in the differential "
                "design samples."
            )
            continue
        if unique_count <= 1:
            statuses.append(DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT)
            reasons.append(
                "Feature is all-constant across the differential design samples."
            )
            continue
        statuses.append(DIFFERENTIAL_RESULT_STATUS_TESTED)
        reasons.append(
            "Feature has finite, non-constant values for the differential design "
            "samples."
        )

    return pd.DataFrame(
        {
            "site_key": index_as_strings(matrix.index),
            DIFFERENTIAL_RESULT_STATUS_COLUMN: statuses,
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN: reasons,
            "analysed_value_count": np.full(
                int(matrix.shape[0]),
                analysed_value_count,
                dtype=np.int64,
            ),
            "observed_value_count": observed_value_counts,
            "invalid_numeric_value_count": invalid_numeric_counts,
            "unique_observed_value_count": np.asarray(
                unique_observed_counts,
                dtype=np.int64,
            ),
        },
        index=index_snapshot(matrix.index),
    )


def _build_imputation_policy_inputs(
    *,
    dataset_view: DatasetInternalView,
    matrix_index: pd.Index,
    analysis_sample_ids: tuple[str, ...],
    design: ExperimentalDesign,
    contrasts: tuple[Contrast, ...],
    policy: str,
    max_fraction: float,
    minimum_condition_replicates: int,
) -> DifferentialImputationPolicyInputs | None:
    if policy == IMPUTED_VALUE_POLICY_REJECT:
        return None
    if policy != IMPUTED_VALUE_POLICY_WITHHOLD_IMPUTED_FEATURES:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.imputation_policy",
            next_action="use a supported differential imputation policy",
            details={"policy": policy},
            message_prefix="differential workflow boundary validation failed",
        )
    feature_metadata = _imputation_summary_for_differential(
        dataset_view=dataset_view,
        matrix_index=matrix_index,
        sample_ids=analysis_sample_ids,
    )
    imputed_fraction = feature_metadata["imputed_fraction"].to_numpy(dtype=float)
    condition_sample_ids = _condition_sample_ids_for_analysis(
        design=design,
        analysis_sample_ids=analysis_sample_ids,
    )
    condition_observed_counts = _condition_observed_counts_for_analysis(
        dataset_view=dataset_view,
        matrix_index=matrix_index,
        condition_sample_ids=condition_sample_ids,
    )
    statuses: list[str] = []
    reasons: list[str] = []
    for row_position in range(int(matrix_index.size)):
        if _has_insufficient_observed_samples_for_contrasts(
            row_position=row_position,
            condition_observed_counts=condition_observed_counts,
            contrasts=contrasts,
            minimum_condition_replicates=minimum_condition_replicates,
        ):
            statuses.append(DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED)
            reasons.append(
                "Feature has fewer observed values than required in at least one "
                "contrasted condition."
            )
            continue
        if float(imputed_fraction[row_position]) > float(max_fraction):
            statuses.append(DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION)
            reasons.append(
                "Feature exceeds the configured imputed-value fraction threshold."
            )
            continue
        statuses.append(DIFFERENTIAL_RESULT_STATUS_TESTED)
        reasons.append("Feature passed imputation policy checks.")

    result_status = pd.Series(
        statuses,
        index=index_snapshot(matrix_index),
        name="result_status",
    )
    result_status_reason = pd.Series(
        reasons,
        index=index_snapshot(matrix_index),
        name=DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    )
    testable_feature_ids = tuple(
        str(feature_id)
        for feature_id, status in zip(
            result_status.index,
            result_status.to_numpy(dtype=str),
            strict=True,
        )
        if status == DIFFERENTIAL_RESULT_STATUS_TESTED
    )
    return DifferentialImputationPolicyInputs(
        feature_metadata=feature_metadata,
        result_status=result_status,
        result_status_reason=result_status_reason,
        testable_feature_ids=testable_feature_ids,
        policy=policy,
        max_fraction=max_fraction,
    )


def _condition_sample_ids_for_analysis(
    *,
    design: ExperimentalDesign,
    analysis_sample_ids: tuple[str, ...],
) -> dict[str, tuple[str, ...]]:
    analysis_sample_set = set(analysis_sample_ids)
    grouped: dict[str, list[str]] = {}
    for record in design.samples:
        sample_id = str(record.sample_id)
        if sample_id not in analysis_sample_set:
            continue
        grouped.setdefault(record.condition, []).append(sample_id)
    return {condition: tuple(sample_ids) for condition, sample_ids in grouped.items()}


def _condition_observed_counts_for_analysis(
    *,
    dataset_view: DatasetInternalView,
    matrix_index: pd.Index,
    condition_sample_ids: dict[str, tuple[str, ...]],
) -> dict[str, pd.Series]:
    observed_counts: dict[str, pd.Series] = {}
    for condition, sample_ids in condition_sample_ids.items():
        if not sample_ids:
            continue
        summary = _imputation_summary_for_differential(
            dataset_view=dataset_view,
            matrix_index=matrix_index,
            sample_ids=sample_ids,
        )
        observed_counts[condition] = series_copy(
            dataframe_column(summary, "observed_cell_count"),
            deep=True,
        )
    return observed_counts


def _imputation_summary_for_differential(
    *,
    dataset_view: DatasetInternalView,
    matrix_index: pd.Index,
    sample_ids: tuple[str, ...],
) -> pd.DataFrame:
    try:
        summary = dataset_view.imputation_observation_summary(
            feature_ids=tuple(index_as_strings(matrix_index)),
            sample_ids=sample_ids,
        )
    except DatasetValidationError as exc:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.imputation_metadata_alignment",
            next_action=(
                "ensure imputation observation metadata is aligned to the "
                "differential feature and sample labels"
            ),
            details={
                "feature_count": int(matrix_index.size),
                "sample_count": len(sample_ids),
                "error": str(exc),
            },
            message_prefix="differential workflow boundary validation failed",
        ) from exc
    if summary is None:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.imputation_metadata",
            next_action=(
                "build the analysis-ready dataset through a supported imputation "
                "preprocessing path that preserves the observed-cell mask"
            ),
            message_prefix="differential workflow boundary validation failed",
        )
    required_columns = (
        "feature_id",
        "observed_cell_count",
        "imputed_cell_count",
        "total_analysed_cell_count",
        "imputed_fraction",
    )
    missing_columns = [
        column for column in required_columns if column not in summary.columns
    ]
    if missing_columns:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.imputation_metadata_contract",
            next_action=(
                "ensure the dataset imputation summary includes feature-level "
                "observed and imputed cell counts"
            ),
            details={"missing_columns": missing_columns},
            message_prefix="differential workflow boundary validation failed",
        )
    if not summary.index.equals(matrix_index):
        raise WorkflowBoundaryError(
            seam="differential.interpreter.imputation_metadata_alignment",
            next_action=(
                "ensure imputation observation summary order exactly matches the "
                "differential execution matrix"
            ),
            message_prefix="differential workflow boundary validation failed",
        )
    return dataframe_copy(summary, deep=True)


def _has_insufficient_observed_samples_for_contrasts(
    *,
    row_position: int,
    condition_observed_counts: dict[str, pd.Series],
    contrasts: tuple[Contrast, ...],
    minimum_condition_replicates: int,
) -> bool:
    for contrast in contrasts:
        for condition in (
            contrast.numerator_condition,
            contrast.denominator_condition,
        ):
            observed_counts = condition_observed_counts.get(condition)
            if observed_counts is None:
                return True
            observed_count = int(observed_counts.iloc[row_position])
            if observed_count < int(minimum_condition_replicates):
                return True
    return False


def _failed_model_fit_feature_ids(
    computation_request: DifferentialComputationRequest,
) -> tuple[str, ...]:
    matrix = computation_request.matrix
    design_frame = computation_request.design.frame
    matrix_aligned = dataframe_loc(matrix, columns=list(design_frame.index))
    design_decomposition = cast(
        DifferentialDesignDecomposition,
        computation_request.design_decomposition,
    )
    try:
        design_decomposition.assert_matches_design(
            design_frame.to_numpy(dtype=float),
            field_name="differential.design",
        )
    except DifferentialDesignDecompositionError:
        return ()

    matrix_values: NDArray[np.float64] = np.asarray(
        matrix_aligned.to_numpy(dtype=float),
        dtype=np.float64,
    )
    response: NDArray[np.float64] = np.transpose(matrix_values)
    try:
        residual_variance = design_decomposition.fit(response).residual_variance
    except DifferentialDesignDecompositionError:
        return ()
    failed_mask = ~np.isfinite(residual_variance) | (residual_variance <= 0.0)
    return tuple(
        str(feature_id)
        for feature_id, failed in zip(matrix_aligned.index, failed_mask, strict=True)
        if bool(failed)
    )


def _with_failed_model_fit_status(
    *,
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs,
    failed_feature_ids: tuple[str, ...],
) -> DifferentialFeatureEligibilityInputs:
    failed_feature_id_set = set(failed_feature_ids)
    feature_metadata = dataframe_copy(
        feature_eligibility_inputs.feature_metadata,
        deep=True,
    )
    status_values: NDArray[np.object_] = np.array(
        dataframe_column(feature_metadata, DIFFERENTIAL_RESULT_STATUS_COLUMN).to_numpy(
            dtype=object
        ),
        dtype=object,
        copy=True,
    )
    reason_values: NDArray[np.object_] = np.array(
        dataframe_column(
            feature_metadata,
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
        ).to_numpy(dtype=object),
        dtype=object,
        copy=True,
    )
    fit_failure_reason = (
        "Feature model fit failed before multiple-testing correction; residual "
        "variance was zero or non-finite."
    )
    for row_position, label in enumerate(feature_metadata.index):
        if str(label) not in failed_feature_id_set:
            continue
        status_values[int(row_position)] = DIFFERENTIAL_RESULT_STATUS_WITHHELD_OTHER
        reason_values[int(row_position)] = fit_failure_reason
    feature_metadata[DIFFERENTIAL_RESULT_STATUS_COLUMN] = status_values
    feature_metadata[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] = reason_values
    result_status = pd.Series(
        status_values.astype(str),
        index=index_snapshot(feature_eligibility_inputs.result_status.index),
        name=DIFFERENTIAL_RESULT_STATUS_COLUMN,
    )
    return DifferentialFeatureEligibilityInputs(
        feature_metadata=feature_metadata,
        result_status=result_status,
        testable_feature_ids=tuple(
            feature_id
            for feature_id in feature_eligibility_inputs.testable_feature_ids
            if feature_id not in failed_feature_id_set
        ),
        attach_to_result_tables=True,
    )


__all__ = [
    "DifferentialComputationEligibilityResolver",
    "DifferentialExecutionEligibilityResolution",
    "DifferentialPostFitEligibilityResolution",
    "DifferentialPostFitEligibilityResolver",
    "DifferentialPreFitEligibilityResolution",
    "DifferentialPreFitEligibilityResolver",
    "differential_status_counts",
    "feature_ids",
    "filter_matrix_for_feature_ids",
]
