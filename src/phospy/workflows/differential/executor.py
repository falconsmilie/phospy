"""Internal executor for differential workflow requests."""

from __future__ import annotations

from typing import cast

import numpy as np
import pandas as pd

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.datasets.preprocessing.batch_correction import (
    BatchCorrectionReport,
)
from phospy.science.differential.executor import (
    DifferentialAnalysisExecutor as DifferentialComputationExecutor,
)
from phospy.science.differential.models import (
    DifferentialAnalysisRequest as DifferentialComputationRequest,
)
from phospy.science.differential.models import (
    DifferentialAnalysisResult,
    DifferentialComputationResult,
    DifferentialContrastDefinition,
    DifferentialModelDiagnostics,
    EmpiricalBayesPriorDiagnostics,
    MeanVarianceTrendDiagnostics,
)
from phospy.workflows.differential.models import (
    DifferentialExecutionDesignInputs,
    DifferentialFeatureEligibilityInputs,
    DifferentialImputationPolicyInputs,
    InterpretedDifferentialAnalysisRequest,
)


class DifferentialAnalysisExecutor:
    """Run differential computation for interpreted execution inputs."""

    def __init__(
        self,
        *,
        computation_executor: DifferentialComputationExecutor | None = None,
    ) -> None:
        self._computation_executor = (
            computation_executor or DifferentialComputationExecutor()
        )

    def run(
        self, request: InterpretedDifferentialAnalysisRequest
    ) -> DifferentialAnalysisResult:
        if not isinstance(
            cast(object, request), InterpretedDifferentialAnalysisRequest
        ):
            raise WorkflowBoundaryError(
                seam="differential.executor.interpreted_request_type",
                next_action=(
                    "pass interpreter output into DifferentialAnalysisExecutor.run"
                ),
                message_prefix="differential workflow boundary validation failed",
            )
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
        result = self._computation_executor.run(computation_request)
        residual_variance = result.residual_variance
        posterior_residual_variance = result.posterior_residual_variance
        prior_residual_variance = result.prior_residual_variance
        prior_degrees_of_freedom_series_value = (
            result.prior_degrees_of_freedom_series_value
        )
        prior_diagnostics = result.prior_diagnostics
        mean_variance_trend_diagnostics = result.mean_variance_trend_diagnostics
        contrast_source_tables = dict(result._contrast_tables.items())  # pyright: ignore[reportPrivateUsage] - owned contrast tables are forwarded without copying
        full_index = request.result_identity_metadata.index
        if not result.residual_variance.index.equals(full_index):
            residual_variance = _expand_series_to_full_index(
                result.residual_variance,
                full_index=full_index,
            )
            posterior_residual_variance = _expand_series_to_full_index(
                result.posterior_residual_variance,
                full_index=full_index,
            )
            prior_residual_variance = _expand_series_to_full_index(
                result.prior_residual_variance,
                full_index=full_index,
            )
            prior_degrees_of_freedom_series_value = _expand_series_to_full_index(
                result.prior_degrees_of_freedom_series_value,
                full_index=full_index,
            )
            prior_diagnostics = _expand_prior_diagnostics_to_full_index(
                result.prior_diagnostics,
                full_index=full_index,
            )
            mean_variance_trend_diagnostics = _expand_trend_diagnostics_to_full_index(
                result.mean_variance_trend_diagnostics,
                full_index=full_index,
            )
            contrast_source_tables = {
                contrast_name: _expand_stat_table_to_full_index(
                    table,
                    full_index=full_index,
                )
                for contrast_name, table in contrast_source_tables.items()
            }
        contrast_tables = {
            contrast_name: _attach_result_identity_metadata(
                table=table,
                identity_metadata=request.result_identity_metadata,
                contrast_name=contrast_name,
                imputation_policy_inputs=imputation_policy_inputs,
                feature_eligibility_inputs=feature_eligibility_inputs,
            )
            for contrast_name, table in contrast_source_tables.items()
        }
        diagnostics = _build_model_diagnostics(request=request, result=result)
        return DifferentialAnalysisResult._from_owned(  # pyright: ignore[reportPrivateUsage] - trusted internal ownership-preserving constructor
            residual_variance=residual_variance,
            posterior_residual_variance=posterior_residual_variance,
            prior_residual_variance=prior_residual_variance,
            prior_degrees_of_freedom_series_value=(
                prior_degrees_of_freedom_series_value
            ),
            prior_variance=result.prior_variance,
            prior_degrees_of_freedom=result.prior_degrees_of_freedom,
            residual_degrees_of_freedom=result.residual_degrees_of_freedom,
            empirical_bayes_method=result.empirical_bayes_method,
            empirical_bayes_robust=result.empirical_bayes_robust,
            empirical_bayes_trend=result.empirical_bayes_trend,
            prior_diagnostics=prior_diagnostics,
            mean_variance_trend_diagnostics=mean_variance_trend_diagnostics,
            diagnostics=diagnostics,
            policy_provenance=request.policy_provenance,
            contrast_tables=contrast_tables,
            workflow_provenance=request.workflow_provenance,
            input_dataset_preprocessing_report=request.dataset_preprocessing_report,
            feature_eligibility=(
                None
                if feature_eligibility_inputs is None
                else feature_eligibility_inputs.feature_metadata
            ),
        )


def _attach_result_identity_metadata(
    *,
    table: pd.DataFrame,
    identity_metadata: pd.DataFrame,
    contrast_name: str,
    imputation_policy_inputs: DifferentialImputationPolicyInputs | None = None,
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs | None = None,
) -> pd.DataFrame:
    if not table.index.equals(identity_metadata.index):
        raise WorkflowBoundaryError(
            seam="differential.executor.result_identity_alignment",
            next_action=(
                "ensure interpreted result_identity_metadata index exactly matches "
                "differential contrast table index"
            ),
            details={"contrast_name": contrast_name},
            message_prefix="differential workflow boundary validation failed",
        )
    enriched = pd.DataFrame(identity_metadata, copy=True)
    if (
        feature_eligibility_inputs is not None
        and feature_eligibility_inputs.attach_to_result_tables
    ):
        _attach_feature_eligibility_metadata(
            enriched=enriched,
            feature_eligibility_inputs=feature_eligibility_inputs,
            contrast_name=contrast_name,
        )
    elif imputation_policy_inputs is not None:
        _attach_imputation_policy_metadata(
            enriched=enriched,
            imputation_policy_inputs=imputation_policy_inputs,
            contrast_name=contrast_name,
        )
    for column_name in ("logFC", "t", "P.Value", "adj.P.Val"):
        contrast_column = table[column_name]
        enriched[column_name] = contrast_column.to_numpy(dtype=float)
    return enriched


def _build_model_diagnostics(
    *,
    request: InterpretedDifferentialAnalysisRequest,
    result: DifferentialComputationResult,
) -> DifferentialModelDiagnostics:
    design_frame = request.computation_request.design.frame
    policy = request.policy_provenance
    execution_design = request.execution_design
    design_columns = tuple(str(label) for label in design_frame.columns)
    contrast_definitions = (
        policy.contrasts
        if policy is not None
        else _contrast_definitions_from_matrix(request)
    )
    batch_or_covariate_terms = _batch_or_covariate_terms(execution_design)
    unsupported_assumptions = _unsupported_assumptions(
        request=request,
        batch_or_covariate_terms=batch_or_covariate_terms,
    )
    warnings = _diagnostic_warnings(
        request=request,
        batch_or_covariate_terms=batch_or_covariate_terms,
    )
    return DifferentialModelDiagnostics(
        model_type="moderated_ols_fixed_effect",
        design_columns=design_columns,
        contrast_definitions=contrast_definitions,
        rank=int(request.design_rank),
        n_samples=int(design_frame.shape[0]),
        n_sites=int(request.result_identity_metadata.shape[0]),
        residual_degrees_of_freedom=float(result.residual_degrees_of_freedom),
        variance_method="ordinary_least_squares_residual_variance",
        moderation_method=_moderation_method(
            result.empirical_bayes_method,
            robust=bool(result.empirical_bayes_robust),
            trend=bool(result.empirical_bayes_trend),
        ),
        multiple_testing_method=request.config.multiple_testing.method,
        imputation_policy=request.config.imputed_value_policy,
        missing_value_policy=(
            "reject_missing_values_before_differential_execution"
            if policy is None
            else policy.missing_values.policy
        ),
        intensity_scale=(
            "not_recorded"
            if policy is None
            else policy.statistical_testing.input_intensity_scale
        ),
        normalisation_state=request.normalisation_state,
        batch_or_covariate_terms=batch_or_covariate_terms,
        unsupported_assumptions=unsupported_assumptions,
        warnings=warnings,
    )


def _contrast_definitions_from_matrix(
    request: InterpretedDifferentialAnalysisRequest,
) -> tuple[DifferentialContrastDefinition, ...]:
    contrast_frame = cast(pd.DataFrame, request.computation_request.contrasts.frame)
    contrast_values = contrast_frame.to_numpy(dtype=float)
    coefficient_names = tuple(contrast_frame.index)
    definitions: list[DifferentialContrastDefinition] = []
    for column_index, contrast_name in enumerate(contrast_frame.columns):
        definitions.append(
            DifferentialContrastDefinition(
                name=str(contrast_name),
                numerator_condition="not_recorded",
                denominator_condition="not_recorded",
                coefficients=tuple(
                    (
                        str(coefficient_name),
                        float(contrast_values[row_index, column_index]),
                    )
                    for row_index, coefficient_name in enumerate(coefficient_names)
                ),
            )
        )
    return tuple(definitions)


def _batch_or_covariate_terms(
    execution_design: DifferentialExecutionDesignInputs | None,
) -> tuple[str, ...]:
    if execution_design is None:
        return ()
    terms: list[str] = []
    for covariate in execution_design.covariate_columns:
        terms.extend(covariate.columns)
    block_metadata = execution_design.block_column_metadata
    if block_metadata is not None:
        terms.extend(column for _, column in block_metadata.columns)
    return _unique_text(tuple(str(term) for term in terms))


def _unsupported_assumptions(
    *,
    request: InterpretedDifferentialAnalysisRequest,
    batch_or_covariate_terms: tuple[str, ...],
) -> tuple[str, ...]:
    policy = request.policy_provenance
    assumptions: list[str] = []
    if policy is not None:
        assumptions.extend(policy.unsupported_design.intentionally_rejected_features)
        assumptions.extend(policy.design.limitations)
    if batch_or_covariate_terms:
        assumptions.append(
            "fixed-effect covariates are ordinary model terms, not full batch "
            "correction or mixed-effect modelling"
        )
    batch_report = _batch_correction_report(request)
    if batch_report is not None and batch_report.status == "applied":
        assumptions.append(
            "upstream batch-correction assumptions are carried as preprocessing "
            "provenance and are not revalidated or rerun by differential analysis"
        )
    if request.ruv_readiness_enabled or request.ruv_readiness_ready:
        assumptions.append(
            "ruv_readiness metadata is report-only for differential analysis and "
            "does not enable RUV, SPS/RUV-III, duplicateCorrelation, or mixed "
            "effects"
        )
    return _unique_text(tuple(assumptions))


def _diagnostic_warnings(
    *,
    request: InterpretedDifferentialAnalysisRequest,
    batch_or_covariate_terms: tuple[str, ...],
) -> tuple[str, ...]:
    warnings: list[str] = []
    if batch_or_covariate_terms:
        warnings.append(
            "Fixed-effect covariates in this differential model are ordinary "
            "design terms; they are not full batch correction or mixed-effect "
            "modelling."
        )
    batch_report = _batch_correction_report(request)
    if batch_report is not None and batch_report.status == "applied":
        warnings.append(
            "Input dataset records upstream batch correction "
            f"method={batch_report.method!r}; differential analysis does not "
            "rerun that correction or establish limma/PhosR batch-correction "
            "parity."
        )
    if request.ruv_readiness_enabled or request.ruv_readiness_ready:
        warnings.append(
            "Dataset RUV-readiness metadata is diagnostic/report-only for "
            "differential analysis; no RUV, SPS/RUV-III, duplicateCorrelation, "
            "or mixed-effect model was fit."
        )
    return _unique_text(tuple(warnings))


def _batch_correction_report(
    request: InterpretedDifferentialAnalysisRequest,
) -> BatchCorrectionReport | None:
    preprocessing_report = request.dataset_preprocessing_report
    if preprocessing_report is None:
        return None
    return preprocessing_report.batch_correction


def _moderation_method(method: str, *, robust: bool, trend: bool) -> str:
    parts = ["empirical_bayes", str(method)]
    if robust and str(method) != "robust":
        parts.append("robust")
    if trend:
        parts.append("trend")
    return "_".join(parts)


def _unique_text(values: tuple[str, ...]) -> tuple[str, ...]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return tuple(ordered)


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
            "status_counts": _status_counts(feature_eligibility_inputs.result_status)
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
        columns=_index_snapshot(computation_request.matrix.columns),
    )
    return DifferentialComputationRequest(
        matrix=filtered_matrix,
        design=computation_request.design,
        contrasts=computation_request.contrasts,
        empirical_bayes=computation_request.empirical_bayes,
        multiple_testing_method=computation_request.multiple_testing_method,
    )


def _status_counts(result_status: pd.Series) -> dict[str, int]:
    counts = result_status.astype(str).value_counts(sort=False)
    return {str(status): int(count) for status, count in counts.items()}


def _attach_imputation_policy_metadata(
    *,
    enriched: pd.DataFrame,
    imputation_policy_inputs: DifferentialImputationPolicyInputs,
    contrast_name: str,
) -> None:
    feature_metadata = imputation_policy_inputs.feature_metadata
    result_status = imputation_policy_inputs.result_status
    if not feature_metadata.index.equals(
        enriched.index
    ) or not result_status.index.equals(enriched.index):
        raise WorkflowBoundaryError(
            seam="differential.executor.imputation_metadata_alignment",
            next_action=(
                "ensure interpreted imputation policy metadata aligns to public "
                "differential result rows"
            ),
            details={"contrast_name": contrast_name},
            message_prefix="differential workflow boundary validation failed",
        )
    imputed_cell_count = feature_metadata["imputed_cell_count"]
    observed_cell_count = feature_metadata["observed_cell_count"]
    imputed_fraction = feature_metadata["imputed_fraction"]
    enriched["imputed_cell_count"] = imputed_cell_count.to_numpy(dtype=np.int64)
    enriched["observed_cell_count"] = observed_cell_count.to_numpy(dtype=np.int64)
    enriched["imputed_fraction"] = imputed_fraction.to_numpy(dtype=float)
    enriched["imputation_policy"] = imputation_policy_inputs.policy
    enriched["imputation_fraction_threshold"] = imputation_policy_inputs.max_fraction
    enriched["result_status"] = result_status.astype(str).to_numpy()
    enriched["result_status_reason"] = (
        imputation_policy_inputs.result_status_reason.astype(str).to_numpy()
    )


def _attach_feature_eligibility_metadata(
    *,
    enriched: pd.DataFrame,
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs,
    contrast_name: str,
) -> None:
    feature_metadata = feature_eligibility_inputs.feature_metadata
    result_status = feature_eligibility_inputs.result_status
    if not feature_metadata.index.equals(
        enriched.index
    ) or not result_status.index.equals(enriched.index):
        raise WorkflowBoundaryError(
            seam="differential.executor.feature_eligibility_alignment",
            next_action=(
                "ensure feature eligibility metadata aligns to public "
                "differential result rows"
            ),
            details={"contrast_name": contrast_name},
            message_prefix="differential workflow boundary validation failed",
        )
    columns = (
        "analysed_value_count",
        "observed_value_count",
        "invalid_numeric_value_count",
        "unique_observed_value_count",
        "imputed_cell_count",
        "observed_cell_count",
        "imputed_fraction",
        "imputation_policy",
        "imputation_fraction_threshold",
        "result_status",
        "result_status_reason",
    )
    for column_name in columns:
        if column_name not in feature_metadata.columns:
            continue
        column = feature_metadata[column_name]
        enriched[column_name] = column.to_numpy()


def _expand_stat_table_to_full_index(
    table: pd.DataFrame,
    *,
    full_index: pd.Index,
) -> pd.DataFrame:
    statistic_columns = ("logFC", "t", "P.Value", "adj.P.Val")
    expanded_values = np.full(
        (int(full_index.size), len(statistic_columns)),
        np.nan,
        dtype=float,
    )
    full_positions_by_feature_id = {
        str(feature_id): position for position, feature_id in enumerate(full_index)
    }
    source_values = np.column_stack(
        [table[column_name].to_numpy(dtype=float) for column_name in statistic_columns]
    )
    for source_row_position, feature_id in enumerate(table.index):
        target_row_position = full_positions_by_feature_id[str(feature_id)]
        expanded_values[target_row_position, :] = source_values[source_row_position, :]
    return pd.DataFrame(
        expanded_values,
        index=_index_snapshot(full_index),
        columns=pd.Index(statistic_columns),
    )


def _expand_series_to_full_index(
    series: pd.Series,
    *,
    full_index: pd.Index,
) -> pd.Series:
    expanded_values = np.full(int(full_index.size), np.nan, dtype=float)
    full_positions_by_feature_id = {
        str(feature_id): position for position, feature_id in enumerate(full_index)
    }
    source_values = series.to_numpy(dtype=float)
    for source_position, feature_id in enumerate(series.index):
        target_position = full_positions_by_feature_id[str(feature_id)]
        expanded_values[target_position] = source_values[source_position]
    return pd.Series(
        expanded_values,
        index=_index_snapshot(full_index),
        name=series.name,
    )


def _index_snapshot(index: pd.Index) -> pd.Index:
    return pd.Index(list(index), name=index.name)


def _expand_prior_diagnostics_to_full_index(
    diagnostics: EmpiricalBayesPriorDiagnostics,
    *,
    full_index: pd.Index,
) -> EmpiricalBayesPriorDiagnostics:
    return EmpiricalBayesPriorDiagnostics(
        method=diagnostics.method,
        robust=diagnostics.robust,
        trend=diagnostics.trend,
        winsor_tail_p=diagnostics.winsor_tail_p,
        base_prior_variance=diagnostics.base_prior_variance,
        base_prior_degrees_of_freedom=diagnostics.base_prior_degrees_of_freedom,
        robust_outlier_count=diagnostics.robust_outlier_count,
        robust_outlier_fraction=diagnostics.robust_outlier_fraction,
        winsorized_low_count=diagnostics.winsorized_low_count,
        winsorized_high_count=diagnostics.winsorized_high_count,
        prior_variance=_expand_series_to_full_index(
            diagnostics.prior_variance,
            full_index=full_index,
        ),
        prior_degrees_of_freedom=_expand_series_to_full_index(
            diagnostics.prior_degrees_of_freedom,
            full_index=full_index,
        ),
        _assume_owned=True,
    )


def _expand_trend_diagnostics_to_full_index(
    diagnostics: MeanVarianceTrendDiagnostics | None,
    *,
    full_index: pd.Index,
) -> MeanVarianceTrendDiagnostics | None:
    if diagnostics is None:
        return None
    return MeanVarianceTrendDiagnostics(
        mean_intensity=_expand_series_to_full_index(
            diagnostics.mean_intensity,
            full_index=full_index,
        ),
        log_residual_variance=_expand_series_to_full_index(
            diagnostics.log_residual_variance,
            full_index=full_index,
        ),
        fitted_log_prior_variance=_expand_series_to_full_index(
            diagnostics.fitted_log_prior_variance,
            full_index=full_index,
        ),
        _assume_owned=True,
    )


__all__ = ["DifferentialAnalysisExecutor"]
