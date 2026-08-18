"""Differential workflow public result assembly."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd

from phospy.contracts.configs.differential import (
    PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
)
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.datasets.preprocessing.batch_correction import (
    BatchCorrectionReport,
)
from phospy.science.differential.internal_view import (
    DifferentialComputationResultInternalView,
)
from phospy.science.differential.models import (
    ContrastMatrix,
    DifferentialAnalysisResult,
    DifferentialComputationResult,
    DifferentialContrastDefinition,
    DifferentialModelDiagnostics,
    DifferentialPolicyProvenance,
    EmpiricalBayesPriorDiagnostics,
    MeanVarianceTrendDiagnostics,
)
from phospy.science.differential.models.duplicate_correlation import (
    DuplicateCorrelationWorkflowProvenance,
)
from phospy.workflows.differential.caveats import (
    finalize_differential_result_caveats,
)
from phospy.workflows.differential.eligibility import (
    DifferentialExecutionEligibilityResolution,
)
from phospy.workflows.differential.imputation_inference import (
    imputation_inference_columns,
)
from phospy.workflows.differential.models import (
    DifferentialExecutionDesignInputs,
    DifferentialFeatureEligibilityInputs,
    DifferentialImputationPolicyInputs,
    InterpretedDifferentialAnalysisRequest,
)
from phospy.workflows.differential.provenance import (
    finalize_differential_policy_provenance,
)


class DifferentialResultAssembler:
    """Assemble public differential workflow results from fitted outputs."""

    def run(
        self,
        *,
        request: InterpretedDifferentialAnalysisRequest,
        computation_result: DifferentialComputationResult,
        eligibility: DifferentialExecutionEligibilityResolution,
        workflow_provenance: Mapping[str, object],
        duplicate_correlation: DuplicateCorrelationWorkflowProvenance | None = None,
    ) -> DifferentialAnalysisResult:
        _require_fitted_decomposition_identity(
            request=request,
            computation_result=computation_result,
        )
        residual_variance = computation_result.residual_variance
        posterior_residual_variance = computation_result.posterior_residual_variance
        prior_residual_variance = computation_result.prior_residual_variance
        prior_degrees_of_freedom_series_value = (
            computation_result.prior_degrees_of_freedom_series_value
        )
        prior_diagnostics = computation_result.prior_diagnostics
        mean_variance_trend_diagnostics = (
            computation_result.mean_variance_trend_diagnostics
        )
        contrast_source_tables: Mapping[str, pd.DataFrame] = (
            DifferentialComputationResultInternalView(
                computation_result
            ).contrast_tables
        )
        full_index = request.result_identity_metadata.index
        if not computation_result.residual_variance.index.equals(full_index):
            residual_variance = _expand_series_to_full_index(
                computation_result.residual_variance,
                full_index=full_index,
            )
            posterior_residual_variance = _expand_series_to_full_index(
                computation_result.posterior_residual_variance,
                full_index=full_index,
            )
            prior_residual_variance = _expand_series_to_full_index(
                computation_result.prior_residual_variance,
                full_index=full_index,
            )
            prior_degrees_of_freedom_series_value = _expand_series_to_full_index(
                computation_result.prior_degrees_of_freedom_series_value,
                full_index=full_index,
            )
            prior_diagnostics = _expand_prior_diagnostics_to_full_index(
                computation_result.prior_diagnostics,
                full_index=full_index,
            )
            mean_variance_trend_diagnostics = _expand_trend_diagnostics_to_full_index(
                computation_result.mean_variance_trend_diagnostics,
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
                imputation_policy_inputs=request.imputation_policy_inputs,
                feature_eligibility_inputs=eligibility.feature_eligibility_inputs,
            )
            for contrast_name, table in contrast_source_tables.items()
        }
        policy_provenance = finalize_differential_policy_provenance(
            policy_provenance=request.policy_provenance,
            imputation_policy_inputs=request.imputation_policy_inputs,
            feature_eligibility_inputs=eligibility.feature_eligibility_inputs,
            duplicate_correlation=duplicate_correlation,
        )
        diagnostics = _build_model_diagnostics(
            request=request,
            result=computation_result,
            policy_provenance=policy_provenance,
        )
        caveats = finalize_differential_result_caveats(
            caveats=request.caveats,
            imputation_policy_inputs=request.imputation_policy_inputs,
            feature_eligibility_inputs=eligibility.feature_eligibility_inputs,
        )
        return DifferentialAnalysisResult.from_trusted_owned(
            residual_variance=residual_variance,
            posterior_residual_variance=posterior_residual_variance,
            prior_residual_variance=prior_residual_variance,
            prior_degrees_of_freedom_series_value=(
                prior_degrees_of_freedom_series_value
            ),
            prior_variance=computation_result.prior_variance,
            prior_degrees_of_freedom=computation_result.prior_degrees_of_freedom,
            residual_degrees_of_freedom=computation_result.residual_degrees_of_freedom,
            empirical_bayes_method=computation_result.empirical_bayes_method,
            empirical_bayes_robust=computation_result.empirical_bayes_robust,
            empirical_bayes_trend=computation_result.empirical_bayes_trend,
            prior_diagnostics=prior_diagnostics,
            mean_variance_trend_diagnostics=mean_variance_trend_diagnostics,
            diagnostics=diagnostics,
            policy_provenance=policy_provenance,
            contrast_tables=contrast_tables,
            workflow_provenance=workflow_provenance,
            caveats=caveats,
            input_dataset_preprocessing_report=request.dataset_preprocessing_report,
            feature_eligibility=(
                None
                if eligibility.feature_eligibility_inputs is None
                else eligibility.feature_eligibility_inputs.feature_metadata
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
    policy_provenance: DifferentialPolicyProvenance | None = None,
) -> DifferentialModelDiagnostics:
    design_frame = request.computation_request.design.frame
    policy = (
        request.policy_provenance if policy_provenance is None else policy_provenance
    )
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
    duplicate_correlation_requested = (
        request.execution_config.paired_design_policy
        == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
    )
    return DifferentialModelDiagnostics(
        model_type=(
            "moderated_gls_duplicate_correlation"
            if duplicate_correlation_requested
            else "moderated_ols_fixed_effect"
        ),
        design_columns=design_columns,
        contrast_definitions=contrast_definitions,
        rank=int(result.design_decomposition.rank),
        n_samples=int(design_frame.shape[0]),
        n_sites=int(request.result_identity_metadata.shape[0]),
        residual_degrees_of_freedom=float(
            result.design_decomposition.residual_degrees_of_freedom
        ),
        decomposition_method=result.design_decomposition.decomposition_method,
        solver=result.design_decomposition.solver,
        column_scale_method=result.design_decomposition.column_scale_method,
        rank_tolerance_policy=result.design_decomposition.rank_tolerance_policy,
        rank_tolerance=result.design_decomposition.rank_tolerance,
        condition_number=result.design_decomposition.condition_number,
        max_condition_number=result.design_decomposition.max_condition_number,
        singular_values=result.design_decomposition.singular_values,
        variance_method=(
            "compound_symmetry_gls_residual_variance"
            if duplicate_correlation_requested
            else "ordinary_least_squares_residual_variance"
        ),
        moderation_method=_moderation_method(
            result.empirical_bayes_method,
            robust=bool(result.empirical_bayes_robust),
            trend=bool(result.empirical_bayes_trend),
        ),
        multiple_testing_method=request.execution_config.multiple_testing_method,
        imputation_policy=request.execution_config.imputed_value_policy,
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


def _require_fitted_decomposition_identity(
    *,
    request: InterpretedDifferentialAnalysisRequest,
    computation_result: DifferentialComputationResult,
) -> None:
    if computation_result.design_decomposition is not request.design_decomposition:
        raise WorkflowBoundaryError(
            seam="differential.executor.fitted_decomposition_identity",
            next_action=(
                "fit differential statistics with the same design decomposition "
                "that was validated and interpreted"
            ),
            message_prefix="differential workflow boundary validation failed",
        )
    if (
        request.computation_request.design_decomposition
        is not request.design_decomposition
    ):
        raise WorkflowBoundaryError(
            seam="differential.executor.computation_decomposition_identity",
            next_action=(
                "pass the interpreted design decomposition into the computation "
                "request without rebuilding it"
            ),
            message_prefix="differential workflow boundary validation failed",
        )


def _contrast_definitions_from_matrix(
    request: InterpretedDifferentialAnalysisRequest,
) -> tuple[DifferentialContrastDefinition, ...]:
    contrasts = request.computation_request.contrasts
    if not isinstance(contrasts, ContrastMatrix):
        raise WorkflowBoundaryError(
            seam="differential.result_assembly.contrast_matrix",
            next_action=(
                "pass a validated ContrastMatrix from the differential workflow "
                "validator into result assembly"
            ),
            message_prefix="differential workflow boundary validation failed",
        )
    contrast_frame = contrasts.frame
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
        if (
            request.execution_config.paired_design_policy
            == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
        ):
            assumptions.append(
                "ruv_readiness metadata is report-only for differential analysis "
                "and does not enable RUV, SPS/RUV-III, or mixed effects"
            )
        else:
            assumptions.append(
                "ruv_readiness metadata is report-only for differential analysis and "
                "does not enable RUV, SPS/RUV-III, duplicate_correlation, or mixed "
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
        if (
            request.execution_config.paired_design_policy
            == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
        ):
            warnings.append(
                "Dataset RUV-readiness metadata is diagnostic/report-only for "
                "differential analysis; no RUV, SPS/RUV-III, or mixed-effect "
                "model was fit."
            )
        else:
            warnings.append(
                "Dataset RUV-readiness metadata is diagnostic/report-only for "
                "differential analysis; no RUV, SPS/RUV-III, duplicate_correlation, "
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
    for column_name, values in imputation_inference_columns(
        feature_metadata=feature_metadata,
        result_status=result_status,
    ).items():
        enriched[column_name] = values


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
        "contains_imputed_cells",
        "observed_only_fit",
        "residual_df_adjusted_for_imputation",
        "inferential_status",
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


__all__ = ["DifferentialResultAssembler"]
