"""Internal executor for differential workflow requests."""

from __future__ import annotations

from typing import cast

import numpy as np
import pandas as pd

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.differential.executor import (
    DifferentialAnalysisExecutor as DifferentialComputationExecutor,
)
from phospy.science.differential.models import (
    DifferentialAnalysisRequest as DifferentialComputationRequest,
)
from phospy.science.differential.models import (
    DifferentialAnalysisResult,
    EmpiricalBayesPriorDiagnostics,
    MeanVarianceTrendDiagnostics,
)
from phospy.workflows.differential.models import (
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
        if imputation_policy_inputs is not None:
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
        if imputation_policy_inputs is not None:
            full_index = request.result_identity_metadata.index
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
            )
            for contrast_name, table in contrast_source_tables.items()
        }
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
            policy_provenance=request.policy_provenance,
            contrast_tables=contrast_tables,
            workflow_provenance=request.workflow_provenance,
            input_dataset_preprocessing_report=request.dataset_preprocessing_report,
        )


def _attach_result_identity_metadata(
    *,
    table: pd.DataFrame,
    identity_metadata: pd.DataFrame,
    contrast_name: str,
    imputation_policy_inputs: DifferentialImputationPolicyInputs | None = None,
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
    if imputation_policy_inputs is not None:
        _attach_imputation_policy_metadata(
            enriched=enriched,
            imputation_policy_inputs=imputation_policy_inputs,
            contrast_name=contrast_name,
        )
    for column_name in ("logFC", "t", "P.Value", "adj.P.Val"):
        contrast_column = table[column_name]
        enriched[column_name] = contrast_column.to_numpy(dtype=float)
    return enriched


def _filter_computation_request_for_imputation_policy(
    *,
    computation_request: DifferentialComputationRequest,
    imputation_policy_inputs: DifferentialImputationPolicyInputs,
) -> DifferentialComputationRequest:
    testable_feature_ids = list(imputation_policy_inputs.testable_feature_ids)
    if not testable_feature_ids:
        raise WorkflowBoundaryError(
            seam="differential.executor.imputation_policy_testable_features",
            next_action=(
                "raise differential.imputed_value_max_fraction, require more "
                "observed values per condition, or use a non-imputed dataset"
            ),
            details={"policy": imputation_policy_inputs.policy},
            message_prefix="differential workflow boundary validation failed",
        )
    row_positions_by_feature_id = {
        str(feature_id): position
        for position, feature_id in enumerate(computation_request.matrix.index)
    }
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
    )


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
