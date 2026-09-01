"""Differential workflow public result assembly."""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any, NoReturn, cast

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
    ProteinAwareDifferentialComputationResultInternalView,
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
    ProteinAwareDifferentialDiagnostics,
)
from phospy.science.differential.models.diagnostics import (
    PROTEIN_AWARE_DIFFERENTIAL_CLAIM_STATUS_EXPERIMENTAL,
    PROTEIN_AWARE_DIFFERENTIAL_CONDITION_NUMBER_SUMMARY_SCOPE,
    PROTEIN_AWARE_DIFFERENTIAL_CONDITION_NUMBER_SUMMARY_STATISTIC,
    PROTEIN_AWARE_DIFFERENTIAL_MODEL_TYPE,
)
from phospy.science.differential.models.duplicate_correlation import (
    DuplicateCorrelationWorkflowProvenance,
)
from phospy.science.differential.models.protein_aware import (
    ProteinAwareDifferentialComputationResult,
)
from phospy.science.differential.models.tables import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE,
)
from phospy.science.differential.protein_covariate_adjusted import (
    PROTEIN_AWARE_CENTERING_POLICY,
    PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME,
)
from phospy.workflows.differential.caveats import (
    build_protein_aware_result_caveats,
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
    ProteinAwareDifferentialResolvedInputs,
)
from phospy.workflows.differential.provenance import (
    build_protein_aware_policy_provenance,
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

    def run_protein_aware(
        self,
        *,
        request: InterpretedDifferentialAnalysisRequest,
        resolved_inputs: ProteinAwareDifferentialResolvedInputs,
        computation_result: ProteinAwareDifferentialComputationResult,
        workflow_provenance: Mapping[str, object],
    ) -> DifferentialAnalysisResult:
        _require_protein_aware_assembly_alignment(
            request=request,
            resolved_inputs=resolved_inputs,
            computation_result=computation_result,
        )
        full_index = request.result_identity_metadata.index
        feature_eligibility_inputs = (
            _protein_aware_feature_eligibility_after_computation(
                resolved_inputs=resolved_inputs,
                computation_result=computation_result,
            )
        )

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
            for contrast_name, table in ProteinAwareDifferentialComputationResultInternalView(
                computation_result
            ).contrast_tables.items()
        }
        contrast_tables = {
            contrast_name: _attach_result_identity_metadata(
                table=table,
                identity_metadata=request.result_identity_metadata,
                contrast_name=contrast_name,
                imputation_policy_inputs=request.imputation_policy_inputs,
                feature_eligibility_inputs=feature_eligibility_inputs,
            )
            for contrast_name, table in contrast_source_tables.items()
        }
        protein_aware_diagnostics = _build_protein_aware_diagnostics(
            request=request,
            resolved_inputs=resolved_inputs,
            computation_result=computation_result,
            feature_eligibility_inputs=feature_eligibility_inputs,
        )
        protein_aware_policy_provenance = build_protein_aware_policy_provenance(
            request=request,
            resolved_inputs=resolved_inputs,
            protein_aware_diagnostics=protein_aware_diagnostics,
        )
        policy_provenance = finalize_differential_policy_provenance(
            policy_provenance=request.policy_provenance,
            imputation_policy_inputs=request.imputation_policy_inputs,
            feature_eligibility_inputs=feature_eligibility_inputs,
            duplicate_correlation=None,
            protein_aware=protein_aware_policy_provenance,
        )
        diagnostics = _build_protein_aware_model_diagnostics(
            request=request,
            resolved_inputs=resolved_inputs,
            computation_result=computation_result,
            policy_provenance=policy_provenance,
            protein_aware_diagnostics=protein_aware_diagnostics,
        )
        caveats = finalize_differential_result_caveats(
            caveats=(
                *request.caveats,
                *build_protein_aware_result_caveats(
                    protein_aware=protein_aware_policy_provenance
                ),
            ),
            imputation_policy_inputs=request.imputation_policy_inputs,
            feature_eligibility_inputs=feature_eligibility_inputs,
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
            protein_aware_diagnostics=protein_aware_diagnostics,
            policy_provenance=policy_provenance,
            contrast_tables=contrast_tables,
            workflow_provenance=workflow_provenance,
            caveats=caveats,
            input_dataset_preprocessing_report=request.dataset_preprocessing_report,
            feature_eligibility=feature_eligibility_inputs.feature_metadata,
        )


def _require_protein_aware_assembly_alignment(
    *,
    request: InterpretedDifferentialAnalysisRequest,
    resolved_inputs: ProteinAwareDifferentialResolvedInputs,
    computation_result: ProteinAwareDifferentialComputationResult,
) -> None:
    full_index_labels = tuple(
        str(value) for value in request.result_identity_metadata.index.tolist()
    )
    if full_index_labels != resolved_inputs.full_site_ids:
        _raise_protein_aware_assembly_error(
            seam="full_index",
            next_action=(
                "assemble protein-aware results against the same full phosphosite "
                "index resolved before fitting"
            ),
            details={
                "result_index": list(full_index_labels),
                "resolved_full_site_ids": list(resolved_inputs.full_site_ids),
            },
        )
    if (
        request.execution_config.paired_design_policy
        == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
    ):
        _raise_protein_aware_assembly_error(
            seam="duplicate_correlation",
            next_action=(
                "do not assemble protein-aware differential results for "
                "duplicate-correlation requests"
            ),
        )
    request_sample_order = tuple(
        str(value) for value in request.computation_request.matrix.columns.tolist()
    )
    if resolved_inputs.sample_order != request_sample_order:
        _raise_protein_aware_assembly_error(
            seam="request_sample_order",
            next_action=(
                "assemble protein-aware results with the same sample order interpreted "
                "for ordinary differential computation"
            ),
            details={
                "request_sample_order": list(request_sample_order),
                "resolved_sample_order": list(resolved_inputs.sample_order),
            },
        )
    request_design_frame = pd.DataFrame(
        request.computation_request.design.frame,
        copy=False,
    )
    resolved_design_frame = pd.DataFrame(resolved_inputs.base_design.frame, copy=False)
    if not request_design_frame.equals(resolved_design_frame):
        _raise_protein_aware_assembly_error(
            seam="request_base_design",
            next_action=(
                "assemble protein-aware results with the exact ordinary design matrix "
                "interpreted for the request"
            ),
        )
    request_contrast_frame = pd.DataFrame(
        request.computation_request.contrasts.frame,
        copy=False,
    )
    resolved_contrast_frame = pd.DataFrame(
        resolved_inputs.base_contrasts.frame,
        copy=False,
    )
    if not request_contrast_frame.equals(resolved_contrast_frame):
        _raise_protein_aware_assembly_error(
            seam="request_base_contrasts",
            next_action=(
                "assemble protein-aware results with the exact ordinary contrast "
                "matrix interpreted for the request"
            ),
        )
    if computation_result.method_id != resolved_inputs.method_id:
        _raise_protein_aware_assembly_error(
            seam="method",
            next_action=(
                "carry one protein-aware method identifier from input resolution "
                "through computation and result assembly"
            ),
            details={
                "resolved_method_id": str(resolved_inputs.method_id),
                "computation_method_id": str(computation_result.method_id),
            },
        )
    tested_site_ids = tuple(str(value) for value in computation_result.tested_site_ids)
    expected_tested_order = tuple(
        site_id
        for site_id in resolved_inputs.tested_site_ids
        if site_id in tested_site_ids
    )
    if tested_site_ids != expected_tested_order:
        _raise_protein_aware_assembly_error(
            seam="tested_site_order",
            next_action=(
                "preserve resolved protein-aware tested phosphosite order in "
                "computation outputs"
            ),
            details={
                "computation_tested_site_ids": list(tested_site_ids),
                "expected_tested_site_ids": list(expected_tested_order),
            },
        )
    failure_site_ids = tuple(
        str(value) for value in computation_result.site_failure_diagnostics.index
    )
    unexpected_failure_ids = [
        site_id
        for site_id in failure_site_ids
        if site_id not in set(resolved_inputs.tested_site_ids)
    ]
    if unexpected_failure_ids:
        _raise_protein_aware_assembly_error(
            seam="failure_site_ids",
            next_action=(
                "report only resolved protein-aware tested phosphosites as "
                "computation-time failures"
            ),
            details={"unexpected_failure_site_ids": unexpected_failure_ids[:5]},
        )
    completed_ids = set(tested_site_ids).union(failure_site_ids)
    missing_completed_ids = [
        site_id
        for site_id in resolved_inputs.tested_site_ids
        if site_id not in completed_ids
    ]
    if missing_completed_ids:
        _raise_protein_aware_assembly_error(
            seam="completion_coverage",
            next_action=(
                "return either tested statistics or protein-aware failure "
                "diagnostics for every resolved tested phosphosite"
            ),
            details={"missing_site_ids": missing_completed_ids[:5]},
        )
    result_contrasts = set(computation_result.contrast_tables)
    expected_contrasts = set(
        str(value) for value in resolved_inputs.base_contrasts.frame.columns.tolist()
    )
    if result_contrasts != expected_contrasts:
        _raise_protein_aware_assembly_error(
            seam="contrast_tables",
            next_action=(
                "return exactly one stat-only table for each resolved base contrast"
            ),
            details={
                "result_contrasts": sorted(result_contrasts),
                "expected_contrasts": sorted(expected_contrasts),
            },
        )


def _protein_aware_feature_eligibility_after_computation(
    *,
    resolved_inputs: ProteinAwareDifferentialResolvedInputs,
    computation_result: ProteinAwareDifferentialComputationResult,
) -> DifferentialFeatureEligibilityInputs:
    feature_metadata = pd.DataFrame(
        resolved_inputs.feature_eligibility_inputs.feature_metadata,
        copy=True,
    )
    result_status = pd.Series(
        resolved_inputs.feature_eligibility_inputs.result_status,
        copy=True,
    )
    result_status = result_status.astype(str)
    if DIFFERENTIAL_RESULT_STATUS_COLUMN not in feature_metadata.columns:
        feature_metadata[DIFFERENTIAL_RESULT_STATUS_COLUMN] = result_status.to_numpy()
    if DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN not in feature_metadata.columns:
        feature_metadata[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] = ""

    successful_site_ids = set(computation_result.tested_site_ids)
    for site_id in computation_result.tested_site_ids:
        feature_metadata.loc[site_id, DIFFERENTIAL_RESULT_STATUS_COLUMN] = (
            DIFFERENTIAL_RESULT_STATUS_TESTED
        )
        feature_metadata.loc[site_id, DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] = ""
        result_status.loc[site_id] = DIFFERENTIAL_RESULT_STATUS_TESTED
        if "protein_aware_tested" in feature_metadata.columns:
            feature_metadata.loc[site_id, "protein_aware_tested"] = True

    failure_diagnostics = computation_result.site_failure_diagnostics
    for site_id, row in failure_diagnostics.iterrows():
        site_key = str(site_id)
        if site_key in successful_site_ids:
            _raise_protein_aware_assembly_error(
                seam="failure_overlap",
                next_action=(
                    "do not report a phosphosite as both successfully tested and "
                    "withheld by protein-aware computation"
                ),
                details={"site_key": site_key},
            )
        status = str(row[DIFFERENTIAL_RESULT_STATUS_COLUMN])
        reason = str(row[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN])
        feature_metadata.loc[site_key, DIFFERENTIAL_RESULT_STATUS_COLUMN] = status
        feature_metadata.loc[
            site_key,
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
        ] = reason
        result_status.loc[site_key] = status
        if "protein_aware_tested" in feature_metadata.columns:
            feature_metadata.loc[site_key, "protein_aware_tested"] = False
        if "protein_aware_failure_message" in feature_metadata.columns:
            feature_metadata.loc[site_key, "protein_aware_failure_message"] = str(
                row.get("failure_message", "")
            )

    return DifferentialFeatureEligibilityInputs(
        feature_metadata=feature_metadata,
        result_status=pd.Series(
            feature_metadata[DIFFERENTIAL_RESULT_STATUS_COLUMN].astype(str),
            index=feature_metadata.index.copy(),
            name=DIFFERENTIAL_RESULT_STATUS_COLUMN,
        ),
        testable_feature_ids=tuple(
            str(value) for value in computation_result.tested_site_ids
        ),
        attach_to_result_tables=True,
    )


def _build_protein_aware_diagnostics(
    *,
    request: InterpretedDifferentialAnalysisRequest,
    resolved_inputs: ProteinAwareDifferentialResolvedInputs,
    computation_result: ProteinAwareDifferentialComputationResult,
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs,
) -> ProteinAwareDifferentialDiagnostics:
    per_site_diagnostics = _protein_aware_per_site_diagnostics(
        resolved_inputs=resolved_inputs,
        computation_result=computation_result,
        feature_eligibility_inputs=feature_eligibility_inputs,
    )
    condition_summary = _condition_number_summary(
        computation_result.augmented_design_diagnostics,
    )
    total_site_count = int(len(resolved_inputs.full_site_ids))
    tested_site_count = int(len(computation_result.tested_site_ids))
    return ProteinAwareDifferentialDiagnostics(
        method_id=str(resolved_inputs.method_id),
        claim_status=PROTEIN_AWARE_DIFFERENTIAL_CLAIM_STATUS_EXPERIMENTAL,
        preparation_policy=resolved_inputs.preparation_policy,
        protein_mapping_policy=resolved_inputs.protein_mapping_policy,
        protein_covariate_centered=True,
        protein_covariate_standardized=False,
        protein_covariate_centering_policy=PROTEIN_AWARE_CENTERING_POLICY,
        protein_covariate_imputation_policy="none",
        fallback_policy="no_fallback_to_ordinary_differential_lane",
        total_site_count=total_site_count,
        ordinary_eligible_site_count=_count_from_pairs(
            resolved_inputs.eligibility_counts,
            key="ordinary_testable_site_count",
            default=total_site_count,
        ),
        protein_preparation_eligible_site_count=_count_from_pairs(
            resolved_inputs.eligibility_counts,
            key="protein_preparation_candidate_site_count",
            default=int(len(resolved_inputs.candidate_matched_pairs.index)),
        ),
        tested_site_count=tested_site_count,
        withheld_site_count=total_site_count - tested_site_count,
        status_counts=_status_count_items(
            feature_eligibility_inputs.feature_metadata[
                DIFFERENTIAL_RESULT_STATUS_COLUMN
            ],
        ),
        reason_counts=_reason_count_items(feature_eligibility_inputs.feature_metadata),
        execution_sample_order=resolved_inputs.sample_order,
        base_design_rank=int(request.design_decomposition.rank),
        base_residual_degrees_of_freedom=float(request.residual_degrees_of_freedom),
        expected_augmented_rank=int(request.design_decomposition.rank) + 1,
        common_augmented_rank=_common_integer_column(
            computation_result.augmented_design_diagnostics,
            column_name="rank",
        ),
        common_augmented_residual_degrees_of_freedom=float(
            computation_result.residual_degrees_of_freedom
        ),
        distinct_matched_protein_row_count=_distinct_text_count(
            resolved_inputs.candidate_matched_pairs,
            column_name="total_protein_row_key",
        ),
        fitted_protein_row_count=int(
            computation_result.augmented_design_diagnostics.shape[0]
        ),
        condition_number_summary_scope=(
            PROTEIN_AWARE_DIFFERENTIAL_CONDITION_NUMBER_SUMMARY_SCOPE
        ),
        condition_number_summary_statistic=(
            PROTEIN_AWARE_DIFFERENTIAL_CONDITION_NUMBER_SUMMARY_STATISTIC
        ),
        min_augmented_condition_number=condition_summary["min"],
        median_augmented_condition_number=condition_summary["median"],
        max_augmented_condition_number=condition_summary["max"],
        per_site_diagnostics=per_site_diagnostics,
        _assume_owned=True,
    )


def _protein_aware_per_site_diagnostics(
    *,
    resolved_inputs: ProteinAwareDifferentialResolvedInputs,
    computation_result: ProteinAwareDifferentialComputationResult,
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs,
) -> pd.DataFrame:
    full_index = pd.Index(resolved_inputs.full_site_ids, name="site_key")
    metadata = pd.DataFrame(feature_eligibility_inputs.feature_metadata, copy=False)
    site_success = pd.DataFrame(computation_result.site_diagnostics, copy=False)
    augmented_success = pd.DataFrame(
        computation_result.augmented_design_diagnostics,
        copy=False,
    )
    site_failures = pd.DataFrame(
        computation_result.site_failure_diagnostics,
        copy=False,
    )
    augmented_failures = pd.DataFrame(
        computation_result.augmented_design_failure_diagnostics,
        copy=False,
    )
    rows: list[dict[str, object]] = []
    for site_id in full_index:
        site_key = str(site_id)
        metadata_row = cast(pd.Series, metadata.loc[site_key])
        status = str(metadata_row[DIFFERENTIAL_RESULT_STATUS_COLUMN])
        reason = str(metadata_row[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN])
        total_row_key = _optional_text(
            metadata_row.get("total_protein_row_key"),
        )
        protein_identifier = _optional_text(metadata_row.get("protein_identifier"))
        success_row = (
            cast(pd.Series, site_success.loc[site_key])
            if site_key in site_success.index
            else None
        )
        if success_row is not None:
            total_row_key = _optional_text(success_row.get("total_protein_row_key"))
            protein_identifier = _optional_text(
                success_row.get("protein_identifier", protein_identifier)
            )
        failure_row = (
            cast(pd.Series, site_failures.loc[site_key])
            if site_key in site_failures.index
            else None
        )
        if failure_row is not None:
            total_row_key = _optional_text(
                failure_row.get("total_protein_row_key", total_row_key)
            )
        augmented_row = _augmented_design_row(
            total_row_key=total_row_key,
            augmented_success=augmented_success,
            augmented_failures=augmented_failures,
        )
        centered_variance = _optional_numeric(
            metadata_row.get("protein_covariate_centered_variance")
        )
        centered_std = (
            math.sqrt(centered_variance)
            if centered_variance is not None and centered_variance >= 0.0
            else None
        )
        rows.append(
            {
                "site_key": site_key,
                "protein_identifier": protein_identifier,
                "total_protein_row_key": total_row_key,
                "protein_aware_preparation_eligibility": _optional_text(
                    metadata_row.get("protein_aware_preparation_eligibility")
                ),
                "protein_aware_preparation_reasons": _json_ready_reasons(
                    metadata_row.get("protein_aware_preparation_reasons")
                ),
                "protein_covariate_raw_mean": _first_numeric(
                    metadata_row.get("protein_covariate_raw_mean"),
                    None
                    if augmented_row is None
                    else augmented_row.get("protein_raw_mean"),
                ),
                "protein_covariate_raw_standard_deviation": _first_numeric(
                    metadata_row.get("protein_covariate_raw_standard_deviation"),
                    (
                        None
                        if augmented_row is None
                        else augmented_row.get("protein_raw_standard_deviation")
                    ),
                ),
                "protein_covariate_centered_standard_deviation": _first_numeric(
                    centered_std,
                    (
                        None
                        if augmented_row is None
                        else _sqrt_optional_numeric(
                            augmented_row.get("protein_centered_variance")
                        )
                    ),
                ),
                "protein_augmented_design_rank": _first_numeric(
                    (None if augmented_row is None else augmented_row.get("rank")),
                    metadata_row.get("protein_augmented_design_rank"),
                ),
                "protein_augmented_design_residual_degrees_of_freedom": _first_numeric(
                    (
                        None
                        if augmented_row is None
                        else augmented_row.get("residual_degrees_of_freedom")
                    ),
                    metadata_row.get(
                        "protein_augmented_design_residual_degrees_of_freedom"
                    ),
                ),
                "protein_augmented_design_condition_number": _first_numeric(
                    (
                        None
                        if augmented_row is None
                        else augmented_row.get("condition_number")
                    ),
                    metadata_row.get("protein_augmented_design_condition_number"),
                ),
                "protein_contrast_estimability_status": (
                    _protein_contrast_estimability_status(status)
                ),
                DIFFERENTIAL_RESULT_STATUS_COLUMN: status,
                DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN: reason,
                "protein_covariate_coefficient": (
                    float(computation_result.protein_coefficient.loc[site_key])
                    if status == DIFFERENTIAL_RESULT_STATUS_TESTED
                    else np.nan
                ),
                "protein_aware_failure_message": _protein_failure_message(
                    metadata_row=metadata_row,
                    failure_row=failure_row,
                ),
            }
        )
    return pd.DataFrame(rows, index=full_index)


def _build_protein_aware_model_diagnostics(
    *,
    request: InterpretedDifferentialAnalysisRequest,
    resolved_inputs: ProteinAwareDifferentialResolvedInputs,
    computation_result: ProteinAwareDifferentialComputationResult,
    policy_provenance: DifferentialPolicyProvenance | None,
    protein_aware_diagnostics: ProteinAwareDifferentialDiagnostics,
) -> DifferentialModelDiagnostics:
    design_frame = resolved_inputs.base_design.frame
    augmented_design_diagnostics = computation_result.augmented_design_diagnostics
    design_columns = (
        *tuple(str(label) for label in design_frame.columns),
        PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME,
    )
    batch_or_covariate_terms = _unique_text(
        (
            *_batch_or_covariate_terms(request.execution_design),
            PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME,
        )
    )
    warnings = _unique_text(
        (
            *_diagnostic_warnings(
                request=request,
                batch_or_covariate_terms=batch_or_covariate_terms,
            ),
            "Protein-aware DifferentialModelDiagnostics.condition_number is the "
            "median condition number across independently fitted total-protein-row "
            "augmented designs; max_condition_number and rank_tolerance are base "
            "design values; singular_values is empty because there is no single "
            "shared augmented design.",
        )
    )
    unsupported_assumptions = _unique_text(
        (
            *_unsupported_assumptions(
                request=request,
                batch_or_covariate_terms=batch_or_covariate_terms,
            ),
            "protein-aware differential fitting adjusts each phosphosite for a "
            "matched mean-centered total-protein covariate and does not subtract, "
            "normalise, impute, or fall back to ordinary differential results",
        )
    )
    median_condition_number = (
        protein_aware_diagnostics.median_augmented_condition_number
    )
    return DifferentialModelDiagnostics(
        model_type=PROTEIN_AWARE_DIFFERENTIAL_MODEL_TYPE,
        design_columns=design_columns,
        contrast_definitions=_protein_aware_contrast_definitions(
            request=request,
            policy_provenance=policy_provenance,
        ),
        rank=protein_aware_diagnostics.common_augmented_rank,
        n_samples=int(design_frame.shape[0]),
        n_sites=int(len(resolved_inputs.full_site_ids)),
        residual_degrees_of_freedom=float(
            computation_result.residual_degrees_of_freedom
        ),
        decomposition_method=_common_text_column(
            augmented_design_diagnostics,
            column_name="decomposition_method",
            default=request.design_decomposition.decomposition_method,
        ),
        solver=_common_text_column(
            augmented_design_diagnostics,
            column_name="solver",
            default=request.design_decomposition.solver,
        ),
        column_scale_method=_common_text_column(
            augmented_design_diagnostics,
            column_name="column_scale_method",
            default=request.design_decomposition.column_scale_method,
        ),
        rank_tolerance_policy=_common_text_column(
            augmented_design_diagnostics,
            column_name="rank_tolerance_policy",
            default=request.design_decomposition.rank_tolerance_policy,
        ),
        rank_tolerance=request.design_decomposition.rank_tolerance,
        condition_number=(
            0.0 if median_condition_number is None else median_condition_number
        ),
        max_condition_number=request.design_decomposition.max_condition_number,
        singular_values=(),
        variance_method=(
            "protein_covariate_adjusted_ordinary_least_squares_residual_variance"
        ),
        moderation_method=_moderation_method(
            computation_result.empirical_bayes_method,
            robust=bool(computation_result.empirical_bayes_robust),
            trend=bool(computation_result.empirical_bayes_trend),
        ),
        multiple_testing_method=str(
            resolved_inputs.computation_request.multiple_testing_method
        ),
        imputation_policy=str(request.execution_config.imputed_value_policy),
        missing_value_policy=(
            "reject_missing_values_before_differential_execution"
            if policy_provenance is None
            else policy_provenance.missing_values.policy
        ),
        intensity_scale=(
            "not_recorded"
            if policy_provenance is None
            else policy_provenance.statistical_testing.input_intensity_scale
        ),
        normalisation_state=request.normalisation_state,
        batch_or_covariate_terms=batch_or_covariate_terms,
        unsupported_assumptions=unsupported_assumptions,
        warnings=warnings,
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


def _protein_aware_contrast_definitions(
    *,
    request: InterpretedDifferentialAnalysisRequest,
    policy_provenance: DifferentialPolicyProvenance | None,
) -> tuple[DifferentialContrastDefinition, ...]:
    base_definitions = (
        policy_provenance.contrasts
        if policy_provenance is not None
        else _contrast_definitions_from_matrix(request)
    )
    definitions: list[DifferentialContrastDefinition] = []
    for definition in base_definitions:
        coefficients = tuple(
            (str(coefficient), float(weight))
            for coefficient, weight in definition.coefficients
        )
        definitions.append(
            DifferentialContrastDefinition(
                name=definition.name,
                numerator_condition=definition.numerator_condition,
                denominator_condition=definition.denominator_condition,
                coefficients=(
                    *coefficients,
                    (PROTEIN_AWARE_COVARIATE_COEFFICIENT_NAME, 0.0),
                ),
                description=definition.description,
            )
        )
    return tuple(definitions)


def _condition_number_summary(
    augmented_design_diagnostics: pd.DataFrame,
) -> dict[str, float | None]:
    if "condition_number" not in augmented_design_diagnostics.columns:
        return {"min": None, "median": None, "max": None}
    values = pd.to_numeric(
        augmented_design_diagnostics["condition_number"],
        errors="coerce",
    ).to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    if not int(finite.size):
        return {"min": None, "median": None, "max": None}
    return {
        "min": float(np.min(finite)),
        "median": float(np.median(finite)),
        "max": float(np.max(finite)),
    }


def _common_integer_column(frame: pd.DataFrame, *, column_name: str) -> int:
    if column_name not in frame.columns:
        _raise_protein_aware_assembly_error(
            seam=f"augmented_design_diagnostics.{column_name}",
            next_action=(
                "return augmented design diagnostics with common rank and residual "
                "degrees of freedom fields"
            ),
        )
    values = pd.to_numeric(frame[column_name], errors="coerce").to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    if not int(finite.size):
        _raise_protein_aware_assembly_error(
            seam=f"augmented_design_diagnostics.{column_name}",
            next_action=f"record finite {column_name} values for fitted designs",
        )
    first = int(finite[0])
    if any(int(value) != first for value in finite):
        _raise_protein_aware_assembly_error(
            seam=f"augmented_design_diagnostics.{column_name}",
            next_action=(
                "record a common augmented design rank for the protein-aware result"
            ),
            details={"values": [int(value) for value in finite.tolist()]},
        )
    return first


def _common_text_column(
    frame: pd.DataFrame,
    *,
    column_name: str,
    default: str,
) -> str:
    if column_name not in frame.columns:
        return str(default)
    values = tuple(
        str(value).strip()
        for value in frame[column_name].tolist()
        if str(value).strip()
    )
    if not values:
        return str(default)
    unique_values = tuple(dict.fromkeys(values))
    if len(unique_values) == 1:
        return unique_values[0]
    return "mixed:" + ",".join(unique_values)


def _distinct_text_count(frame: pd.DataFrame, *, column_name: str) -> int:
    if column_name not in frame.columns:
        return 0
    values = [
        text
        for text in (_optional_text(value) for value in frame[column_name].tolist())
        if text is not None
    ]
    return int(len(dict.fromkeys(values)))


def _count_from_pairs(
    pairs: tuple[tuple[str, int], ...],
    *,
    key: str,
    default: int,
) -> int:
    counts = {str(name): int(count) for name, count in pairs}
    return int(counts.get(key, default))


def _status_count_items(status: pd.Series) -> tuple[tuple[str, int], ...]:
    values = [str(value) for value in status.tolist()]
    return tuple((value, values.count(value)) for value in dict.fromkeys(values))


def _reason_count_items(feature_metadata: pd.DataFrame) -> tuple[tuple[str, int], ...]:
    if DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN not in feature_metadata.columns:
        return ()
    values = [
        str(value).strip()
        for value in feature_metadata[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN].tolist()
        if str(value).strip()
    ]
    return tuple((value, values.count(value)) for value in dict.fromkeys(values))


def _optional_text(value: object) -> str | None:
    if _is_missing_scalar(value):
        return None
    text = str(value).strip()
    return text if text else None


def _optional_numeric(value: object) -> float | None:
    if _is_missing_scalar(value):
        return None
    try:
        numeric = float(cast(Any, value))
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


def _first_numeric(*values: object) -> float | None:
    for value in values:
        numeric = _optional_numeric(value)
        if numeric is not None:
            return numeric
    return None


def _sqrt_optional_numeric(value: object) -> float | None:
    numeric = _optional_numeric(value)
    if numeric is None or numeric < 0.0:
        return None
    return math.sqrt(numeric)


def _json_ready_reasons(value: object) -> tuple[str, ...]:
    if _is_missing_scalar(value):
        return ()
    if isinstance(value, str):
        text = value.strip()
        return () if not text else (text,)
    if isinstance(value, tuple):
        items = cast(tuple[object, ...], value)
        return tuple(str(item) for item in items if str(item).strip())
    if isinstance(value, list):
        items = cast(list[object], value)
        return tuple(str(item) for item in items if str(item).strip())
    return (str(value),)


def _is_missing_scalar(value: object) -> bool:
    if value is None or value is pd.NA or value is pd.NaT:
        return True
    if isinstance(value, float):
        return math.isnan(value)
    if isinstance(value, np.floating):
        return math.isnan(float(cast(Any, value)))
    return False


def _augmented_design_row(
    *,
    total_row_key: str | None,
    augmented_success: pd.DataFrame,
    augmented_failures: pd.DataFrame,
) -> pd.Series | None:
    if total_row_key is None:
        return None
    if total_row_key in augmented_success.index:
        return cast(pd.Series, augmented_success.loc[total_row_key])
    if total_row_key in augmented_failures.index:
        return cast(pd.Series, augmented_failures.loc[total_row_key])
    return None


def _protein_contrast_estimability_status(status: str) -> str:
    if status == DIFFERENTIAL_RESULT_STATUS_TESTED:
        return "estimable"
    if status == DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE:
        return "non_estimable"
    return "not_applicable"


def _protein_failure_message(
    *,
    metadata_row: pd.Series,
    failure_row: pd.Series | None,
) -> str:
    if failure_row is not None:
        text = _optional_text(failure_row.get("failure_message"))
        if text is not None:
            return text
    text = _optional_text(metadata_row.get("protein_aware_failure_message"))
    return "" if text is None else text


def _raise_protein_aware_assembly_error(
    *,
    seam: str,
    next_action: str,
    details: dict[str, object] | None = None,
) -> NoReturn:
    raise WorkflowBoundaryError(
        seam=f"differential.protein_aware_result_assembly.{seam}",
        next_action=next_action,
        details=details,
        message_prefix="differential workflow boundary validation failed",
    )


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
