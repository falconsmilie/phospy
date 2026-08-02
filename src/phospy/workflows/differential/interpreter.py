"""Internal interpreter for differential workflow requests."""

from __future__ import annotations

import pandas as pd

from phospy.contracts.configs import DifferentialAnalysisConfig
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.differential.models import (
    ContrastMatrix,
    DesignMatrix,
)
from phospy.science.differential.models import (
    DifferentialAnalysisRequest as DifferentialComputationRequest,
)
from phospy.validation.identity_contracts import (
    enforce_display_id_column,
    enforce_site_key_column,
    enforce_site_key_column_raw_matches_index,
)
from phospy.validation.workflows.differential import (
    ExperimentalDesignContractValidator,
)
from phospy.workflows._pandas_typing import (
    dataframe_copy,
    dataframe_loc,
    dataframe_reindex,
    index_as_strings,
    series_as_strings,
)
from phospy.workflows.differential.caveats import build_differential_result_caveats
from phospy.workflows.differential.design_assembly import (
    DifferentialExecutionDesignAssembler,
)
from phospy.workflows.differential.eligibility import (
    DifferentialPreFitEligibilityResolver,
    differential_status_counts,
    filter_matrix_for_feature_ids,
)
from phospy.workflows.differential.models import (
    InterpretedDifferentialAnalysisRequest,
    ResolvedDifferentialExecutionConfig,
    ValidatedDifferentialAnalysisRequest,
)
from phospy.workflows.differential.provenance import (
    build_differential_policy_provenance,
)
from phospy.workflows.differential.replicates import (
    TechnicalReplicateAggregator,
)
from phospy.workflows.intensity_scale_evidence import (
    with_input_dataset_reference_context,
    with_input_intensity_scale_evidence,
)


class DifferentialAnalysisInterpreter:
    """Resolve a validated differential request into execution-ready inputs."""

    def __init__(
        self,
        *,
        design_validator: ExperimentalDesignContractValidator | None = None,
        technical_replicate_aggregator: TechnicalReplicateAggregator | None = None,
        pre_fit_eligibility_resolver: DifferentialPreFitEligibilityResolver
        | None = None,
        execution_design_assembler: DifferentialExecutionDesignAssembler | None = None,
    ) -> None:
        self._design_validator = (
            design_validator or ExperimentalDesignContractValidator()
        )
        self._technical_replicate_aggregator = (
            technical_replicate_aggregator or TechnicalReplicateAggregator()
        )
        self._pre_fit_eligibility_resolver = (
            pre_fit_eligibility_resolver or DifferentialPreFitEligibilityResolver()
        )
        self._execution_design_assembler = (
            execution_design_assembler or DifferentialExecutionDesignAssembler()
        )

    def run(
        self, request: ValidatedDifferentialAnalysisRequest
    ) -> InterpretedDifferentialAnalysisRequest:
        aggregation_plan = request.technical_replicate_aggregation_plan
        resolved_dataset = request.dataset
        resolved_design = request.design
        resolved_contrasts = request.contrasts
        resolved_analysis_sample_ids = request.analysis_sample_ids
        resolved_design_matrix = request.design_matrix
        resolved_contrast_matrix = request.contrast_matrix
        resolved_design_decomposition = request.design_decomposition
        resolved_workflow_provenance = request.workflow_provenance
        resolved_design_build_result = request.design_build_result
        execution_config = _resolve_execution_config(request.config)

        if aggregation_plan is not None and aggregation_plan.requires_aggregation:
            technical_replicate_resolution = self._technical_replicate_aggregator.run(
                dataset=request.dataset,
                design=request.design,
                aggregation_plan=aggregation_plan,
            )
            resolved_dataset = technical_replicate_resolution.dataset
            resolved_design = technical_replicate_resolution.design
            resolved_workflow_provenance = (
                technical_replicate_resolution.workflow_provenance
            )
            resolved_design_contract = self._design_validator.run(
                dataset=resolved_dataset,
                design=resolved_design,
                contrasts=request.contrasts,
                allow_design_subset=execution_config.allow_design_subset,
                minimum_condition_replicates=(
                    execution_config.minimum_condition_replicates
                ),
                paired_design_policy=execution_config.paired_design_policy,
            )
            resolved_contrasts = resolved_design_contract.contrasts
            resolved_analysis_sample_ids = resolved_design_contract.analysis_sample_ids
            resolved_design_matrix = DesignMatrix(resolved_design_contract.design_frame)
            resolved_contrast_matrix = ContrastMatrix(
                resolved_design_contract.contrast_frame
            )
            resolved_design_decomposition = (
                resolved_design_contract.design_decomposition
            )
            resolved_design_build_result = resolved_design_contract.design_build_result

        analysis_sample_ids = resolved_analysis_sample_ids
        resolved_dataset_view = (
            request.dataset_view
            if resolved_dataset is request.dataset and request.dataset_view is not None
            else DatasetInternalView(resolved_dataset)
        )
        resolved_site_metadata = resolved_dataset_view.site_metadata
        matrix = dataframe_loc(
            resolved_dataset_view.phospho,
            columns=list(analysis_sample_ids),
        )
        matrix = _prefer_site_key_index_for_differential_results(
            matrix=matrix,
            site_metadata=resolved_site_metadata,
        )
        design_aligned = resolved_design_matrix.frame
        contrasts_aligned = resolved_contrast_matrix.frame

        matrix_samples = pd.Index(matrix.columns)
        design_samples = pd.Index(design_aligned.index)
        if not matrix_samples.equals(design_samples):
            raise WorkflowBoundaryError(
                seam="differential.interpreter.sample_label_alignment",
                next_action=(
                    "ensure validated design sample IDs exactly match the "
                    "analysis matrix sample order"
                ),
                details={
                    "matrix_samples": index_as_strings(matrix_samples),
                    "design_samples": index_as_strings(design_samples),
                },
                message_prefix="differential workflow boundary validation failed",
            )
        matrix_aligned = dataframe_copy(matrix, deep=False)
        result_identity_metadata = _build_result_identity_metadata(
            site_metadata=resolved_site_metadata,
            expected_index=matrix_aligned.index,
        )
        pre_fit_eligibility = self._pre_fit_eligibility_resolver.run(
            dataset_view=resolved_dataset_view,
            matrix=matrix_aligned,
            analysis_sample_ids=analysis_sample_ids,
            design=resolved_design,
            contrasts=resolved_contrasts,
            policy=execution_config.imputed_value_policy,
            max_fraction=execution_config.imputed_value_max_fraction,
            minimum_condition_replicates=(
                execution_config.minimum_condition_replicates
            ),
        )
        imputation_policy_inputs = pre_fit_eligibility.imputation_policy_inputs
        feature_eligibility_inputs = pre_fit_eligibility.feature_eligibility_inputs
        if not feature_eligibility_inputs.testable_feature_ids:
            raise WorkflowBoundaryError(
                seam="differential.interpreter.feature_eligibility",
                next_action=(
                    "provide at least one feature with finite, non-constant values "
                    "that satisfies the configured differential eligibility policy"
                ),
                details={
                    "status_counts": differential_status_counts(
                        feature_eligibility_inputs.result_status
                    )
                },
                message_prefix="differential workflow boundary validation failed",
            )
        matrix_for_computation = filter_matrix_for_feature_ids(
            matrix=matrix_aligned,
            feature_ids=feature_eligibility_inputs.testable_feature_ids,
        )

        rank = int(resolved_design_decomposition.rank)
        residual_dof = float(resolved_design_decomposition.residual_degrees_of_freedom)
        contrast_values = contrasts_aligned.to_numpy(dtype=float)
        invalid_contrast_positions = (
            resolved_design_decomposition.invalid_contrast_positions(contrast_values)
        )
        if invalid_contrast_positions:
            raise WorkflowBoundaryError(
                seam="differential.interpreter.non_estimable_contrast",
                next_action=(
                    "update contrasts so each column is estimable under the "
                    "resolved design matrix"
                ),
                details={
                    "contrast_names": [
                        str(contrasts_aligned.columns[int(position)])
                        for position in invalid_contrast_positions
                    ],
                },
                message_prefix="differential workflow boundary validation failed",
            )

        execution_design = self._execution_design_assembler.run(
            design=resolved_design,
            contrasts=resolved_contrasts,
            design_aligned=design_aligned,
            contrasts_aligned=contrasts_aligned,
            design_build_result=resolved_design_build_result,
            paired_design_policy=execution_config.paired_design_policy,
            design_decomposition=resolved_design_decomposition,
        )
        computation_request = DifferentialComputationRequest(
            matrix=matrix_for_computation,
            design=execution_design.design_matrix,
            contrasts=execution_design.contrast_matrix,
            design_decomposition=resolved_design_decomposition,
            empirical_bayes=execution_config.empirical_bayes,
            multiple_testing_method=execution_config.multiple_testing_method,
        )
        resolved_workflow_provenance = with_input_intensity_scale_evidence(
            resolved_workflow_provenance,
            dataset=resolved_dataset,
        )
        resolved_workflow_provenance = with_input_dataset_reference_context(
            resolved_workflow_provenance,
            dataset=resolved_dataset,
        )
        provenance_request = ValidatedDifferentialAnalysisRequest(
            dataset=resolved_dataset,
            design=resolved_design,
            contrasts=resolved_contrasts,
            analysis_sample_ids=resolved_analysis_sample_ids,
            design_matrix=resolved_design_matrix,
            contrast_matrix=resolved_contrast_matrix,
            design_decomposition=resolved_design_decomposition,
            config=request.config,
            technical_replicate_aggregation_plan=aggregation_plan,
            workflow_provenance=resolved_workflow_provenance,
            dataset_preprocessing_report=resolved_dataset.preprocessing_report,
            design_build_result=resolved_design_build_result,
        )
        policy_provenance = build_differential_policy_provenance(
            request=provenance_request,
            design_decomposition=resolved_design_decomposition,
        )
        ruv_readiness_enabled = bool(
            resolved_dataset.processing_state.ruv_readiness.enabled
        )
        ruv_readiness_ready = bool(
            resolved_dataset.processing_state.ruv_readiness.ready
        )
        caveats = build_differential_result_caveats(
            dataset=resolved_dataset,
            config=request.config,
            policy_provenance=policy_provenance,
            imputation_policy_inputs=imputation_policy_inputs,
            feature_eligibility_inputs=feature_eligibility_inputs,
            ruv_readiness_enabled=ruv_readiness_enabled,
            ruv_readiness_ready=ruv_readiness_ready,
        )
        return InterpretedDifferentialAnalysisRequest(
            computation_request=computation_request,
            result_identity_metadata=result_identity_metadata,
            config=request.config,
            execution_config=execution_config,
            design_rank=rank,
            residual_degrees_of_freedom=residual_dof,
            design_decomposition=resolved_design_decomposition,
            policy_provenance=policy_provenance,
            workflow_provenance=resolved_workflow_provenance,
            caveats=caveats,
            dataset_preprocessing_report=resolved_dataset.preprocessing_report,
            execution_design=execution_design,
            imputation_policy_inputs=imputation_policy_inputs,
            feature_eligibility_inputs=feature_eligibility_inputs,
            normalisation_state=_normalisation_state_label(resolved_dataset),
            ruv_readiness_enabled=ruv_readiness_enabled,
            ruv_readiness_ready=ruv_readiness_ready,
        )


def _build_result_identity_metadata(
    *,
    site_metadata: pd.DataFrame,
    expected_index: pd.Index,
) -> pd.DataFrame:
    required_columns = (
        "site_key",
        "display_id",
        "organism",
        "protein_namespace",
        "protein_identifier",
        "gene_symbol",
        "site",
    )
    optional_columns = (
        "protein_id",
        "protein_accession",
        "isoform_id",
    )
    missing = [
        column for column in required_columns if column not in site_metadata.columns
    ]
    if missing:
        joined = ", ".join(missing)
        raise WorkflowBoundaryError(
            seam="differential.interpreter.result_identity_columns",
            next_action=(
                "ensure resolved dataset.site_metadata includes site_key, "
                "display_id, organism, protein_namespace, protein_identifier, "
                "gene_symbol, and site"
            ),
            details={"missing_columns": missing},
            message_prefix=(
                "differential workflow boundary validation failed: "
                f"missing required result identity columns: {joined}"
            ),
        )
    try:
        aligned = dataframe_copy(
            dataframe_loc(site_metadata, rows=expected_index),
            deep=True,
        )
    except KeyError as exc:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.result_identity_alignment",
            next_action=(
                "ensure resolved dataset.site_metadata is indexed by the exact "
                "site_key labels used by the differential matrix"
            ),
            details={"expected_index_count": int(expected_index.size)},
            message_prefix="differential workflow boundary validation failed",
        ) from exc
    if not aligned.index.equals(expected_index):
        raise WorkflowBoundaryError(
            seam="differential.interpreter.result_identity_alignment",
            next_action=(
                "ensure resolved dataset.site_metadata index order exactly matches "
                "the differential matrix site_key order"
            ),
            details={
                "expected_index_count": int(expected_index.size),
                "actual_index_count": int(aligned.index.size),
            },
            message_prefix="differential workflow boundary validation failed",
        )
    identity = dataframe_copy(aligned, deep=True)
    try:
        enforce_site_key_column_raw_matches_index(
            site_metadata=identity,
            field_name="differential workflow request dataset.site_metadata",
            error_type=ValueError,
        )
        site_key_values = enforce_site_key_column(
            site_metadata=identity,
            field_name="differential workflow request dataset.site_metadata",
            error_type=ValueError,
        )
    except ValueError as exc:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.result_identity_site_key",
            next_action=(
                "ensure resolved dataset.site_metadata.site_key contains valid "
                "protein-scoped site keys"
            ),
            details={"error": str(exc)},
            message_prefix="differential workflow boundary validation failed",
        ) from exc
    try:
        display_id_values = enforce_display_id_column(
            site_metadata=identity,
            field_name="differential workflow request dataset.site_metadata",
            error_type=ValueError,
        )
    except ValueError as exc:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.result_identity_display_id",
            next_action=(
                "ensure resolved dataset.site_metadata.display_id is present for "
                "every differential matrix row"
            ),
            details={"error": str(exc)},
            message_prefix="differential workflow boundary validation failed",
        ) from exc
    site_key_list = series_as_strings(site_key_values)
    display_id_list = series_as_strings(display_id_values)
    identity.index = pd.Index(site_key_list, name="site_key")
    identity["site_key"] = site_key_list
    identity["display_id"] = display_id_list
    selected_columns = required_columns + tuple(
        column for column in optional_columns if column in aligned.columns
    )
    return dataframe_loc(identity, columns=list(selected_columns))


def _prefer_site_key_index_for_differential_results(
    *,
    matrix: pd.DataFrame,
    site_metadata: pd.DataFrame,
) -> pd.DataFrame:
    if "site_key" not in site_metadata.columns:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.result_identity_columns",
            next_action=(
                "ensure resolved dataset.site_metadata includes site_key before "
                "differential execution"
            ),
            message_prefix="differential workflow boundary validation failed",
        )
    aligned_metadata = dataframe_reindex(site_metadata, matrix.index)
    try:
        site_keys = enforce_site_key_column(
            site_metadata=aligned_metadata,
            field_name="differential workflow request dataset.site_metadata",
            error_type=ValueError,
        )
        enforce_site_key_column_raw_matches_index(
            site_metadata=aligned_metadata,
            field_name="differential workflow request dataset.site_metadata",
            error_type=ValueError,
        )
    except ValueError as exc:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.result_identity_site_key",
            next_action=(
                "ensure resolved dataset.site_metadata.site_key contains valid "
                "protein-scoped site keys"
            ),
            details={"error": str(exc)},
            message_prefix="differential workflow boundary validation failed",
        ) from exc
    remapped = dataframe_copy(matrix, deep=False)
    remapped.index = pd.Index(series_as_strings(site_keys), name="site_key")
    return remapped


def _resolve_execution_config(
    config: DifferentialAnalysisConfig,
) -> ResolvedDifferentialExecutionConfig:
    return ResolvedDifferentialExecutionConfig(
        technical_replicate_policy=config.technical_replicate_policy,
        paired_design_policy=config.paired_design_policy,
        imputed_value_policy=config.imputed_value_policy,
        imputed_value_max_fraction=config.imputed_value_max_fraction,
        allow_design_subset=config.allow_design_subset,
        allow_suspicious_declared_input_scale=(
            config.allow_suspicious_declared_input_scale
        ),
        minimum_condition_replicates=config.minimum_condition_replicates,
        empirical_bayes=config.empirical_bayes,
        multiple_testing_method=config.multiple_testing.method,
    )


def _normalisation_state_label(dataset: object) -> str:
    processing_state = getattr(dataset, "processing_state", None)
    normalisation = getattr(processing_state, "normalisation", None)
    policy = getattr(normalisation, "policy", None)
    if policy is None:
        return "not_recorded"
    value = getattr(policy, "value", policy)
    return str(value)


__all__ = ["DifferentialAnalysisInterpreter"]
