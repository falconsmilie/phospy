"""Internal interpreter for differential workflow requests."""

from __future__ import annotations

from typing import cast

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from phospy.contracts.configs.differential import (
    IMPUTED_VALUE_POLICY_REJECT,
    IMPUTED_VALUE_POLICY_WITHHOLD_IMPUTED_FEATURES,
)
from phospy.errors.validation import DatasetValidationError
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.design.matrix_builder import (
    DesignMatrixBuildResult,
    describe_fixed_effect_design,
)
from phospy.science.design.models import (
    FIXED_EFFECT_COVARIATE_KIND_BATCH,
    FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL,
    FIXED_EFFECT_COVARIATE_KIND_CONTINUOUS,
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    Contrast,
    ExperimentalDesign,
    PairedDesignPolicy,
)
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_INVALID_NUMERIC_VALUES,
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
    dataframe_column,
    dataframe_copy,
    dataframe_loc,
    dataframe_reindex,
    index_as_strings,
    index_snapshot,
    series_as_strings,
    series_copy,
)
from phospy.workflows.differential.caveats import build_differential_result_caveats
from phospy.workflows.differential.models import (
    DifferentialBlockColumnMetadata,
    DifferentialConditionContrastVector,
    DifferentialCovariateColumnMetadata,
    DifferentialExecutionDesignInputs,
    DifferentialFeatureEligibilityInputs,
    DifferentialImputationPolicyInputs,
    InterpretedDifferentialAnalysisRequest,
    ValidatedDifferentialAnalysisRequest,
)
from phospy.workflows.differential.provenance import (
    build_differential_policy_provenance,
)
from phospy.workflows.differential.replicates import (
    TechnicalReplicateAggregator,
)
from phospy.workflows.intensity_scale_evidence import (
    with_input_intensity_scale_evidence,
)


class DifferentialAnalysisInterpreter:
    """Resolve a validated differential request into execution-ready inputs."""

    def __init__(
        self,
        *,
        design_validator: ExperimentalDesignContractValidator | None = None,
        technical_replicate_aggregator: TechnicalReplicateAggregator | None = None,
    ) -> None:
        self._design_validator = (
            design_validator or ExperimentalDesignContractValidator()
        )
        self._technical_replicate_aggregator = (
            technical_replicate_aggregator or TechnicalReplicateAggregator()
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
        resolved_workflow_provenance = request.workflow_provenance
        resolved_design_build_result = request.design_build_result

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
                allow_design_subset=request.config.allow_design_subset,
                minimum_condition_replicates=request.config.minimum_condition_replicates,
                paired_design_policy=request.config.paired_design_policy,
            )
            resolved_contrasts = resolved_design_contract.contrasts
            resolved_analysis_sample_ids = resolved_design_contract.analysis_sample_ids
            resolved_design_matrix = DesignMatrix(resolved_design_contract.design_frame)
            resolved_contrast_matrix = ContrastMatrix(
                resolved_design_contract.contrast_frame
            )
            resolved_design_build_result = resolved_design_contract.design_build_result

        analysis_sample_ids = resolved_analysis_sample_ids
        resolved_dataset_view = DatasetInternalView(resolved_dataset)
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
        matrix_aligned = dataframe_copy(matrix, deep=True)
        result_identity_metadata = _build_result_identity_metadata(
            site_metadata=resolved_site_metadata,
            expected_index=matrix_aligned.index,
        )
        imputation_policy_inputs = _build_imputation_policy_inputs(
            dataset_view=resolved_dataset_view,
            matrix_index=matrix_aligned.index,
            analysis_sample_ids=analysis_sample_ids,
            design=resolved_design,
            contrasts=resolved_contrasts,
            policy=request.config.imputed_value_policy,
            max_fraction=request.config.imputed_value_max_fraction,
            minimum_condition_replicates=request.config.minimum_condition_replicates,
        )
        feature_eligibility_inputs = _build_feature_eligibility_inputs(
            matrix=matrix_aligned,
            imputation_policy_inputs=imputation_policy_inputs,
        )
        if not feature_eligibility_inputs.testable_feature_ids:
            raise WorkflowBoundaryError(
                seam="differential.interpreter.feature_eligibility",
                next_action=(
                    "provide at least one feature with finite, non-constant values "
                    "that satisfies the configured differential eligibility policy"
                ),
                details={
                    "status_counts": _status_counts(
                        feature_eligibility_inputs.result_status
                    )
                },
                message_prefix="differential workflow boundary validation failed",
            )
        matrix_for_computation = _filter_matrix_for_feature_ids(
            matrix=matrix_aligned,
            feature_ids=feature_eligibility_inputs.testable_feature_ids,
        )

        design_values: NDArray[np.float64] = np.asarray(
            design_aligned.to_numpy(dtype=float),
            dtype=np.float64,
        )
        design_shape = design_values.shape
        sample_count = int(design_shape[0])
        coefficient_count = int(design_shape[1])
        rank = int(np.linalg.matrix_rank(design_values))
        if rank < coefficient_count:
            raise WorkflowBoundaryError(
                seam="differential.interpreter.design_rank",
                next_action=(
                    "remove collinear design terms or simplify the design matrix so "
                    "it is full column rank"
                ),
                details={"rank": rank, "columns": coefficient_count},
                message_prefix="differential workflow boundary validation failed",
            )

        residual_dof = float(sample_count - rank)
        if residual_dof <= 0.0:
            raise WorkflowBoundaryError(
                seam="differential.interpreter.residual_dof",
                next_action=(
                    "increase sample count or reduce design terms so residual "
                    "degrees of freedom stays positive"
                ),
                details={
                    "samples": sample_count,
                    "rank": rank,
                    "residual_dof": residual_dof,
                },
                message_prefix="differential workflow boundary validation failed",
            )

        design_transpose = np.transpose(design_values)
        contrast_values: NDArray[np.float64] = np.asarray(
            contrasts_aligned.to_numpy(dtype=float),
            dtype=np.float64,
        )
        contrast_transpose = np.transpose(contrast_values)
        design_crossproduct = cast(
            NDArray[np.float64],
            np.matmul(design_transpose, design_values),
        )
        xtx_inv = cast(
            NDArray[np.float64],
            np.linalg.pinv(design_crossproduct),
        )
        contrast_covariance = cast(
            NDArray[np.float64],
            np.matmul(np.matmul(contrast_transpose, xtx_inv), contrast_values),
        )
        contrast_scale = cast(
            NDArray[np.float64],
            np.sqrt(np.diagonal(contrast_covariance)),
        )
        if np.any(~np.isfinite(contrast_scale)) or np.any(contrast_scale <= 0.0):
            raise WorkflowBoundaryError(
                seam="differential.interpreter.non_estimable_contrast",
                next_action=(
                    "update contrasts so each column is estimable under the "
                    "resolved design matrix"
                ),
                details={
                    "contrast_names": index_as_strings(contrasts_aligned.columns),
                },
                message_prefix="differential workflow boundary validation failed",
            )

        execution_design = _build_execution_design_inputs(
            design=resolved_design,
            contrasts=resolved_contrasts,
            design_aligned=design_aligned,
            contrasts_aligned=contrasts_aligned,
            design_build_result=resolved_design_build_result,
            paired_design_policy=request.config.paired_design_policy,
        )
        computation_request = DifferentialComputationRequest(
            matrix=matrix_for_computation,
            design=execution_design.design_matrix,
            contrasts=execution_design.contrast_matrix,
            empirical_bayes=request.config.empirical_bayes,
            multiple_testing_method=request.config.multiple_testing.method,
        )
        resolved_workflow_provenance = with_input_intensity_scale_evidence(
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
            config=request.config,
            technical_replicate_aggregation_plan=aggregation_plan,
            workflow_provenance=resolved_workflow_provenance,
            dataset_preprocessing_report=resolved_dataset.preprocessing_report,
            design_build_result=resolved_design_build_result,
        )
        policy_provenance = build_differential_policy_provenance(
            request=provenance_request,
            design_rank=rank,
            residual_degrees_of_freedom=residual_dof,
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
            design_rank=rank,
            residual_degrees_of_freedom=residual_dof,
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


def _build_execution_design_inputs(
    *,
    design: ExperimentalDesign,
    contrasts: tuple[Contrast, ...],
    design_aligned: pd.DataFrame,
    contrasts_aligned: pd.DataFrame,
    design_build_result: DesignMatrixBuildResult | None,
    paired_design_policy: PairedDesignPolicy,
) -> DifferentialExecutionDesignInputs:
    sample_order = tuple(str(label) for label in design_aligned.index)
    coefficient_labels = tuple(str(label) for label in design_aligned.columns)
    if design_build_result is not None:
        _validate_design_build_result_alignment(
            design_build_result=design_build_result,
            sample_order=sample_order,
            coefficient_labels=coefficient_labels,
        )
    formula = (
        design_build_result.formula
        if design_build_result is not None
        else describe_fixed_effect_design(
            design,
            paired_design_policy=paired_design_policy,
        )
    )
    condition_labels = (
        design_build_result.condition_labels
        if design_build_result is not None
        else design.condition_labels()
    )
    design_matrix = DesignMatrix(dataframe_copy(design_aligned, deep=True))
    contrast_matrix = ContrastMatrix(dataframe_copy(contrasts_aligned, deep=True))
    covariate_columns = _build_covariate_column_metadata(
        design=design,
        design_build_result=design_build_result,
        coefficient_labels=coefficient_labels,
    )
    block_column_metadata = _build_block_column_metadata(
        design_build_result=design_build_result,
        coefficient_labels=coefficient_labels,
        paired_design_policy=paired_design_policy,
    )
    return DifferentialExecutionDesignInputs(
        design_matrix=design_matrix,
        contrast_matrix=contrast_matrix,
        condition_contrast_vectors=_build_condition_contrast_vectors(
            contrasts=contrasts,
            contrasts_aligned=contrasts_aligned,
        ),
        covariate_columns=covariate_columns,
        formula=formula,
        description=_execution_design_description(
            formula=formula,
            covariate_columns=covariate_columns,
            block_column_metadata=block_column_metadata,
        ),
        sample_order=sample_order,
        paired_design_policy=paired_design_policy,
        block_column_metadata=block_column_metadata,
        condition_labels=condition_labels,
        coefficient_labels=coefficient_labels,
    )


def _validate_design_build_result_alignment(
    *,
    design_build_result: DesignMatrixBuildResult,
    sample_order: tuple[str, ...],
    coefficient_labels: tuple[str, ...],
) -> None:
    build_samples = tuple(str(label) for label in design_build_result.sample_labels)
    build_coefficients = tuple(
        str(label) for label in design_build_result.coefficient_labels
    )
    if build_samples != sample_order or build_coefficients != coefficient_labels:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.design_build_metadata_alignment",
            next_action=(
                "ensure validated design build metadata describes the exact "
                "execution design matrix"
            ),
            details={
                "build_samples": list(build_samples),
                "execution_samples": list(sample_order),
                "build_coefficients": list(build_coefficients),
                "execution_coefficients": list(coefficient_labels),
            },
            message_prefix="differential workflow boundary validation failed",
        )


def _build_covariate_column_metadata(
    *,
    design: ExperimentalDesign,
    design_build_result: DesignMatrixBuildResult | None,
    coefficient_labels: tuple[str, ...],
) -> tuple[DifferentialCovariateColumnMetadata, ...]:
    modelled_covariates = tuple(
        covariate for covariate in design.fixed_effects if covariate.include_in_model
    )
    if not modelled_covariates:
        return ()
    if design_build_result is None:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.design_build_metadata_missing",
            next_action=(
                "pass validator-produced design build metadata into the interpreter "
                "for fixed-effect differential designs"
            ),
            details={
                "covariates": [covariate.name for covariate in modelled_covariates],
            },
            message_prefix="differential workflow boundary validation failed",
        )

    encoded_covariates = set(design_build_result.encoded_covariates)
    missing_encoded_covariates = [
        covariate.name
        for covariate in modelled_covariates
        if covariate.name not in encoded_covariates
    ]
    if missing_encoded_covariates:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.covariate_encoding_metadata",
            next_action=(
                "ensure design build metadata includes every modelled fixed-effect "
                "covariate"
            ),
            details={"missing_covariates": missing_encoded_covariates},
            message_prefix="differential workflow boundary validation failed",
        )

    coefficient_set = set(coefficient_labels)
    metadata: list[DifferentialCovariateColumnMetadata] = []
    for covariate in modelled_covariates:
        columns = tuple(design_build_result.covariate_columns.get(covariate.name, ()))
        if not columns:
            raise WorkflowBoundaryError(
                seam="differential.interpreter.covariate_column_metadata",
                next_action=(
                    "ensure design build metadata records execution columns for "
                    "every modelled fixed-effect covariate"
                ),
                details={"covariate": covariate.name},
                message_prefix="differential workflow boundary validation failed",
            )
        missing_columns = [
            column for column in columns if column not in coefficient_set
        ]
        if missing_columns:
            raise WorkflowBoundaryError(
                seam="differential.interpreter.covariate_column_alignment",
                next_action=(
                    "ensure covariate encoding metadata columns are present in the "
                    "execution design matrix"
                ),
                details={
                    "covariate": covariate.name,
                    "missing_columns": missing_columns,
                },
                message_prefix="differential workflow boundary validation failed",
            )

        levels: tuple[str, ...] = ()
        reference_level: str | None = None
        unused_levels: tuple[str, ...] = ()
        if covariate.kind == FIXED_EFFECT_COVARIATE_KIND_CONTINUOUS:
            pass
        elif covariate.kind in {
            FIXED_EFFECT_COVARIATE_KIND_BATCH,
            FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL,
        }:
            if (
                covariate.name not in design_build_result.categorical_levels
                or covariate.name not in design_build_result.reference_levels
                or covariate.name not in design_build_result.unused_levels
            ):
                raise WorkflowBoundaryError(
                    seam="differential.interpreter.categorical_covariate_metadata",
                    next_action=(
                        "ensure design build metadata records categorical levels "
                        "for every modelled categorical fixed-effect covariate"
                    ),
                    details={"covariate": covariate.name},
                    message_prefix="differential workflow boundary validation failed",
                )
            levels = tuple(design_build_result.categorical_levels[covariate.name])
            reference_level = design_build_result.reference_levels[covariate.name]
            unused_levels = tuple(design_build_result.unused_levels[covariate.name])
        else:
            raise WorkflowBoundaryError(
                seam="differential.interpreter.unsupported_covariate_kind",
                next_action=(
                    "validate fixed-effect covariate kinds before interpretation"
                ),
                details={"covariate": covariate.name, "kind": covariate.kind},
                message_prefix="differential workflow boundary validation failed",
            )
        metadata.append(
            DifferentialCovariateColumnMetadata(
                name=covariate.name,
                kind=covariate.kind,
                columns=columns,
                levels=levels,
                reference_level=reference_level,
                unused_levels=unused_levels,
            )
        )
    return tuple(metadata)


def _build_block_column_metadata(
    *,
    design_build_result: DesignMatrixBuildResult | None,
    coefficient_labels: tuple[str, ...],
    paired_design_policy: PairedDesignPolicy,
) -> DifferentialBlockColumnMetadata | None:
    if paired_design_policy != PAIRED_DESIGN_POLICY_FIXED_BLOCK:
        return None
    if design_build_result is None:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.design_build_metadata_missing",
            next_action=(
                "pass validator-produced design build metadata into the interpreter "
                "for fixed-block differential designs"
            ),
            message_prefix="differential workflow boundary validation failed",
        )
    if (
        not design_build_result.block_levels
        or design_build_result.block_reference_level is None
    ):
        raise WorkflowBoundaryError(
            seam="differential.interpreter.block_column_metadata",
            next_action=(
                "ensure design build metadata records fixed-block levels and "
                "reference level for fixed-block differential designs"
            ),
            message_prefix="differential workflow boundary validation failed",
        )

    coefficient_set = set(coefficient_labels)
    columns = tuple(
        (level, column)
        for level in design_build_result.block_levels
        for column in (design_build_result.block_columns.get(level),)
        if column is not None
    )
    missing_columns = [column for _, column in columns if column not in coefficient_set]
    if missing_columns:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.block_column_alignment",
            next_action=(
                "ensure fixed-block encoding metadata columns are present in the "
                "execution design matrix"
            ),
            details={"missing_columns": missing_columns},
            message_prefix="differential workflow boundary validation failed",
        )

    return DifferentialBlockColumnMetadata(
        levels=design_build_result.block_levels,
        reference_level=design_build_result.block_reference_level,
        columns=columns,
    )


def _build_condition_contrast_vectors(
    *,
    contrasts: tuple[Contrast, ...],
    contrasts_aligned: pd.DataFrame,
) -> tuple[DifferentialConditionContrastVector, ...]:
    contrast_vectors: list[DifferentialConditionContrastVector] = []
    for contrast in contrasts:
        if contrast.name not in contrasts_aligned.columns:
            raise WorkflowBoundaryError(
                seam="differential.interpreter.contrast_vector_missing",
                next_action=(
                    "ensure validated contrast matrix includes every requested "
                    "condition contrast"
                ),
                details={"contrast": contrast.name},
                message_prefix="differential workflow boundary validation failed",
            )
        vector = dataframe_column(contrasts_aligned, contrast.name)
        vector_values = vector.to_numpy(dtype=float)
        coefficients = tuple(
            (str(coefficient_name), float(vector_values[row_index]))
            for row_index, coefficient_name in enumerate(contrasts_aligned.index)
        )
        contrast_vectors.append(
            DifferentialConditionContrastVector(
                name=contrast.name,
                numerator_condition=contrast.numerator_condition,
                denominator_condition=contrast.denominator_condition,
                coefficients=coefficients,
            )
        )
    return tuple(contrast_vectors)


def _execution_design_description(
    *,
    formula: str,
    covariate_columns: tuple[DifferentialCovariateColumnMetadata, ...],
    block_column_metadata: DifferentialBlockColumnMetadata | None,
) -> str:
    if not covariate_columns and block_column_metadata is None:
        return "condition-only fixed-effect design"
    return f"fixed-effect design: {formula}"


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


def _filter_matrix_for_feature_ids(
    *,
    matrix: pd.DataFrame,
    feature_ids: tuple[str, ...],
) -> pd.DataFrame:
    return dataframe_copy(
        dataframe_loc(matrix, rows=list(feature_ids)),
        deep=True,
    )


def _status_counts(result_status: pd.Series) -> dict[str, int]:
    counts = result_status.astype(str).value_counts(sort=False)
    return {str(status): int(count) for status, count in counts.items()}


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
    remapped = dataframe_copy(matrix, deep=True)
    remapped.index = pd.Index(series_as_strings(site_keys), name="site_key")
    return remapped


def _normalisation_state_label(dataset: object) -> str:
    processing_state = getattr(dataset, "processing_state", None)
    normalisation = getattr(processing_state, "normalisation", None)
    policy = getattr(normalisation, "policy", None)
    if policy is None:
        return "not_recorded"
    value = getattr(policy, "value", policy)
    return str(value)


__all__ = ["DifferentialAnalysisInterpreter"]
