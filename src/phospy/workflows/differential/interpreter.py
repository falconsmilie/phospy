"""Internal interpreter for differential workflow requests."""

from __future__ import annotations

from typing import cast

import numpy as np
import pandas as pd

from phospy.contracts.configs.differential import (
    IMPUTED_VALUE_POLICY_REJECT,
    IMPUTED_VALUE_POLICY_WITHHOLD_IMPUTED_FEATURES,
)
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
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED,
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
from phospy.workflows.differential.models import (
    DifferentialBlockColumnMetadata,
    DifferentialConditionContrastVector,
    DifferentialCovariateColumnMetadata,
    DifferentialExecutionDesignInputs,
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
        matrix = resolved_dataset_view.phospho.loc[:, list(analysis_sample_ids)]
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
                    "matrix_samples": matrix_samples.astype(str).tolist(),
                    "design_samples": design_samples.astype(str).tolist(),
                },
                message_prefix="differential workflow boundary validation failed",
            )
        matrix_aligned = cast(pd.DataFrame, matrix.copy(deep=True))
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

        design_values = design_aligned.to_numpy(dtype=float)
        design_shape = cast(tuple[int, int], design_values.shape)
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

        xtx_inv = np.linalg.pinv(design_values.T @ design_values)
        contrast_values = contrasts_aligned.to_numpy(dtype=float)
        contrast_covariance = contrast_values.T @ xtx_inv @ contrast_values
        contrast_scale = np.sqrt(np.diag(contrast_covariance))
        if np.any(~np.isfinite(contrast_scale)) or np.any(contrast_scale <= 0.0):
            raise WorkflowBoundaryError(
                seam="differential.interpreter.non_estimable_contrast",
                next_action=(
                    "update contrasts so each column is estimable under the "
                    "resolved design matrix"
                ),
                details={
                    "contrast_names": contrasts_aligned.columns.astype(str).tolist(),
                },
                message_prefix="differential workflow boundary validation failed",
            )

        execution_design = _build_execution_design_inputs(
            design=resolved_design,
            contrasts=resolved_contrasts,
            design_aligned=design_aligned,
            contrasts_aligned=contrasts_aligned,
            design_build_result=resolved_design_build_result,
            paired_design_policy=cast(
                PairedDesignPolicy,
                request.config.paired_design_policy,
            ),
        )
        computation_request = DifferentialComputationRequest(
            matrix=cast(pd.DataFrame, matrix_aligned),
            design=execution_design.design_matrix,
            contrasts=execution_design.contrast_matrix,
            empirical_bayes=request.config.empirical_bayes,
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
        return InterpretedDifferentialAnalysisRequest(
            computation_request=computation_request,
            result_identity_metadata=result_identity_metadata,
            config=request.config,
            design_rank=rank,
            residual_degrees_of_freedom=residual_dof,
            policy_provenance=policy_provenance,
            workflow_provenance=resolved_workflow_provenance,
            dataset_preprocessing_report=resolved_dataset.preprocessing_report,
            execution_design=execution_design,
            imputation_policy_inputs=imputation_policy_inputs,
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
    design_matrix = DesignMatrix(cast(pd.DataFrame, design_aligned.copy(deep=True)))
    contrast_matrix = ContrastMatrix(
        cast(pd.DataFrame, contrasts_aligned.copy(deep=True))
    )
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
        vector = contrasts_aligned.loc[:, contrast.name]
        coefficients = tuple(
            (str(coefficient_name), float(vector.loc[coefficient_name]))
            for coefficient_name in contrasts_aligned.index
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
        aligned = cast(pd.DataFrame, site_metadata.loc[expected_index].copy(deep=True))
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
    identity = cast(pd.DataFrame, aligned.copy(deep=True))
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
    identity.index = pd.Index(site_key_values.tolist(), name="site_key")
    identity.loc[:, "site_key"] = site_key_values.tolist()
    identity.loc[:, "display_id"] = display_id_values.tolist()
    selected_columns = required_columns + tuple(
        column for column in optional_columns if column in aligned.columns
    )
    return cast(pd.DataFrame, identity.loc[:, list(selected_columns)])


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
    observed_mask = dataset_view.imputation_observed_mask
    if observed_mask is None:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.imputation_metadata",
            next_action=(
                "build the analysis-ready dataset through a supported imputation "
                "preprocessing path that preserves the observed-cell mask"
            ),
            message_prefix="differential workflow boundary validation failed",
        )
    try:
        aligned_mask = cast(
            pd.DataFrame,
            observed_mask.loc[list(matrix_index), list(analysis_sample_ids)],
        )
    except KeyError as exc:
        raise WorkflowBoundaryError(
            seam="differential.interpreter.imputation_metadata_alignment",
            next_action=(
                "ensure imputation observation metadata is aligned to the "
                "differential feature and sample labels"
            ),
            details={
                "feature_count": int(matrix_index.size),
                "sample_count": len(analysis_sample_ids),
            },
            message_prefix="differential workflow boundary validation failed",
        ) from exc
    if not aligned_mask.index.equals(matrix_index) or tuple(
        aligned_mask.columns.astype(str)
    ) != tuple(analysis_sample_ids):
        raise WorkflowBoundaryError(
            seam="differential.interpreter.imputation_metadata_alignment",
            next_action=(
                "ensure imputation observation metadata order exactly matches the "
                "differential execution matrix"
            ),
            message_prefix="differential workflow boundary validation failed",
        )

    observed_values = aligned_mask.to_numpy(dtype=bool)
    sample_count = int(observed_values.shape[1])
    observed_counts = observed_values.sum(axis=1).astype(np.int64)
    imputed_counts = (sample_count - observed_counts).astype(np.int64)
    imputed_fraction = imputed_counts.astype(float) / float(sample_count)
    feature_metadata = pd.DataFrame(
        {
            "imputed_cell_count": imputed_counts,
            "observed_cell_count": observed_counts,
            "imputed_fraction": imputed_fraction,
        },
        index=matrix_index.copy(),
    )
    feature_metadata.index.name = matrix_index.name

    condition_sample_ids = _condition_sample_ids_for_analysis(
        design=design,
        analysis_sample_ids=analysis_sample_ids,
    )
    statuses: list[str] = []
    for row_position in range(int(aligned_mask.shape[0])):
        observed_row = aligned_mask.iloc[row_position, :]
        if _has_insufficient_observed_samples_for_contrasts(
            observed_row=observed_row,
            condition_sample_ids=condition_sample_ids,
            contrasts=contrasts,
            minimum_condition_replicates=minimum_condition_replicates,
        ):
            statuses.append(DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED)
            continue
        if float(imputed_fraction[row_position]) > float(max_fraction):
            statuses.append(DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION)
            continue
        statuses.append(DIFFERENTIAL_RESULT_STATUS_TESTED)

    result_status = pd.Series(
        statuses,
        index=matrix_index.copy(),
        name="result_status",
    )
    testable_feature_ids = tuple(
        str(feature_id)
        for feature_id in result_status.index[
            result_status == DIFFERENTIAL_RESULT_STATUS_TESTED
        ].tolist()
    )
    return DifferentialImputationPolicyInputs(
        feature_metadata=feature_metadata,
        result_status=result_status,
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


def _has_insufficient_observed_samples_for_contrasts(
    *,
    observed_row: pd.Series,
    condition_sample_ids: dict[str, tuple[str, ...]],
    contrasts: tuple[Contrast, ...],
    minimum_condition_replicates: int,
) -> bool:
    for contrast in contrasts:
        for condition in (
            contrast.numerator_condition,
            contrast.denominator_condition,
        ):
            sample_ids = condition_sample_ids.get(condition, ())
            if not sample_ids:
                return True
            observed_count = int(
                observed_row.loc[list(sample_ids)].to_numpy(dtype=bool).sum()
            )
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
    aligned_metadata = site_metadata.reindex(matrix.index)
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
    remapped = matrix.copy(deep=True)
    remapped.index = pd.Index(site_keys.tolist(), name="site_key")
    return remapped


__all__ = ["DifferentialAnalysisInterpreter"]
