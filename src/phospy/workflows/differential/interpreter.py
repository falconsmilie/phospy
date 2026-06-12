"""Internal interpreter for differential workflow requests."""

from __future__ import annotations

from typing import cast

import numpy as np
import pandas as pd

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.design.matrix_builder import (
    DesignMatrixBuildResult,
    describe_fixed_effect_design,
)
from phospy.science.design.models import (
    FIXED_EFFECT_COVARIATE_KIND_BATCH,
    FIXED_EFFECT_COVARIATE_KIND_CATEGORICAL,
    FIXED_EFFECT_COVARIATE_KIND_CONTINUOUS,
    Contrast,
    ExperimentalDesign,
)
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
from phospy.workflows.differential.models import (
    DifferentialConditionContrastVector,
    DifferentialCovariateColumnMetadata,
    DifferentialExecutionDesignInputs,
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
            )
            resolved_contrasts = resolved_design_contract.contrasts
            resolved_analysis_sample_ids = resolved_design_contract.analysis_sample_ids
            resolved_design_matrix = DesignMatrix(resolved_design_contract.design_frame)
            resolved_contrast_matrix = ContrastMatrix(
                resolved_design_contract.contrast_frame
            )
            resolved_design_build_result = resolved_design_contract.design_build_result

        analysis_sample_ids = resolved_analysis_sample_ids
        matrix = resolved_dataset._borrow_phospho_frame().loc[
            :, list(analysis_sample_ids)
        ]
        matrix = _prefer_site_key_index_for_differential_results(
            matrix=matrix,
            site_metadata=resolved_dataset._borrow_site_metadata_frame(),  # pyright: ignore[reportPrivateUsage] - workflow boundary reads trusted internal dataset snapshots
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
            site_metadata=resolved_dataset._borrow_site_metadata_frame(),  # pyright: ignore[reportPrivateUsage] - workflow boundary reads trusted internal dataset snapshots
            expected_index=matrix_aligned.index,
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
        )


def _build_execution_design_inputs(
    *,
    design: ExperimentalDesign,
    contrasts: tuple[Contrast, ...],
    design_aligned: pd.DataFrame,
    contrasts_aligned: pd.DataFrame,
    design_build_result: DesignMatrixBuildResult | None,
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
        else describe_fixed_effect_design(design)
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
        ),
        sample_order=sample_order,
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
) -> str:
    if not covariate_columns:
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
