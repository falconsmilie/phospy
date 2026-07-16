"""Differential execution design assembly."""

from __future__ import annotations

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
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    Contrast,
    ExperimentalDesign,
    PairedDesignPolicy,
)
from phospy.science.differential.linear_model import DifferentialDesignDecomposition
from phospy.science.differential.models import ContrastMatrix, DesignMatrix
from phospy.workflows._pandas_typing import dataframe_column, dataframe_copy
from phospy.workflows.differential.models import (
    DifferentialBlockColumnMetadata,
    DifferentialConditionContrastVector,
    DifferentialCovariateColumnMetadata,
    DifferentialExecutionDesignInputs,
)


class DifferentialExecutionDesignAssembler:
    """Build execution-ready differential design metadata."""

    def run(
        self,
        *,
        design: ExperimentalDesign,
        contrasts: tuple[Contrast, ...],
        design_aligned: pd.DataFrame,
        contrasts_aligned: pd.DataFrame,
        design_build_result: DesignMatrixBuildResult | None,
        paired_design_policy: PairedDesignPolicy,
        design_decomposition: DifferentialDesignDecomposition,
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
            design_decomposition=design_decomposition,
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
                next_action="validate fixed-effect covariate kinds before interpretation",
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


__all__ = ["DifferentialExecutionDesignAssembler"]
