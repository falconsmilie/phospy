"""Experimental-design validation for differential workflows."""

from __future__ import annotations

from dataclasses import dataclass
from typing import cast

import pandas as pd

from phospy.contracts.configs.differential import (
    IMPUTED_VALUE_POLICY_REJECT,
    PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
    PAIRED_DESIGN_POLICY_FIXED_BLOCK,
    PAIRED_DESIGN_POLICY_REJECT,
    DifferentialImputedValuePolicy,
    PairedDesignPolicy,
)
from phospy.errors.validation import WorkflowValidationError
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.design.matrix_builder import (
    DesignMatrixBuilder,
    DesignMatrixBuildResult,
)
from phospy.science.design.models import (
    Contrast,
    ExperimentalDesign,
)
from phospy.science.differential.linear_model import DifferentialDesignDecomposition
from phospy.science.transformations.models import (
    IntensityScaleEstablishmentMode,
    IntensityScaleKind,
)
from phospy.validation.workflows.differential_design_rules import (
    ContrastFrameBuilder,
    DuplicateCorrelationDesignValidator,
    ExperimentalDesignConditionReplicateValidator,
    ExperimentalDesignContrastSetValidator,
    ExperimentalDesignFixedEffectValidator,
    ExperimentalDesignInputValidator,
    ExperimentalDesignSampleAlignmentValidator,
    FixedBlockDesignValidator,
    ResolvedDifferentialDesignMatrixValidator,
)
from phospy.validation.workflows.quantitative import (
    DIFFERENTIAL_LOG_ABUNDANCE_INPUT_CONTRACT,
    WorkflowQuantitativeInputValidator,
)

_DIFFERENTIAL_LOGFC_SCALE_ERROR_MESSAGE = (
    "Differential analysis reports logFC and therefore requires established "
    "log2-scale phospho intensities. Build the dataset with log2 preprocessing "
    "enabled, or provide an analysis-ready dataset with validated log2 "
    "intensity-scale state."
)
_DIFFERENTIAL_IMPUTED_DATA_ERROR_MESSAGE = (
    "Differential analysis does not currently treat imputed cells as observed "
    "measurements. Use a non-imputed dataset, filter features before imputation, "
    "or wait for/enable an explicit imputation-aware differential policy."
)
_DIFFERENTIAL_SUSPICIOUS_DECLARED_SCALE_ERROR_PREFIX = (
    "differential analysis rejects suspicious declared log2 intensity scale by default"
)


@dataclass(frozen=True, slots=True)
class ValidatedExperimentalDesignContract:
    """Validated design contract resolved into matrix-ready DataFrames."""

    design: ExperimentalDesign
    contrasts: tuple[Contrast, ...]
    analysis_sample_ids: tuple[str, ...]
    condition_labels: tuple[str, ...]
    design_frame: pd.DataFrame
    contrast_frame: pd.DataFrame
    design_decomposition: DifferentialDesignDecomposition
    design_build_result: DesignMatrixBuildResult | None = None


class ExperimentalDesignContractValidator:
    """Validate typed experimental design and contrast definitions."""

    def __init__(
        self,
        *,
        design_matrix_builder: DesignMatrixBuilder | None = None,
        input_validator: ExperimentalDesignInputValidator | None = None,
        contrast_validator: ExperimentalDesignContrastSetValidator | None = None,
        sample_alignment_validator: (
            ExperimentalDesignSampleAlignmentValidator | None
        ) = None,
        fixed_effect_validator: ExperimentalDesignFixedEffectValidator | None = None,
        condition_replicate_validator: (
            ExperimentalDesignConditionReplicateValidator | None
        ) = None,
        fixed_block_validator: FixedBlockDesignValidator | None = None,
        duplicate_correlation_validator: (
            DuplicateCorrelationDesignValidator | None
        ) = None,
        contrast_frame_builder: ContrastFrameBuilder | None = None,
        resolved_design_validator: (
            ResolvedDifferentialDesignMatrixValidator | None
        ) = None,
    ) -> None:
        self._design_matrix_builder = design_matrix_builder or DesignMatrixBuilder()
        self._input_validator = input_validator or ExperimentalDesignInputValidator()
        self._contrast_validator = (
            contrast_validator or ExperimentalDesignContrastSetValidator()
        )
        self._sample_alignment_validator = (
            sample_alignment_validator or ExperimentalDesignSampleAlignmentValidator()
        )
        self._fixed_effect_validator = (
            fixed_effect_validator or ExperimentalDesignFixedEffectValidator()
        )
        self._condition_replicate_validator = (
            condition_replicate_validator
            or ExperimentalDesignConditionReplicateValidator()
        )
        self._fixed_block_validator = (
            fixed_block_validator or FixedBlockDesignValidator()
        )
        self._duplicate_correlation_validator = (
            duplicate_correlation_validator or DuplicateCorrelationDesignValidator()
        )
        self._contrast_frame_builder = contrast_frame_builder or ContrastFrameBuilder()
        self._resolved_design_validator = (
            resolved_design_validator or ResolvedDifferentialDesignMatrixValidator()
        )

    def run(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        design: ExperimentalDesign,
        contrasts: tuple[Contrast, ...],
        allow_design_subset: bool,
        minimum_condition_replicates: int,
        paired_design_policy: PairedDesignPolicy = PAIRED_DESIGN_POLICY_REJECT,
        dataset_view: DatasetInternalView | None = None,
    ) -> ValidatedExperimentalDesignContract:
        input_validation = self._input_validator.run(
            dataset=dataset,
            design=design,
            allow_design_subset=allow_design_subset,
            minimum_condition_replicates=minimum_condition_replicates,
            paired_design_policy=paired_design_policy,
        )

        normalized_contrasts = self._contrast_validator.run(contrasts)
        sample_alignment = self._sample_alignment_validator.run(
            dataset=dataset,
            design=design,
            allow_design_subset=allow_design_subset,
            fixed_block_requested=input_validation.fixed_block_requested,
            dataset_view=dataset_view,
        )
        self._fixed_effect_validator.run(design=design)
        condition_validation = self._condition_replicate_validator.run(
            design=design,
            contrasts=normalized_contrasts,
            minimum_condition_replicates=minimum_condition_replicates,
        )
        known_conditions = condition_validation.condition_labels

        if input_validation.fixed_block_requested:
            self._fixed_block_validator.run(
                records=design.samples,
                contrasts=normalized_contrasts,
            )
            design_build_result = self._design_matrix_builder.run(
                design=design,
                condition_labels=known_conditions,
                paired_design_policy=PAIRED_DESIGN_POLICY_FIXED_BLOCK,
            )
        else:
            design_build_result = self._design_matrix_builder.run(
                design=design,
                condition_labels=known_conditions,
                paired_design_policy=(
                    PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
                    if input_validation.duplicate_correlation_requested
                    else PAIRED_DESIGN_POLICY_REJECT
                ),
            )
        design_frame = design_build_result.frame
        contrast_frame = self._contrast_frame_builder.run(
            coefficient_labels=tuple(str(label) for label in design_frame.columns),
            condition_labels=known_conditions,
            contrasts=normalized_contrasts,
        )
        design_decomposition = self._resolved_design_validator.run(
            design_frame=design_frame,
            contrast_frame=contrast_frame,
        )
        if input_validation.duplicate_correlation_requested:
            self._duplicate_correlation_validator.run(
                records=design.samples,
                analysis_sample_ids=sample_alignment.design_sample_ids,
                design_frame=design_frame,
                contrast_frame=contrast_frame,
                design_decomposition=design_decomposition,
                dataset=dataset,
            )
        return ValidatedExperimentalDesignContract(
            design=design,
            contrasts=normalized_contrasts,
            analysis_sample_ids=sample_alignment.design_sample_ids,
            condition_labels=known_conditions,
            design_frame=design_frame,
            contrast_frame=contrast_frame,
            design_decomposition=design_decomposition,
            design_build_result=design_build_result,
        )


class DifferentialDatasetEligibilityValidator:
    """Validate dataset quantitative-scale eligibility for differential logFC."""

    def __init__(
        self,
        *,
        quantitative_input_validator: WorkflowQuantitativeInputValidator | None = None,
    ) -> None:
        self._quantitative_input_validator = (
            quantitative_input_validator or WorkflowQuantitativeInputValidator()
        )

    def run(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        imputed_value_policy: DifferentialImputedValuePolicy = (
            IMPUTED_VALUE_POLICY_REJECT
        ),
        allow_suspicious_declared_input_scale: bool = False,
        dataset_view: DatasetInternalView | None = None,
    ) -> None:
        if not isinstance(cast(object, dataset), AnalysisReadyPhosphoDataset):
            raise WorkflowValidationError(
                "differential workflow request dataset must be AnalysisReadyPhosphoDataset"
            )
        if not isinstance(cast(object, allow_suspicious_declared_input_scale), bool):
            raise WorkflowValidationError(
                "differential workflow request allow_suspicious_declared_input_scale "
                "must be a bool"
            )
        if (
            imputed_value_policy != IMPUTED_VALUE_POLICY_REJECT
            and dataset.imputation_observation_metadata is None
        ):
            raise WorkflowValidationError(
                "differential imputed_value_policy="
                f"{imputed_value_policy!r} requires dataset-owned imputation "
                "observation metadata. Build the analysis-ready dataset through "
                "a supported imputation preprocessing path that preserves the "
                "observed-cell mask."
            )
        if (
            dataset.processing_state.missing_data.imputed
            and imputed_value_policy == IMPUTED_VALUE_POLICY_REJECT
        ):
            raise WorkflowValidationError(_DIFFERENTIAL_IMPUTED_DATA_ERROR_MESSAGE)
        phospho_scale = dataset.intensity_scale_state.phospho
        if (
            not dataset.intensity_scale_state.is_established
            or phospho_scale.kind is not IntensityScaleKind.LOG2
        ):
            raise WorkflowValidationError(_DIFFERENTIAL_LOGFC_SCALE_ERROR_MESSAGE)
        self._quantitative_input_validator.run(
            dataset=dataset,
            contract=DIFFERENTIAL_LOG_ABUNDANCE_INPUT_CONTRACT,
            context="differential workflow request dataset",
            dataset_view=dataset_view,
        )
        provenance = dataset.intensity_scale_state.establishment_provenance
        if (
            provenance is not None
            and provenance.mode is IntensityScaleEstablishmentMode.DECLARED
            and provenance.diagnostic_warnings
            and not allow_suspicious_declared_input_scale
        ):
            raise WorkflowValidationError(
                _suspicious_declared_scale_error_message(
                    first_warning=provenance.diagnostic_warnings[0]
                )
            )


__all__ = [
    "DifferentialDatasetEligibilityValidator",
    "ExperimentalDesignContractValidator",
    "ValidatedExperimentalDesignContract",
]


def _suspicious_declared_scale_error_message(*, first_warning: str) -> str:
    return (
        f"{_DIFFERENTIAL_SUSPICIOUS_DECLARED_SCALE_ERROR_PREFIX}; "
        f"first diagnostic warning: {first_warning}. "
        "recommended fix: rebuild dataset with correct input scale; "
        "apply supported log2 transformation; or explicitly set differential override "
        "if the declaration is scientifically trusted."
    )
