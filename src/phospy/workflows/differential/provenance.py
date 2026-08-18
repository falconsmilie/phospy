"""Structured differential-policy provenance construction."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import cast

import pandas as pd

from phospy.contracts.configs.differential import (
    IMPUTED_VALUE_POLICY_REJECT,
    IMPUTED_VALUE_POLICY_WITHHOLD_IMPUTED_FEATURES,
    PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
)
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.provenance import (
    RowAttritionRecord,
    RowAttritionReport,
    fingerprint_matrix,
    fingerprint_table,
)
from phospy.science.design.matrix_builder import (
    DesignMatrixBuildResult,
    describe_fixed_effect_design,
)
from phospy.science.design.models import ExperimentalDesign
from phospy.science.differential.compound_symmetry_gls import (
    COMPOUND_SYMMETRY_GLS_STATUS_FIT,
    CompoundSymmetryGLSFit,
)
from phospy.science.differential.linear_model import DifferentialDesignDecomposition
from phospy.science.differential.models import (
    DesignMatrix,
    DifferentialContrastDefinition,
    DifferentialDesignMatrixSummary,
    DifferentialEmpiricalBayesProvenance,
    DifferentialFixedEffectCovariateProvenance,
    DifferentialMissingValuePolicyProvenance,
    DifferentialPolicyProvenance,
    DifferentialReplicatePolicyProvenance,
    DifferentialStatisticalTestingProvenance,
    DifferentialTechnicalReplicateGroup,
    DifferentialUnsupportedDesignPolicyProvenance,
)
from phospy.science.differential.models import (
    DifferentialAnalysisRequest as DifferentialComputationRequest,
)
from phospy.science.differential.models.duplicate_correlation import (
    DUPLICATE_CORRELATION_BLOCK_TREATMENT_CONSENSUS_CORRELATION,
    DUPLICATE_CORRELATION_COVARIANCE_STRUCTURE_COMPOUND_SYMMETRY,
    DUPLICATE_CORRELATION_ESTIMATOR_FEATURE_WISE_REML,
    DUPLICATE_CORRELATION_ESTIMATOR_POLICY_VERSION,
    DUPLICATE_CORRELATION_GLS_FIT_STATUS_FIT,
    DUPLICATE_CORRELATION_WORKFLOW_PROVENANCE_VERSION,
    DuplicateCorrelationConsensusResult,
    DuplicateCorrelationWorkflowProvenance,
)
from phospy.workflows.differential.imputation_inference import (
    imputation_inference_summary_payload,
    summarize_differential_imputation_inference,
)
from phospy.workflows.differential.models import (
    DifferentialFeatureEligibilityInputs,
    DifferentialImputationPolicyInputs,
    InterpretedDifferentialAnalysisRequest,
    ValidatedDifferentialAnalysisRequest,
)
from phospy.workflows.differential.reliability import (
    resolved_minimum_condition_replicates,
    resolved_reliability_profile,
)
from phospy.workflows.intensity_scale_evidence import (
    input_intensity_scale_evidence_from_dataset,
)

_DIFFERENTIAL_TEST_STATISTIC = "moderated_t"
_DIFFERENTIAL_P_VALUE_METHOD = "two_sided_t_distribution_survival_function"
_DIFFERENTIAL_LOG2_LOGFC_INTERPRETATION = (
    "fitted condition contrast on the established log2 phosphosite intensity scale"
)
_DIFFERENTIAL_NON_LOG2_LOGFC_INTERPRETATION = (
    "fitted condition contrast on the declared input scale; not a log2 fold-change"
)
_DIFFERENTIAL_MISSING_VALUE_POLICY = (
    "reject_missing_values_before_differential_execution"
)
_DIFFERENTIAL_MISSING_VALUE_STAGE = "analysis_ready_dataset_boundary"
_DIFFERENTIAL_ADJUSTED_P_VALUE_SCOPE = (
    "adjustment_over_tested_features_only_per_contrast"
)
_DIFFERENTIAL_REJECT_IMPUTATION_LIMITATIONS: tuple[str, ...] = (
    "upstream-imputed datasets are rejected before differential execution",
)
_DIFFERENTIAL_WITHHOLD_IMPUTATION_LIMITATIONS: tuple[str, ...] = (
    (
        "features above differential.imputed_value_max_fraction are withheld "
        "from model fitting and receive missing test statistics"
    ),
    (
        "features with insufficient originally observed samples in any contrast "
        "condition are withheld from model fitting"
    ),
    (
        "Adjusted p-values are computed over tested features within each contrast "
        "only; withheld features are not part of the denominator"
    ),
    (
        "tested features are fit on the analysis-ready matrix; this policy is "
        "not observed-only fitting and does not use feature-specific residual "
        "degrees of freedom"
    ),
)
_DIFFERENTIAL_UNSUPPORTED_DESIGN_FEATURES: tuple[str, ...] = (
    (
        "correlated repeated-measure differential modelling beyond explicit "
        "fixed_block and duplicate_correlation policies"
    ),
    "mixed-effects differential modelling",
    "random subject-effect differential modelling",
)
_DIFFERENTIAL_UNSUPPORTED_ENFORCEMENT_STAGE = (
    "validation.workflows.differential.ExperimentalDesignContractValidator"
)
_DIFFERENTIAL_RANK_VALIDATION_STATUS = "validated_full_rank"
_DIFFERENTIAL_ESTIMABILITY_VALIDATION_STATUS = "validated_estimable"
_DIFFERENTIAL_CONDITIONING_VALIDATION_STATUS = "validated_scaled_svd_conditioning"
_DIFFERENTIAL_BLOCK_ID_FIELD_NAME = "block_id"
_DIFFERENTIAL_REJECT_BLOCK_CONDITION_COVERAGE_RULE = (
    "block terms are not constructed under paired_design_policy='reject'; "
    "explicit block_id values are rejected before design-matrix construction"
)
_DIFFERENTIAL_FIXED_BLOCK_CONDITION_COVERAGE_RULE = (
    "for every requested condition contrast, every block must contain both "
    "numerator and denominator conditions; incomplete or partially covered "
    "blocks are rejected before execution"
)
_DIFFERENTIAL_DUPLICATE_CORRELATION_CONDITION_COVERAGE_RULE = (
    "block_id values are retained as covariance groups, not fixed-effect "
    "columns; contrasts are validated against the non-block fixed-effects "
    "design, and at least one block must contain repeated observations"
)
_DIFFERENTIAL_UNPAIRED_LIMITATIONS: tuple[str, ...] = (
    "paired_design_policy='reject' does not construct fixed-block terms",
    (
        "explicit block_id metadata is rejected unless "
        "paired_design_policy='fixed_block' or "
        "paired_design_policy='duplicate_correlation'"
    ),
    (
        "unpaired condition and covariate workflows do not fit "
        "duplicate_correlation, mixed-effects, or random subject-effect models"
    ),
)
_DIFFERENTIAL_DUPLICATE_CORRELATION_LIMITATIONS: tuple[str, ...] = (
    (
        "duplicate_correlation estimates one consensus compound-symmetry "
        "within-block correlation"
    ),
    "duplicate_correlation does not fit feature-specific random effects",
    "duplicate_correlation does not add block_id levels as fixed-effect columns",
    (
        "singleton and incomplete blocks may contribute to fixed-effect fitting, "
        "but only repeated blocks inform within-block correlation"
    ),
)
_DIFFERENTIAL_DUPLICATE_CORRELATION_UNSUPPORTED_DESIGN_FEATURES: tuple[str, ...] = (
    "feature-specific final random-effects differential fits",
    "random-slope or multi-random-effect differential modelling",
    "mixed-effects differential modelling beyond one consensus block correlation",
)
_DIFFERENTIAL_FIXED_BLOCK_LIMITATIONS: tuple[str, ...] = (
    "fixed_block adds block_id levels as ordinary fixed-effect design columns",
    (
        "fixed_block does not estimate within-block correlation and is not "
        "limma duplicateCorrelation"
    ),
    "fixed_block does not fit mixed-effects or random subject-effect models",
    (
        "incomplete blocks are rejected before execution; samples are not "
        "silently dropped"
    ),
)


@dataclass(frozen=True, slots=True)
class _ImputationInferenceProvenanceFields:
    tested_feature_count: int
    withheld_feature_count: int
    tested_imputed_feature_count: int
    tested_imputed_cell_count: int
    observed_only_fit: bool
    residual_df_adjusted_for_imputation: bool
    inferential_status: str
    adjusted_p_value_denominator_feature_count: int


def build_differential_policy_provenance(
    *,
    request: ValidatedDifferentialAnalysisRequest,
    design_decomposition: DifferentialDesignDecomposition,
    imputation_policy_inputs: DifferentialImputationPolicyInputs | None = None,
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs | None = None,
) -> DifferentialPolicyProvenance:
    """Build deterministic structured differential-policy provenance records."""

    if design_decomposition is not request.design_decomposition:
        raise WorkflowBoundaryError(
            seam="differential.provenance.design_decomposition_identity",
            next_action=(
                "assemble differential policy provenance from the same design "
                "decomposition object produced by validation"
            ),
            message_prefix="differential workflow boundary validation failed",
        )

    design_frame = request.design_matrix.frame
    contrast_frame = request.contrast_matrix.frame
    sample_labels = tuple(str(label) for label in design_frame.index)
    coefficient_labels = tuple(str(label) for label in design_frame.columns)
    condition_columns = _condition_columns(
        coefficient_labels=coefficient_labels,
        design=request.design,
    )
    covariate_provenance = _fixed_effect_covariate_provenance(request)
    design_formula = _design_formula(request)
    input_intensity_scale = _input_intensity_scale_label(request)
    input_intensity_scale_evidence = input_intensity_scale_evidence_from_dataset(
        request.dataset
    )
    imputation_inference = _imputation_inference_provenance_fields(
        imputation_policy_inputs=imputation_policy_inputs,
        feature_eligibility_inputs=feature_eligibility_inputs,
    )

    contrast_definitions: list[DifferentialContrastDefinition] = []
    for contrast in request.contrasts:
        vector = contrast_frame.loc[:, contrast.name]
        coefficients = tuple(
            (str(coefficient_name), float(vector.loc[coefficient_name]))
            for coefficient_name in contrast_frame.index
        )
        contrast_definitions.append(
            DifferentialContrastDefinition(
                name=contrast.name,
                numerator_condition=contrast.numerator_condition,
                denominator_condition=contrast.denominator_condition,
                coefficients=coefficients,
                description=_contrast_description(
                    numerator_condition=contrast.numerator_condition,
                    denominator_condition=contrast.denominator_condition,
                ),
            )
        )

    return DifferentialPolicyProvenance(
        design=DifferentialDesignMatrixSummary(
            formula=design_formula,
            sample_labels=sample_labels,
            coefficient_labels=coefficient_labels,
            sample_count=len(sample_labels),
            coefficient_count=len(coefficient_labels),
            rank=int(design_decomposition.rank),
            residual_degrees_of_freedom=float(
                design_decomposition.residual_degrees_of_freedom
            ),
            decomposition_method=design_decomposition.decomposition_method,
            solver=design_decomposition.solver,
            column_scale_method=design_decomposition.column_scale_method,
            rank_tolerance_policy=design_decomposition.rank_tolerance_policy,
            rank_tolerance=design_decomposition.rank_tolerance,
            condition_number=design_decomposition.condition_number,
            max_condition_number=design_decomposition.max_condition_number,
            singular_values=design_decomposition.singular_values,
            description=_design_description(
                formula=design_formula,
                covariates=covariate_provenance,
            ),
            condition_columns=condition_columns,
            covariates=covariate_provenance,
            paired_design_policy=request.config.paired_design_policy,
            block_id_field_name=_DIFFERENTIAL_BLOCK_ID_FIELD_NAME,
            block_count=_block_count(request),
            block_levels=_block_levels(request),
            block_levels_included=_block_levels(request),
            block_reference_level=_block_reference_level(request),
            block_columns=_block_columns(request.design_build_result),
            block_column_names=_block_column_names(request.design_build_result),
            condition_coverage_rule=_condition_coverage_rule(
                request.config.paired_design_policy
            ),
            limitations=_design_limitations(request.config.paired_design_policy),
            rank_validation_status=_DIFFERENTIAL_RANK_VALIDATION_STATUS,
            conditioning_validation_status=(
                _DIFFERENTIAL_CONDITIONING_VALIDATION_STATUS
            ),
            estimability_validation_status=(
                _DIFFERENTIAL_ESTIMABILITY_VALIDATION_STATUS
            ),
        ),
        contrasts=tuple(contrast_definitions),
        replicates=DifferentialReplicatePolicyProvenance(
            minimum_condition_replicates=resolved_minimum_condition_replicates(
                request.config
            ),
            reliability_profile=resolved_reliability_profile(request.config),
            technical_replicate_policy=(
                request.config.technical_replicate_policy.value
            ),
            condition_replicate_counts=_condition_replicate_counts(request.design),
            technical_replicate_groups=_technical_replicate_groups(
                request.workflow_provenance
            ),
        ),
        empirical_bayes=DifferentialEmpiricalBayesProvenance(
            method=request.config.empirical_bayes.method,
            robust=request.config.empirical_bayes.method == "robust",
            trend=request.config.empirical_bayes.trend,
            winsor_tail_p=request.config.empirical_bayes.winsor_tail_p,
        ),
        statistical_testing=DifferentialStatisticalTestingProvenance(
            test_statistic=_DIFFERENTIAL_TEST_STATISTIC,
            p_value_method=_DIFFERENTIAL_P_VALUE_METHOD,
            adjusted_p_value_method=request.config.multiple_testing.method,
            input_intensity_scale=input_intensity_scale,
            input_intensity_scale_evidence_level=(
                input_intensity_scale_evidence.input_intensity_scale_evidence_level
            ),
            input_intensity_scale_source=(
                input_intensity_scale_evidence.input_intensity_scale_source
            ),
            logfc_interpretation=_logfc_interpretation(input_intensity_scale),
            allow_suspicious_declared_input_scale=(
                request.config.allow_suspicious_declared_input_scale
            ),
        ),
        missing_values=DifferentialMissingValuePolicyProvenance(
            policy=_DIFFERENTIAL_MISSING_VALUE_POLICY,
            stage=_DIFFERENTIAL_MISSING_VALUE_STAGE,
            imputed_value_policy=request.config.imputed_value_policy,
            imputed_value_max_fraction=request.config.imputed_value_max_fraction,
            imputation_metadata_required=(
                request.config.imputed_value_policy != IMPUTED_VALUE_POLICY_REJECT
            ),
            adjusted_p_value_scope=_DIFFERENTIAL_ADJUSTED_P_VALUE_SCOPE,
            tested_feature_count=imputation_inference.tested_feature_count,
            withheld_feature_count=imputation_inference.withheld_feature_count,
            tested_imputed_feature_count=(
                imputation_inference.tested_imputed_feature_count
            ),
            tested_imputed_cell_count=(imputation_inference.tested_imputed_cell_count),
            observed_only_fit=imputation_inference.observed_only_fit,
            residual_df_adjusted_for_imputation=(
                imputation_inference.residual_df_adjusted_for_imputation
            ),
            inferential_status=imputation_inference.inferential_status,
            adjusted_p_value_denominator_feature_count=(
                imputation_inference.adjusted_p_value_denominator_feature_count
            ),
            limitations=_imputation_policy_limitations(
                request.config.imputed_value_policy
            ),
        ),
        unsupported_design=DifferentialUnsupportedDesignPolicyProvenance(
            intentionally_rejected_features=_unsupported_design_features(
                request.config.paired_design_policy
            ),
            enforcement_stage=_DIFFERENTIAL_UNSUPPORTED_ENFORCEMENT_STAGE,
        ),
    )


def finalize_differential_policy_provenance(
    *,
    policy_provenance: DifferentialPolicyProvenance | None,
    imputation_policy_inputs: DifferentialImputationPolicyInputs | None,
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs | None,
    duplicate_correlation: DuplicateCorrelationWorkflowProvenance | None = None,
) -> DifferentialPolicyProvenance | None:
    """Refresh imputation inference counts after final row eligibility."""

    if policy_provenance is None:
        return None
    imputation_inference = _imputation_inference_provenance_fields(
        imputation_policy_inputs=imputation_policy_inputs,
        feature_eligibility_inputs=feature_eligibility_inputs,
    )
    return replace(
        policy_provenance,
        missing_values=replace(
            policy_provenance.missing_values,
            tested_feature_count=imputation_inference.tested_feature_count,
            withheld_feature_count=imputation_inference.withheld_feature_count,
            tested_imputed_feature_count=(
                imputation_inference.tested_imputed_feature_count
            ),
            tested_imputed_cell_count=(imputation_inference.tested_imputed_cell_count),
            observed_only_fit=imputation_inference.observed_only_fit,
            residual_df_adjusted_for_imputation=(
                imputation_inference.residual_df_adjusted_for_imputation
            ),
            inferential_status=imputation_inference.inferential_status,
            adjusted_p_value_denominator_feature_count=(
                imputation_inference.adjusted_p_value_denominator_feature_count
            ),
        ),
        duplicate_correlation=(
            duplicate_correlation
            if duplicate_correlation is not None
            else policy_provenance.duplicate_correlation
        ),
    )


def build_duplicate_correlation_workflow_provenance(
    *,
    request: InterpretedDifferentialAnalysisRequest,
    computation_request: DifferentialComputationRequest,
    consensus_result: DuplicateCorrelationConsensusResult,
    gls_fit: CompoundSymmetryGLSFit,
    imputation_policy_inputs: DifferentialImputationPolicyInputs | None = None,
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs | None = None,
) -> DuplicateCorrelationWorkflowProvenance:
    """Build typed result provenance for a completed duplicate-correlation fit."""

    execution_design = request.execution_design
    if execution_design is None or execution_design.block_ids is None:
        raise WorkflowBoundaryError(
            seam="differential.provenance.duplicate_correlation_blocks",
            next_action=(
                "carry validated duplicate-correlation block IDs through "
                "interpretation before provenance assembly"
            ),
            message_prefix="differential workflow boundary validation failed",
        )
    block_structure = consensus_result.block_structure
    convergence_summary = consensus_result.convergence_summary
    boundary_summary = consensus_result.boundary_summary
    if (
        block_structure is None
        or convergence_summary is None
        or boundary_summary is None
    ):
        raise WorkflowBoundaryError(
            seam="differential.provenance.duplicate_correlation_estimator_summary",
            next_action=(
                "return block, convergence, and boundary summaries from the "
                "duplicate-correlation estimator"
            ),
            message_prefix="differential workflow boundary validation failed",
        )
    failed_gls_statuses = tuple(
        status
        for status in gls_fit.feature_fit_statuses
        if str(status) != COMPOUND_SYMMETRY_GLS_STATUS_FIT
    )
    if failed_gls_statuses:
        raise WorkflowBoundaryError(
            seam="differential.provenance.duplicate_correlation_gls_status",
            next_action=(
                "only assemble duplicate-correlation result provenance after a "
                "successful GLS fit"
            ),
            details={"failed_gls_statuses": sorted(set(failed_gls_statuses))},
            message_prefix="differential workflow boundary validation failed",
        )
    block_ids = execution_design.block_ids
    sample_order = execution_design.sample_order
    block_assignment = pd.DataFrame(
        {
            "sample_id": list(sample_order),
            "block_id": list(block_ids),
        },
        index=pd.Index(sample_order, name="sample"),
    )
    block_counts = Counter(block_ids)
    imputation_fields = _imputation_inference_provenance_fields(
        imputation_policy_inputs=imputation_policy_inputs,
        feature_eligibility_inputs=feature_eligibility_inputs,
    )
    design_matrix = cast(DesignMatrix, computation_request.design)
    analysis_matrix_fingerprint = fingerprint_matrix(
        computation_request.matrix,
        name="differential.analysis_matrix",
    )
    return DuplicateCorrelationWorkflowProvenance(
        model="duplicate_correlation",
        provenance_version=DUPLICATE_CORRELATION_WORKFLOW_PROVENANCE_VERSION,
        requested_paired_design_policy=request.config.paired_design_policy,
        normalised_paired_design_policy=request.execution_config.paired_design_policy,
        block_treatment=DUPLICATE_CORRELATION_BLOCK_TREATMENT_CONSENSUS_CORRELATION,
        covariance_structure=DUPLICATE_CORRELATION_COVARIANCE_STRUCTURE_COMPOUND_SYMMETRY,
        estimator=DUPLICATE_CORRELATION_ESTIMATOR_FEATURE_WISE_REML,
        estimator_policy_version=DUPLICATE_CORRELATION_ESTIMATOR_POLICY_VERSION,
        trim_fraction=consensus_result.trim_fraction,
        matrix_authority="workflow approved differential analysis matrix",
        analysis_matrix_fingerprint=analysis_matrix_fingerprint,
        authoritative_matrix_fingerprint=analysis_matrix_fingerprint,
        design_authority="validation.workflows.differential non-block design",
        design_fingerprint=fingerprint_table(
            design_matrix.to_dataframe(),
            name="differential.non_block_fixed_effect_design",
        ),
        block_authority="validation.workflows.differential block_id",
        block_assignment_fingerprint=fingerprint_table(
            block_assignment,
            name="differential.duplicate_correlation_block_assignment",
        ),
        estimator_authority="science.differential duplicate-correlation REML",
        gls_authority="science.differential compound-symmetry GLS",
        failure_authority="validation, REML estimator, and GLS typed failures",
        block_structure=block_structure,
        consensus=consensus_result.to_summary(),
        attempted_feature_count=int(consensus_result.attempted_feature_count or 0),
        trimmed_feature_count_each_tail=(
            consensus_result.trimmed_feature_count_each_tail
        ),
        retained_feature_count_after_trimming=int(
            consensus_result.retained_feature_count_after_trimming or 0
        ),
        failure_reason_counts=consensus_result.failure_reason_counts,
        convergence_summary=convergence_summary,
        boundary_summary=boundary_summary,
        sample_count=int(block_structure.sample_count),
        block_count=int(block_structure.block_count),
        repeated_block_count=int(block_structure.repeated_block_count),
        singleton_block_count=int(block_structure.singleton_block_count),
        minimum_block_size=int(
            block_structure.minimum_block_size or min(block_counts.values())
        ),
        maximum_block_size=int(
            block_structure.maximum_block_size or max(block_counts.values())
        ),
        design_rank=int(consensus_result.design_rank or request.design_rank),
        gls_fit_status=DUPLICATE_CORRELATION_GLS_FIT_STATUS_FIT,
        imputed_values_participated=bool(
            imputation_fields.tested_imputed_cell_count > 0
        ),
        imputed_feature_count=int(imputation_fields.tested_imputed_feature_count),
        imputed_cell_count=int(imputation_fields.tested_imputed_cell_count),
    )


def _input_intensity_scale_label(request: ValidatedDifferentialAnalysisRequest) -> str:
    return str(request.dataset.intensity_scale_state.label)


def _logfc_interpretation(input_intensity_scale: str) -> str:
    if input_intensity_scale == "log2":
        return _DIFFERENTIAL_LOG2_LOGFC_INTERPRETATION
    return _DIFFERENTIAL_NON_LOG2_LOGFC_INTERPRETATION


def _design_formula(request: ValidatedDifferentialAnalysisRequest) -> str:
    if request.design_build_result is not None:
        return request.design_build_result.formula
    return describe_fixed_effect_design(
        request.design,
        paired_design_policy=request.config.paired_design_policy,
    )


def _imputation_policy_limitations(policy: str) -> tuple[str, ...]:
    if policy == IMPUTED_VALUE_POLICY_WITHHOLD_IMPUTED_FEATURES:
        return _DIFFERENTIAL_WITHHOLD_IMPUTATION_LIMITATIONS
    return _DIFFERENTIAL_REJECT_IMPUTATION_LIMITATIONS


def _imputation_inference_provenance_fields(
    *,
    imputation_policy_inputs: DifferentialImputationPolicyInputs | None,
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs | None,
) -> _ImputationInferenceProvenanceFields:
    if imputation_policy_inputs is None:
        return _ImputationInferenceProvenanceFields(
            tested_feature_count=0,
            withheld_feature_count=0,
            tested_imputed_feature_count=0,
            tested_imputed_cell_count=0,
            observed_only_fit=False,
            residual_df_adjusted_for_imputation=False,
            inferential_status="not_applicable",
            adjusted_p_value_denominator_feature_count=0,
        )
    summary = summarize_differential_imputation_inference(
        imputation_policy_inputs=imputation_policy_inputs,
        feature_eligibility_inputs=feature_eligibility_inputs,
    )
    return _ImputationInferenceProvenanceFields(
        tested_feature_count=int(summary.tested_feature_count),
        withheld_feature_count=int(summary.withheld_feature_count),
        tested_imputed_feature_count=int(summary.tested_imputed_feature_count),
        tested_imputed_cell_count=int(summary.tested_imputed_cell_count),
        observed_only_fit=bool(summary.observed_only_fit),
        residual_df_adjusted_for_imputation=bool(
            summary.residual_df_adjusted_for_imputation
        ),
        inferential_status=summary.inferential_status,
        adjusted_p_value_denominator_feature_count=int(
            summary.adjusted_p_value_denominator_feature_count
        ),
    )


def _block_levels(
    request: ValidatedDifferentialAnalysisRequest,
) -> tuple[str, ...]:
    if (
        request.config.paired_design_policy
        == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
    ):
        block_ids = tuple(
            str(record.block_id)
            for record in request.design.samples
            if record.sample_id in set(request.analysis_sample_ids)
            and record.block_id is not None
        )
        return tuple(sorted(set(block_ids)))
    design_build_result = request.design_build_result
    if design_build_result is None:
        return ()
    return design_build_result.block_levels


def _block_count(
    request: ValidatedDifferentialAnalysisRequest,
) -> int:
    return len(_block_levels(request))


def _block_reference_level(
    request: ValidatedDifferentialAnalysisRequest,
) -> str | None:
    if (
        request.config.paired_design_policy
        == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
    ):
        return None
    design_build_result = request.design_build_result
    if design_build_result is None:
        return None
    return design_build_result.block_reference_level


def _block_columns(
    design_build_result: DesignMatrixBuildResult | None,
) -> tuple[tuple[str, str], ...]:
    if design_build_result is None:
        return ()
    return tuple(
        (level, column)
        for level in design_build_result.block_levels
        for column in (design_build_result.block_columns.get(level),)
        if column is not None
    )


def _block_column_names(
    design_build_result: DesignMatrixBuildResult | None,
) -> tuple[str, ...]:
    return tuple(column for _, column in _block_columns(design_build_result))


def _condition_coverage_rule(paired_design_policy: str) -> str:
    if paired_design_policy == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION:
        return _DIFFERENTIAL_DUPLICATE_CORRELATION_CONDITION_COVERAGE_RULE
    if paired_design_policy == "fixed_block":
        return _DIFFERENTIAL_FIXED_BLOCK_CONDITION_COVERAGE_RULE
    return _DIFFERENTIAL_REJECT_BLOCK_CONDITION_COVERAGE_RULE


def _design_limitations(paired_design_policy: str) -> tuple[str, ...]:
    if paired_design_policy == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION:
        return _DIFFERENTIAL_DUPLICATE_CORRELATION_LIMITATIONS
    if paired_design_policy == "fixed_block":
        return _DIFFERENTIAL_FIXED_BLOCK_LIMITATIONS
    return _DIFFERENTIAL_UNPAIRED_LIMITATIONS


def _unsupported_design_features(paired_design_policy: str) -> tuple[str, ...]:
    if paired_design_policy == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION:
        return _DIFFERENTIAL_DUPLICATE_CORRELATION_UNSUPPORTED_DESIGN_FEATURES
    return _DIFFERENTIAL_UNSUPPORTED_DESIGN_FEATURES


def _condition_columns(
    *,
    coefficient_labels: tuple[str, ...],
    design: ExperimentalDesign,
) -> tuple[str, ...]:
    condition_labels = set(design.condition_labels())
    return tuple(label for label in coefficient_labels if label in condition_labels)


def _fixed_effect_covariate_provenance(
    request: ValidatedDifferentialAnalysisRequest,
) -> tuple[DifferentialFixedEffectCovariateProvenance, ...]:
    design_build_result = request.design_build_result
    if design_build_result is None:
        return ()

    records: list[DifferentialFixedEffectCovariateProvenance] = []
    for covariate in request.design.fixed_effects:
        if not covariate.include_in_model:
            continue
        records.append(
            DifferentialFixedEffectCovariateProvenance(
                name=covariate.name,
                kind=covariate.kind,
                columns=tuple(
                    design_build_result.covariate_columns.get(covariate.name, ())
                ),
                levels=tuple(
                    design_build_result.categorical_levels.get(covariate.name, ())
                ),
                reference_level=design_build_result.reference_levels.get(
                    covariate.name
                ),
                unused_levels=tuple(
                    design_build_result.unused_levels.get(covariate.name, ())
                ),
            )
        )
    return tuple(records)


def _design_description(
    *,
    formula: str,
    covariates: tuple[DifferentialFixedEffectCovariateProvenance, ...],
) -> str:
    if not covariates and "+ block" not in formula:
        return "condition-only fixed-effect design"
    return f"fixed-effect design: {formula}"


def _contrast_description(
    *,
    numerator_condition: str,
    denominator_condition: str,
) -> str:
    return (
        f"condition contrast {numerator_condition} - {denominator_condition}; "
        "non-condition coefficients fixed at 0"
    )


def _condition_replicate_counts(
    design: ExperimentalDesign,
) -> tuple[tuple[str, int], ...]:
    condition_order = design.condition_labels()
    records_by_condition: dict[str, list[str | None]] = {
        condition: [] for condition in condition_order
    }
    for record in design.samples:
        records_by_condition[record.condition].append(record.biological_replicate_id)
    counts: list[tuple[str, int]] = []
    for condition in condition_order:
        biological_ids = records_by_condition[condition]
        if biological_ids and all(value is not None for value in biological_ids):
            count = len({str(value) for value in biological_ids if value is not None})
        else:
            count = len(biological_ids)
        counts.append((condition, int(count)))
    return tuple(counts)


def _technical_replicate_groups(
    workflow_provenance: Mapping[str, object] | None,
) -> tuple[DifferentialTechnicalReplicateGroup, ...]:
    if workflow_provenance is None:
        return ()
    groups_raw = workflow_provenance.get("groups")
    if not isinstance(groups_raw, Sequence) or isinstance(
        groups_raw,
        (str, bytes, bytearray),
    ):
        return ()
    groups: list[DifferentialTechnicalReplicateGroup] = []
    for item in groups_raw:
        if not isinstance(item, Mapping):
            continue
        input_sample_ids = item.get("input_sample_ids")
        technical_replicate_ids = item.get("technical_replicate_ids")
        if (
            not isinstance(input_sample_ids, Sequence)
            or isinstance(input_sample_ids, (str, bytes, bytearray))
            or not isinstance(technical_replicate_ids, Sequence)
            or isinstance(technical_replicate_ids, (str, bytes, bytearray))
        ):
            continue
        groups.append(
            DifferentialTechnicalReplicateGroup(
                condition=str(item.get("condition", "")),
                biological_replicate_id=str(item.get("biological_replicate_id", "")),
                output_sample_id=str(item.get("output_sample_id", "")),
                input_sample_ids=tuple(str(value) for value in input_sample_ids),
                technical_replicate_ids=tuple(
                    str(value) for value in technical_replicate_ids
                ),
                n_technical_replicates=int(item.get("n_technical_replicates", 0)),
            )
        )
    return tuple(groups)


class DifferentialWorkflowProvenanceAssembler:
    """Assemble differential execution provenance from eligibility facts."""

    def run(
        self,
        *,
        workflow_provenance: Mapping[str, object] | None,
        input_feature_ids: tuple[str, ...],
        model_fit_feature_ids: tuple[str, ...],
        failed_model_fit_feature_ids: tuple[str, ...],
        multiple_testing_feature_ids: tuple[str, ...],
        imputation_policy_inputs: DifferentialImputationPolicyInputs | None = None,
        feature_eligibility_inputs: DifferentialFeatureEligibilityInputs | None = None,
        duplicate_correlation: DuplicateCorrelationWorkflowProvenance | None = None,
    ) -> Mapping[str, object]:
        payload: dict[str, object] = (
            {} if workflow_provenance is None else dict(workflow_provenance)
        )
        input_count = int(len(input_feature_ids))
        model_fit_count = int(len(model_fit_feature_ids))
        failed_count = int(len(failed_model_fit_feature_ids))
        multiple_testing_count = int(len(multiple_testing_feature_ids))
        payload["row_attrition_metrics"] = {
            "input_sites": input_count,
            "sites_retained_for_model_fitting": model_fit_count,
            "sites_excluded_before_testing": input_count - model_fit_count,
            "sites_with_failed_model_fit": failed_count,
            "sites_included_in_multiple_testing_family": multiple_testing_count,
        }

        records: list[RowAttritionRecord] = []
        if input_count > model_fit_count:
            model_fit_feature_id_set = set(model_fit_feature_ids)
            records.append(
                RowAttritionRecord(
                    stage="differential_feature_eligibility",
                    input_rows=input_count,
                    output_rows=model_fit_count,
                    removed_rows=input_count - model_fit_count,
                    reason="sites_excluded_before_testing",
                    examples=_row_examples(
                        tuple(
                            feature_id
                            for feature_id in input_feature_ids
                            if feature_id not in model_fit_feature_id_set
                        )
                    ),
                )
            )
        if failed_count:
            records.append(
                RowAttritionRecord(
                    stage="differential_model_fit",
                    input_rows=model_fit_count,
                    output_rows=multiple_testing_count,
                    removed_rows=failed_count,
                    reason="failed_model_fit",
                    examples=_row_examples(failed_model_fit_feature_ids),
                )
            )
        if records:
            payload["row_attrition"] = RowAttritionReport.from_records(
                records
            ).to_payload()
        if imputation_policy_inputs is not None:
            summary = summarize_differential_imputation_inference(
                imputation_policy_inputs=imputation_policy_inputs,
                feature_eligibility_inputs=feature_eligibility_inputs,
            )
            payload["imputation_inference"] = imputation_inference_summary_payload(
                summary
            )
        if duplicate_correlation is not None:
            payload["duplicate_correlation"] = duplicate_correlation.to_payload()
        return payload


def _row_examples(values: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(str(value) for value in values[:5])


__all__ = [
    "DifferentialWorkflowProvenanceAssembler",
    "build_differential_policy_provenance",
    "build_duplicate_correlation_workflow_provenance",
    "finalize_differential_policy_provenance",
]
