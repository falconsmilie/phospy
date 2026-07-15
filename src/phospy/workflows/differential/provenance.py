"""Structured differential-policy provenance construction."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from phospy.contracts.configs.differential import (
    IMPUTED_VALUE_POLICY_REJECT,
    IMPUTED_VALUE_POLICY_WITHHOLD_IMPUTED_FEATURES,
)
from phospy.science.design.matrix_builder import (
    DesignMatrixBuildResult,
    describe_fixed_effect_design,
)
from phospy.science.design.models import ExperimentalDesign
from phospy.science.differential.linear_model import DifferentialDesignDecomposition
from phospy.science.differential.models import (
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
from phospy.workflows.differential.models import ValidatedDifferentialAnalysisRequest
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
    "correlated repeated-measure differential modelling beyond explicit fixed blocks",
    "duplicateCorrelation-style correlated-replicate modelling",
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
_DIFFERENTIAL_UNPAIRED_LIMITATIONS: tuple[str, ...] = (
    "paired_design_policy='reject' does not construct fixed-block terms",
    (
        "explicit block_id metadata is rejected unless "
        "paired_design_policy='fixed_block'"
    ),
    (
        "unpaired condition and covariate workflows do not fit "
        "duplicateCorrelation, mixed-effects, or random subject-effect models"
    ),
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


def build_differential_policy_provenance(
    *,
    request: ValidatedDifferentialAnalysisRequest,
    design_decomposition: DifferentialDesignDecomposition,
) -> DifferentialPolicyProvenance:
    """Build deterministic structured differential-policy provenance records."""

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
            block_count=_block_count(request.design_build_result),
            block_levels=_block_levels(request.design_build_result),
            block_levels_included=_block_levels(request.design_build_result),
            block_reference_level=_block_reference_level(request.design_build_result),
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
            minimum_condition_replicates=request.config.minimum_condition_replicates,
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
            limitations=_imputation_policy_limitations(
                request.config.imputed_value_policy
            ),
        ),
        unsupported_design=DifferentialUnsupportedDesignPolicyProvenance(
            intentionally_rejected_features=_DIFFERENTIAL_UNSUPPORTED_DESIGN_FEATURES,
            enforcement_stage=_DIFFERENTIAL_UNSUPPORTED_ENFORCEMENT_STAGE,
        ),
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


def _block_levels(
    design_build_result: DesignMatrixBuildResult | None,
) -> tuple[str, ...]:
    if design_build_result is None:
        return ()
    return design_build_result.block_levels


def _block_count(
    design_build_result: DesignMatrixBuildResult | None,
) -> int:
    return len(_block_levels(design_build_result))


def _block_reference_level(
    design_build_result: DesignMatrixBuildResult | None,
) -> str | None:
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
    if paired_design_policy == "fixed_block":
        return _DIFFERENTIAL_FIXED_BLOCK_CONDITION_COVERAGE_RULE
    return _DIFFERENTIAL_REJECT_BLOCK_CONDITION_COVERAGE_RULE


def _design_limitations(paired_design_policy: str) -> tuple[str, ...]:
    if paired_design_policy == "fixed_block":
        return _DIFFERENTIAL_FIXED_BLOCK_LIMITATIONS
    return _DIFFERENTIAL_UNPAIRED_LIMITATIONS


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


__all__ = ["build_differential_policy_provenance"]
