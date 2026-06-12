"""Structured differential-policy provenance construction."""

from __future__ import annotations

from collections.abc import Mapping

from phospy.science.design.matrix_builder import describe_fixed_effect_design
from phospy.science.design.models import ExperimentalDesign
from phospy.science.differential.models import (
    DifferentialContrastDefinition,
    DifferentialDesignMatrixSummary,
    DifferentialEmpiricalBayesProvenance,
    DifferentialMissingValuePolicyProvenance,
    DifferentialPolicyProvenance,
    DifferentialReplicatePolicyProvenance,
    DifferentialStatisticalTestingProvenance,
    DifferentialTechnicalReplicateGroup,
    DifferentialUnsupportedDesignPolicyProvenance,
)
from phospy.workflows.differential.models import ValidatedDifferentialAnalysisRequest

_DIFFERENTIAL_TEST_STATISTIC = "moderated_t"
_DIFFERENTIAL_P_VALUE_METHOD = "two_sided_t_distribution_survival_function"
_DIFFERENTIAL_MISSING_VALUE_POLICY = (
    "reject_missing_values_before_differential_execution"
)
_DIFFERENTIAL_MISSING_VALUE_STAGE = "analysis_ready_dataset_boundary"
_DIFFERENTIAL_UNSUPPORTED_DESIGN_FEATURES: tuple[str, ...] = (
    "blocking/paired/repeated-measure differential modelling",
)
_DIFFERENTIAL_UNSUPPORTED_ENFORCEMENT_STAGE = (
    "validation.workflows.differential.ExperimentalDesignContractValidator"
)


def build_differential_policy_provenance(
    *,
    request: ValidatedDifferentialAnalysisRequest,
    design_rank: int,
    residual_degrees_of_freedom: float,
) -> DifferentialPolicyProvenance:
    """Build deterministic structured differential-policy provenance records."""

    design_frame = request.design_matrix.frame
    contrast_frame = request.contrast_matrix.frame
    sample_labels = tuple(str(label) for label in design_frame.index)
    coefficient_labels = tuple(str(label) for label in design_frame.columns)

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
            )
        )

    return DifferentialPolicyProvenance(
        design=DifferentialDesignMatrixSummary(
            formula=describe_fixed_effect_design(request.design),
            sample_labels=sample_labels,
            coefficient_labels=coefficient_labels,
            sample_count=len(sample_labels),
            coefficient_count=len(coefficient_labels),
            rank=int(design_rank),
            residual_degrees_of_freedom=float(residual_degrees_of_freedom),
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
        ),
        missing_values=DifferentialMissingValuePolicyProvenance(
            policy=_DIFFERENTIAL_MISSING_VALUE_POLICY,
            stage=_DIFFERENTIAL_MISSING_VALUE_STAGE,
        ),
        unsupported_design=DifferentialUnsupportedDesignPolicyProvenance(
            intentionally_rejected_features=_DIFFERENTIAL_UNSUPPORTED_DESIGN_FEATURES,
            enforcement_stage=_DIFFERENTIAL_UNSUPPORTED_ENFORCEMENT_STAGE,
        ),
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
    if not isinstance(groups_raw, list):
        return ()
    groups: list[DifferentialTechnicalReplicateGroup] = []
    for item in groups_raw:
        if not isinstance(item, Mapping):
            continue
        input_sample_ids = item.get("input_sample_ids")
        technical_replicate_ids = item.get("technical_replicate_ids")
        if not isinstance(input_sample_ids, list) or not isinstance(
            technical_replicate_ids, list
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
