"""Policy provenance models for differential analysis."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import cast

from phospy.errors.input import PhosPyInputError
from phospy.science.differential.models.duplicate_correlation import (
    DuplicateCorrelationWorkflowProvenance,
)


@dataclass(frozen=True, slots=True)
class DifferentialFixedEffectCovariateProvenance:
    """Resolved fixed-effect covariate columns included in the fitted design."""

    name: str
    kind: str
    columns: tuple[str, ...]
    levels: tuple[str, ...] = ()
    reference_level: str | None = None
    unused_levels: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.name:
            raise PhosPyInputError(
                "differential_policy_provenance.design.covariates[].name must be "
                "non-empty"
            )
        if not self.kind:
            raise PhosPyInputError(
                "differential_policy_provenance.design.covariates[].kind must be "
                "non-empty"
            )
        if not self.columns:
            raise PhosPyInputError(
                "differential_policy_provenance.design.covariates[].columns must be "
                "non-empty"
            )
        object.__setattr__(self, "name", str(self.name))
        object.__setattr__(self, "kind", str(self.kind))
        object.__setattr__(self, "columns", tuple(str(value) for value in self.columns))
        object.__setattr__(self, "levels", tuple(str(value) for value in self.levels))
        object.__setattr__(
            self,
            "reference_level",
            None if self.reference_level is None else str(self.reference_level),
        )
        object.__setattr__(
            self,
            "unused_levels",
            tuple(str(value) for value in self.unused_levels),
        )


@dataclass(frozen=True, slots=True)
class DifferentialDesignMatrixSummary:
    """Structured summary of the resolved differential design matrix."""

    formula: str
    sample_labels: tuple[str, ...]
    coefficient_labels: tuple[str, ...]
    sample_count: int
    coefficient_count: int
    rank: int
    residual_degrees_of_freedom: float
    decomposition_method: str = "not_recorded"
    solver: str = "not_recorded"
    column_scale_method: str = "not_recorded"
    rank_tolerance_policy: str = "not_recorded"
    rank_tolerance: float = 0.0
    condition_number: float = 0.0
    max_condition_number: float = 0.0
    singular_values: tuple[float, ...] = ()
    description: str = ""
    condition_columns: tuple[str, ...] = ()
    covariates: tuple[DifferentialFixedEffectCovariateProvenance, ...] = ()
    paired_design_policy: str = "reject"
    block_id_field_name: str = "block_id"
    block_count: int = 0
    block_levels: tuple[str, ...] = ()
    block_levels_included: tuple[str, ...] = ()
    block_reference_level: str | None = None
    block_columns: tuple[tuple[str, str], ...] = ()
    block_column_names: tuple[str, ...] = ()
    condition_coverage_rule: str = ""
    limitations: tuple[str, ...] = ()
    rank_validation_status: str = "not_recorded"
    conditioning_validation_status: str = "not_recorded"
    estimability_validation_status: str = "not_recorded"

    def __post_init__(self) -> None:
        if not self.formula:
            raise PhosPyInputError(
                "differential_policy_provenance.design.formula must be non-empty"
            )
        if not self.sample_labels:
            raise PhosPyInputError(
                "differential_policy_provenance.design.sample_labels must be non-empty"
            )
        if not self.coefficient_labels:
            raise PhosPyInputError(
                "differential_policy_provenance.design.coefficient_labels must be "
                "non-empty"
            )
        if self.sample_count < 1:
            raise PhosPyInputError(
                "differential_policy_provenance.design.sample_count must be >= 1"
            )
        if self.coefficient_count < 1:
            raise PhosPyInputError(
                "differential_policy_provenance.design.coefficient_count must be >= 1"
            )
        if self.rank < 1:
            raise PhosPyInputError(
                "differential_policy_provenance.design.rank must be >= 1"
            )
        if self.rank > self.coefficient_count:
            raise PhosPyInputError(
                "differential_policy_provenance.design.rank cannot exceed "
                "coefficient_count"
            )
        if self.residual_degrees_of_freedom <= 0.0:
            raise PhosPyInputError(
                "differential_policy_provenance.design.residual_degrees_of_freedom "
                "must be > 0.0"
            )
        if self.rank_tolerance < 0.0 or not math.isfinite(self.rank_tolerance):
            raise PhosPyInputError(
                "differential_policy_provenance.design.rank_tolerance must be "
                "finite and >= 0.0"
            )
        if self.condition_number < 0.0 or not math.isfinite(self.condition_number):
            raise PhosPyInputError(
                "differential_policy_provenance.design.condition_number must be "
                "finite and >= 0.0"
            )
        if self.max_condition_number < 0.0 or not math.isfinite(
            self.max_condition_number
        ):
            raise PhosPyInputError(
                "differential_policy_provenance.design.max_condition_number must "
                "be finite and >= 0.0"
            )
        singular_values = tuple(float(value) for value in self.singular_values)
        if any(value < 0.0 or not math.isfinite(value) for value in singular_values):
            raise PhosPyInputError(
                "differential_policy_provenance.design.singular_values must contain "
                "finite values >= 0.0"
            )
        if self.block_count < 0:
            raise PhosPyInputError(
                "differential_policy_provenance.design.block_count must be >= 0"
            )
        if not self.conditioning_validation_status:
            raise PhosPyInputError(
                "differential_policy_provenance.design."
                "conditioning_validation_status must be non-empty"
            )
        if not self.rank_validation_status:
            raise PhosPyInputError(
                "differential_policy_provenance.design.rank_validation_status must "
                "be non-empty"
            )
        if not self.estimability_validation_status:
            raise PhosPyInputError(
                "differential_policy_provenance.design."
                "estimability_validation_status must be non-empty"
            )
        covariates = tuple(self.covariates)
        for covariate in covariates:
            if not isinstance(
                cast(object, covariate),
                DifferentialFixedEffectCovariateProvenance,
            ):
                raise PhosPyInputError(
                    "differential_policy_provenance.design.covariates must contain "
                    "DifferentialFixedEffectCovariateProvenance values"
                )
        object.__setattr__(self, "formula", str(self.formula))
        object.__setattr__(self, "description", str(self.description))
        object.__setattr__(
            self,
            "sample_labels",
            tuple(str(value) for value in self.sample_labels),
        )
        object.__setattr__(
            self,
            "coefficient_labels",
            tuple(str(value) for value in self.coefficient_labels),
        )
        object.__setattr__(
            self,
            "condition_columns",
            tuple(str(value) for value in self.condition_columns),
        )
        object.__setattr__(
            self,
            "decomposition_method",
            str(self.decomposition_method),
        )
        object.__setattr__(self, "solver", str(self.solver))
        object.__setattr__(
            self,
            "column_scale_method",
            str(self.column_scale_method),
        )
        object.__setattr__(
            self,
            "rank_tolerance_policy",
            str(self.rank_tolerance_policy),
        )
        object.__setattr__(self, "rank_tolerance", float(self.rank_tolerance))
        object.__setattr__(self, "condition_number", float(self.condition_number))
        object.__setattr__(
            self,
            "max_condition_number",
            float(self.max_condition_number),
        )
        object.__setattr__(self, "singular_values", singular_values)
        object.__setattr__(self, "covariates", covariates)
        object.__setattr__(
            self,
            "paired_design_policy",
            str(self.paired_design_policy),
        )
        object.__setattr__(
            self,
            "block_id_field_name",
            str(self.block_id_field_name),
        )
        object.__setattr__(self, "block_count", int(self.block_count))
        object.__setattr__(
            self,
            "block_levels",
            tuple(str(value) for value in self.block_levels),
        )
        object.__setattr__(
            self,
            "block_levels_included",
            tuple(str(value) for value in self.block_levels_included),
        )
        object.__setattr__(
            self,
            "block_reference_level",
            (
                None
                if self.block_reference_level is None
                else str(self.block_reference_level)
            ),
        )
        object.__setattr__(
            self,
            "block_columns",
            tuple((str(level), str(column)) for level, column in self.block_columns),
        )
        object.__setattr__(
            self,
            "block_column_names",
            tuple(str(value) for value in self.block_column_names),
        )
        object.__setattr__(
            self,
            "condition_coverage_rule",
            str(self.condition_coverage_rule),
        )
        object.__setattr__(
            self,
            "limitations",
            tuple(str(value) for value in self.limitations),
        )
        object.__setattr__(
            self,
            "rank_validation_status",
            str(self.rank_validation_status),
        )
        object.__setattr__(
            self,
            "conditioning_validation_status",
            str(self.conditioning_validation_status),
        )
        object.__setattr__(
            self,
            "estimability_validation_status",
            str(self.estimability_validation_status),
        )


@dataclass(frozen=True, slots=True)
class DifferentialContrastDefinition:
    """Structured differential contrast definition."""

    name: str
    numerator_condition: str
    denominator_condition: str
    coefficients: tuple[tuple[str, float], ...]
    description: str = ""

    def __post_init__(self) -> None:
        if not self.name:
            raise PhosPyInputError(
                "differential_policy_provenance.contrasts[].name must be non-empty"
            )
        if not self.numerator_condition:
            raise PhosPyInputError(
                "differential_policy_provenance.contrasts[].numerator_condition must "
                "be non-empty"
            )
        if not self.denominator_condition:
            raise PhosPyInputError(
                "differential_policy_provenance.contrasts[].denominator_condition "
                "must be non-empty"
            )
        if not self.coefficients:
            raise PhosPyInputError(
                "differential_policy_provenance.contrasts[].coefficients must be "
                "non-empty"
            )
        object.__setattr__(self, "name", str(self.name))
        object.__setattr__(
            self,
            "numerator_condition",
            str(self.numerator_condition),
        )
        object.__setattr__(
            self,
            "denominator_condition",
            str(self.denominator_condition),
        )
        object.__setattr__(
            self,
            "coefficients",
            tuple((str(name), float(value)) for name, value in self.coefficients),
        )
        object.__setattr__(self, "description", str(self.description))


@dataclass(frozen=True, slots=True)
class DifferentialTechnicalReplicateGroup:
    """Structured technical-replicate lineage for one resolved group."""

    condition: str
    biological_replicate_id: str
    output_sample_id: str
    input_sample_ids: tuple[str, ...]
    technical_replicate_ids: tuple[str, ...]
    n_technical_replicates: int

    def __post_init__(self) -> None:
        if not self.condition:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates.technical_replicate_groups[]"
                ".condition must be non-empty"
            )
        if not self.biological_replicate_id:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates.technical_replicate_groups[]"
                ".biological_replicate_id must be non-empty"
            )
        if not self.output_sample_id:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates.technical_replicate_groups[]"
                ".output_sample_id must be non-empty"
            )
        if not self.input_sample_ids:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates.technical_replicate_groups[]"
                ".input_sample_ids must be non-empty"
            )
        if self.n_technical_replicates < 1:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates.technical_replicate_groups[]"
                ".n_technical_replicates must be >= 1"
            )


@dataclass(frozen=True, slots=True)
class DifferentialReplicatePolicyProvenance:
    """Structured replicate/group requirements for differential analysis."""

    minimum_condition_replicates: int
    reliability_profile: str
    technical_replicate_policy: str
    condition_replicate_counts: tuple[tuple[str, int], ...]
    technical_replicate_groups: tuple[DifferentialTechnicalReplicateGroup, ...] = ()

    def __post_init__(self) -> None:
        if self.minimum_condition_replicates < 1:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates."
                "minimum_condition_replicates must be >= 1"
            )
        if not self.reliability_profile:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates."
                "reliability_profile must be non-empty"
            )
        if not self.technical_replicate_policy:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates."
                "technical_replicate_policy must be non-empty"
            )
        if not self.condition_replicate_counts:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates.condition_replicate_counts "
                "must be non-empty"
            )
        object.__setattr__(
            self,
            "reliability_profile",
            str(self.reliability_profile),
        )


@dataclass(frozen=True, slots=True)
class DifferentialEmpiricalBayesProvenance:
    """Structured empirical-Bayes moderation settings."""

    method: str
    robust: bool
    trend: bool
    winsor_tail_p: tuple[float, float]

    def __post_init__(self) -> None:
        if not self.method:
            raise PhosPyInputError(
                "differential_policy_provenance.empirical_bayes.method must be "
                "non-empty"
            )


@dataclass(frozen=True, slots=True)
class DifferentialStatisticalTestingProvenance:
    """Structured p-value and multiple-testing adjustment settings."""

    test_statistic: str
    p_value_method: str
    adjusted_p_value_method: str
    input_intensity_scale: str = "not_recorded"
    input_intensity_scale_evidence_level: str = "not_recorded"
    input_intensity_scale_source: str = "not_recorded"
    logfc_interpretation: str = "not_recorded"
    allow_suspicious_declared_input_scale: bool = False

    def __post_init__(self) -> None:
        if not self.test_statistic:
            raise PhosPyInputError(
                "differential_policy_provenance.statistical_testing.test_statistic "
                "must be non-empty"
            )
        if not self.p_value_method:
            raise PhosPyInputError(
                "differential_policy_provenance.statistical_testing.p_value_method "
                "must be non-empty"
            )
        if not self.adjusted_p_value_method:
            raise PhosPyInputError(
                "differential_policy_provenance.statistical_testing."
                "adjusted_p_value_method must be non-empty"
            )
        if not self.input_intensity_scale:
            raise PhosPyInputError(
                "differential_policy_provenance.statistical_testing."
                "input_intensity_scale must be non-empty"
            )
        if not self.input_intensity_scale_evidence_level:
            raise PhosPyInputError(
                "differential_policy_provenance.statistical_testing."
                "input_intensity_scale_evidence_level must be non-empty"
            )
        if not self.input_intensity_scale_source:
            raise PhosPyInputError(
                "differential_policy_provenance.statistical_testing."
                "input_intensity_scale_source must be non-empty"
            )
        if not self.logfc_interpretation:
            raise PhosPyInputError(
                "differential_policy_provenance.statistical_testing."
                "logfc_interpretation must be non-empty"
            )
        if not isinstance(
            cast(object, self.allow_suspicious_declared_input_scale), bool
        ):
            raise PhosPyInputError(
                "differential_policy_provenance.statistical_testing."
                "allow_suspicious_declared_input_scale must be a bool"
            )
        object.__setattr__(self, "test_statistic", str(self.test_statistic))
        object.__setattr__(self, "p_value_method", str(self.p_value_method))
        object.__setattr__(
            self,
            "adjusted_p_value_method",
            str(self.adjusted_p_value_method),
        )
        object.__setattr__(
            self,
            "input_intensity_scale",
            str(self.input_intensity_scale),
        )
        object.__setattr__(
            self,
            "input_intensity_scale_evidence_level",
            str(self.input_intensity_scale_evidence_level),
        )
        object.__setattr__(
            self,
            "input_intensity_scale_source",
            str(self.input_intensity_scale_source),
        )
        object.__setattr__(
            self,
            "logfc_interpretation",
            str(self.logfc_interpretation),
        )
        object.__setattr__(
            self,
            "allow_suspicious_declared_input_scale",
            self.allow_suspicious_declared_input_scale,
        )


@dataclass(frozen=True, slots=True)
class DifferentialMissingValuePolicyProvenance:
    """Structured missing-value handling policy for differential execution."""

    policy: str
    stage: str
    imputed_value_policy: str = "reject"
    imputed_value_max_fraction: float = 0.0
    imputation_metadata_required: bool = False
    adjusted_p_value_scope: str = "all_tested_features"
    tested_feature_count: int = 0
    withheld_feature_count: int = 0
    tested_imputed_feature_count: int = 0
    tested_imputed_cell_count: int = 0
    observed_only_fit: bool = False
    residual_df_adjusted_for_imputation: bool = False
    inferential_status: str = "not_applicable"
    adjusted_p_value_denominator_feature_count: int = 0
    limitations: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.policy:
            raise PhosPyInputError(
                "differential_policy_provenance.missing_values.policy must be non-empty"
            )
        if not self.stage:
            raise PhosPyInputError(
                "differential_policy_provenance.missing_values.stage must be non-empty"
            )
        if not self.imputed_value_policy:
            raise PhosPyInputError(
                "differential_policy_provenance.missing_values."
                "imputed_value_policy must be non-empty"
            )
        object.__setattr__(self, "policy", str(self.policy))
        object.__setattr__(self, "stage", str(self.stage))
        object.__setattr__(
            self,
            "imputed_value_policy",
            str(self.imputed_value_policy),
        )
        object.__setattr__(
            self,
            "imputed_value_max_fraction",
            float(self.imputed_value_max_fraction),
        )
        object.__setattr__(
            self,
            "imputation_metadata_required",
            bool(self.imputation_metadata_required),
        )
        object.__setattr__(
            self,
            "adjusted_p_value_scope",
            str(self.adjusted_p_value_scope),
        )
        object.__setattr__(
            self,
            "tested_feature_count",
            _require_non_negative_int(
                self.tested_feature_count,
                field_name=(
                    "differential_policy_provenance.missing_values.tested_feature_count"
                ),
            ),
        )
        object.__setattr__(
            self,
            "withheld_feature_count",
            _require_non_negative_int(
                self.withheld_feature_count,
                field_name=(
                    "differential_policy_provenance.missing_values."
                    "withheld_feature_count"
                ),
            ),
        )
        object.__setattr__(
            self,
            "tested_imputed_feature_count",
            _require_non_negative_int(
                self.tested_imputed_feature_count,
                field_name=(
                    "differential_policy_provenance.missing_values."
                    "tested_imputed_feature_count"
                ),
            ),
        )
        object.__setattr__(
            self,
            "tested_imputed_cell_count",
            _require_non_negative_int(
                self.tested_imputed_cell_count,
                field_name=(
                    "differential_policy_provenance.missing_values."
                    "tested_imputed_cell_count"
                ),
            ),
        )
        if not isinstance(cast(object, self.observed_only_fit), bool):
            raise PhosPyInputError(
                "differential_policy_provenance.missing_values.observed_only_fit "
                "must be a bool"
            )
        if not isinstance(
            cast(object, self.residual_df_adjusted_for_imputation),
            bool,
        ):
            raise PhosPyInputError(
                "differential_policy_provenance.missing_values."
                "residual_df_adjusted_for_imputation must be a bool"
            )
        object.__setattr__(self, "observed_only_fit", self.observed_only_fit)
        object.__setattr__(
            self,
            "residual_df_adjusted_for_imputation",
            self.residual_df_adjusted_for_imputation,
        )
        if not self.inferential_status:
            raise PhosPyInputError(
                "differential_policy_provenance.missing_values."
                "inferential_status must be non-empty"
            )
        object.__setattr__(
            self,
            "inferential_status",
            str(self.inferential_status),
        )
        object.__setattr__(
            self,
            "adjusted_p_value_denominator_feature_count",
            _require_non_negative_int(
                self.adjusted_p_value_denominator_feature_count,
                field_name=(
                    "differential_policy_provenance.missing_values."
                    "adjusted_p_value_denominator_feature_count"
                ),
            ),
        )
        object.__setattr__(
            self,
            "limitations",
            tuple(str(value) for value in self.limitations),
        )


@dataclass(frozen=True, slots=True)
class DifferentialUnsupportedDesignPolicyProvenance:
    """Structured record of unsupported differential-design features."""

    intentionally_rejected_features: tuple[str, ...]
    enforcement_stage: str
    policy: str = "reject_unsupported_design_features_before_execution"

    def __post_init__(self) -> None:
        if not self.intentionally_rejected_features:
            raise PhosPyInputError(
                "differential_policy_provenance.unsupported_design."
                "intentionally_rejected_features must be non-empty"
            )
        if not self.enforcement_stage:
            raise PhosPyInputError(
                "differential_policy_provenance.unsupported_design.enforcement_stage "
                "must be non-empty"
            )
        if not self.policy:
            raise PhosPyInputError(
                "differential_policy_provenance.unsupported_design.policy must be "
                "non-empty"
            )
        object.__setattr__(
            self,
            "intentionally_rejected_features",
            tuple(str(value) for value in self.intentionally_rejected_features),
        )
        object.__setattr__(self, "enforcement_stage", str(self.enforcement_stage))
        object.__setattr__(self, "policy", str(self.policy))


@dataclass(frozen=True, slots=True)
class DifferentialPolicyProvenance:
    """Structured differential-analysis statistical policy provenance."""

    design: DifferentialDesignMatrixSummary
    contrasts: tuple[DifferentialContrastDefinition, ...]
    replicates: DifferentialReplicatePolicyProvenance
    empirical_bayes: DifferentialEmpiricalBayesProvenance
    statistical_testing: DifferentialStatisticalTestingProvenance
    missing_values: DifferentialMissingValuePolicyProvenance
    unsupported_design: DifferentialUnsupportedDesignPolicyProvenance
    duplicate_correlation: DuplicateCorrelationWorkflowProvenance | None = None

    def __post_init__(self) -> None:
        if not self.contrasts:
            raise PhosPyInputError(
                "differential_policy_provenance.contrasts must be non-empty"
            )
        duplicate_correlation = self.duplicate_correlation
        if duplicate_correlation is not None and not isinstance(
            cast(object, duplicate_correlation),
            DuplicateCorrelationWorkflowProvenance,
        ):
            raise PhosPyInputError(
                "differential_policy_provenance.duplicate_correlation must be "
                "DuplicateCorrelationWorkflowProvenance or None"
            )
        if self.design.paired_design_policy != "duplicate_correlation":
            if duplicate_correlation is not None:
                raise PhosPyInputError(
                    "differential_policy_provenance.duplicate_correlation is only "
                    "valid when design.paired_design_policy='duplicate_correlation'"
                )
            return
        if self.design.block_columns or self.design.block_column_names:
            raise PhosPyInputError(
                "differential_policy_provenance duplicate_correlation design must "
                "not include fixed block columns"
            )
        if duplicate_correlation is None:
            return
        if duplicate_correlation.sample_count != self.design.sample_count:
            raise PhosPyInputError(
                "differential_policy_provenance.duplicate_correlation.sample_count "
                "must match design.sample_count"
            )
        if duplicate_correlation.block_count != self.design.block_count:
            raise PhosPyInputError(
                "differential_policy_provenance.duplicate_correlation.block_count "
                "must match design.block_count"
            )
        if duplicate_correlation.design_rank != self.design.rank:
            raise PhosPyInputError(
                "differential_policy_provenance.duplicate_correlation.design_rank "
                "must match design.rank"
            )


def _require_non_negative_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(f"{field_name} must be a non-negative integer")
    if value < 0:
        raise PhosPyInputError(f"{field_name} must be >= 0")
    return int(value)


__all__ = [
    "DifferentialContrastDefinition",
    "DifferentialDesignMatrixSummary",
    "DifferentialEmpiricalBayesProvenance",
    "DifferentialFixedEffectCovariateProvenance",
    "DifferentialMissingValuePolicyProvenance",
    "DifferentialPolicyProvenance",
    "DifferentialReplicatePolicyProvenance",
    "DifferentialStatisticalTestingProvenance",
    "DifferentialTechnicalReplicateGroup",
    "DifferentialUnsupportedDesignPolicyProvenance",
    "DuplicateCorrelationWorkflowProvenance",
]
