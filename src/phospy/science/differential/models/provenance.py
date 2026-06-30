"""Policy provenance models for differential analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import cast

from phospy.errors.input import PhosPyInputError


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
        if self.block_count < 0:
            raise PhosPyInputError(
                "differential_policy_provenance.design.block_count must be >= 0"
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
    technical_replicate_policy: str
    condition_replicate_counts: tuple[tuple[str, int], ...]
    technical_replicate_groups: tuple[DifferentialTechnicalReplicateGroup, ...] = ()

    def __post_init__(self) -> None:
        if self.minimum_condition_replicates < 1:
            raise PhosPyInputError(
                "differential_policy_provenance.replicates."
                "minimum_condition_replicates must be >= 1"
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
    logfc_interpretation: str = "not_recorded"

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
        if not self.logfc_interpretation:
            raise PhosPyInputError(
                "differential_policy_provenance.statistical_testing."
                "logfc_interpretation must be non-empty"
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
            "logfc_interpretation",
            str(self.logfc_interpretation),
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

    def __post_init__(self) -> None:
        if not self.contrasts:
            raise PhosPyInputError(
                "differential_policy_provenance.contrasts must be non-empty"
            )


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
]
