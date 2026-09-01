"""Policy provenance models for differential analysis."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, cast

from phospy.errors.input import PhosPyInputError
from phospy.provenance.immutability import freeze_json_mapping
from phospy.provenance.models import TableFingerprint
from phospy.science.differential.models.duplicate_correlation import (
    DuplicateCorrelationWorkflowProvenance,
)

DIFFERENTIAL_PROTEIN_AWARE_POLICY_METHOD_VERSION = "1"
DIFFERENTIAL_PROTEIN_AWARE_MODEL_FORMULA = "y_s = X beta_s + z_p(s) gamma_s + error"
DIFFERENTIAL_PROTEIN_AWARE_LOGFC_INTERPRETATION = (
    "requested phosphosite condition contrast conditional on the matched "
    "total-protein abundance covariate"
)
DIFFERENTIAL_PROTEIN_AWARE_NUISANCE_COEFFICIENT_NAME = "protein_covariate"
DIFFERENTIAL_PROTEIN_AWARE_CLAIM_STATUS_EXPERIMENTAL = "experimental"
DIFFERENTIAL_PROTEIN_AWARE_CONDITION_NUMBER_SUMMARY_SCOPE = (
    "successfully_fitted_augmented_designs_by_total_protein_row"
)
DIFFERENTIAL_PROTEIN_AWARE_CONDITION_NUMBER_SUMMARY_STATISTIC = (
    "min_median_max_condition_number_across_fitted_augmented_designs"
)
DIFFERENTIAL_PROTEIN_AWARE_PHOSPHO_MATRIX_FINGERPRINT_NAME = (
    "differential.input.phospho_matrix"
)
DIFFERENTIAL_PROTEIN_AWARE_MATCHED_PAIRS_FINGERPRINT_NAME = (
    "differential.input.protein_matched_pairs"
)
DIFFERENTIAL_PROTEIN_AWARE_COVARIATE_MATRIX_FINGERPRINT_NAME = (
    "differential.input.protein_covariate_matrix"
)
DIFFERENTIAL_PROTEIN_AWARE_SITE_ELIGIBILITY_FINGERPRINT_NAME = (
    "differential.input.protein_site_eligibility"
)
DIFFERENTIAL_PROTEIN_AWARE_DESIGN_MATRIX_FINGERPRINT_NAME = (
    "differential.input.design_matrix"
)
DIFFERENTIAL_PROTEIN_AWARE_CONTRAST_MATRIX_FINGERPRINT_NAME = (
    "differential.input.contrast_matrix"
)
_DIFFERENTIAL_PROTEIN_AWARE_LIMITATIONS: tuple[str, ...] = (
    "protein-aware differential analysis is experimental",
    (
        "condition effects are estimated conditional on the matched measured "
        "total-protein covariate and are not causal separation of abundance and "
        "phosphorylation regulation"
    ),
    (
        "the estimator does not claim MSstatsPTM parity and does not fit a joint "
        "phosphosite-total-protein model"
    ),
    "total-protein covariates are mean-centered and not standardized",
    (
        "differential analysis does not automatically normalize or impute "
        "protein covariates"
    ),
    (
        "there is no per-site fallback from the protein-aware lane to ordinary "
        "phosphosite-only differential analysis"
    ),
    "stoichiometry and occupancy are not estimated",
    "duplicate_correlation and mixed-effect protein-aware models are unsupported",
    (
        "upstream phosphosite normalization and protein preprocessing are recorded "
        "as provenance and are not changed by differential analysis"
    ),
)
_DIFFERENTIAL_PROTEIN_AWARE_UNSUPPORTED_CLAIMS: tuple[str, ...] = (
    "MSstatsPTM parity",
    "joint phosphosite-total-protein modelling",
    "causal separation of protein abundance and phosphorylation regulation",
    "stoichiometry or occupancy estimation",
    "automatic protein normalization",
    "protein covariate imputation",
    "duplicate_correlation or mixed-effect protein-aware inference",
    "fallback to ordinary phosphosite-only differential analysis",
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
class DifferentialProteinAwareInputFingerprints:
    """Strict fingerprints for exact protein-aware differential inputs."""

    phospho_matrix: TableFingerprint
    protein_matched_pairs: TableFingerprint
    protein_covariate_matrix: TableFingerprint
    protein_site_eligibility: TableFingerprint
    design_matrix: TableFingerprint
    contrast_matrix: TableFingerprint

    def __post_init__(self) -> None:
        _require_fingerprint(
            self.phospho_matrix,
            field_name=(
                "differential_policy_provenance.protein_aware.input_fingerprints."
                "phospho_matrix"
            ),
            expected_name=DIFFERENTIAL_PROTEIN_AWARE_PHOSPHO_MATRIX_FINGERPRINT_NAME,
        )
        _require_fingerprint(
            self.protein_matched_pairs,
            field_name=(
                "differential_policy_provenance.protein_aware.input_fingerprints."
                "protein_matched_pairs"
            ),
            expected_name=DIFFERENTIAL_PROTEIN_AWARE_MATCHED_PAIRS_FINGERPRINT_NAME,
        )
        _require_fingerprint(
            self.protein_covariate_matrix,
            field_name=(
                "differential_policy_provenance.protein_aware.input_fingerprints."
                "protein_covariate_matrix"
            ),
            expected_name=DIFFERENTIAL_PROTEIN_AWARE_COVARIATE_MATRIX_FINGERPRINT_NAME,
        )
        _require_fingerprint(
            self.protein_site_eligibility,
            field_name=(
                "differential_policy_provenance.protein_aware.input_fingerprints."
                "protein_site_eligibility"
            ),
            expected_name=DIFFERENTIAL_PROTEIN_AWARE_SITE_ELIGIBILITY_FINGERPRINT_NAME,
        )
        _require_fingerprint(
            self.design_matrix,
            field_name=(
                "differential_policy_provenance.protein_aware.input_fingerprints."
                "design_matrix"
            ),
            expected_name=DIFFERENTIAL_PROTEIN_AWARE_DESIGN_MATRIX_FINGERPRINT_NAME,
        )
        _require_fingerprint(
            self.contrast_matrix,
            field_name=(
                "differential_policy_provenance.protein_aware.input_fingerprints."
                "contrast_matrix"
            ),
            expected_name=DIFFERENTIAL_PROTEIN_AWARE_CONTRAST_MATRIX_FINGERPRINT_NAME,
        )


@dataclass(frozen=True, slots=True)
class DifferentialProteinAwarePolicyProvenance:
    """Typed statistical-policy provenance for the protein-aware estimator."""

    method_id: str
    method_version: str
    claim_status: str
    model_formula: str
    logfc_interpretation: str
    nuisance_coefficient_name: str
    preparation_schema_version: int
    preparation_policy: str
    protein_mapping_policy: str
    protein_mapping_policy_parameters: Mapping[str, object]
    protein_reference_context: Mapping[str, object]
    phosphosite_transformation_state: Mapping[str, object]
    total_protein_transformation_state: Mapping[str, object]
    prior_total_protein_correction_state: Mapping[str, object]
    phosphosite_normalisation_state: str
    protein_covariate_centered: bool
    protein_covariate_standardized: bool
    automatic_protein_normalization: bool
    protein_imputation: bool
    phosphosite_only_fallback: bool
    execution_sample_order: tuple[str, ...]
    design_subset_behavior: str
    technical_aggregation_policy: str
    duplicate_correlation_policy: str
    total_site_count: int
    ordinary_eligible_site_count: int
    protein_preparation_eligible_site_count: int
    distinct_matched_protein_row_count: int
    fitted_protein_row_count: int
    tested_site_count: int
    withheld_site_count: int
    status_counts: tuple[tuple[str, int], ...]
    reason_counts: tuple[tuple[str, int], ...]
    base_design_rank: int
    base_residual_degrees_of_freedom: float
    expected_augmented_rank: int
    common_augmented_rank: int
    common_augmented_residual_degrees_of_freedom: float
    condition_number_summary_scope: str
    condition_number_summary_statistic: str
    min_augmented_condition_number: float | None
    median_augmented_condition_number: float | None
    max_augmented_condition_number: float | None
    input_fingerprints: DifferentialProteinAwareInputFingerprints
    limitations: tuple[str, ...] = _DIFFERENTIAL_PROTEIN_AWARE_LIMITATIONS
    unsupported_claims: tuple[str, ...] = _DIFFERENTIAL_PROTEIN_AWARE_UNSUPPORTED_CLAIMS

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "method_id",
            _require_non_empty_text(
                self.method_id,
                field_name="differential_policy_provenance.protein_aware.method_id",
            ),
        )
        object.__setattr__(
            self,
            "method_version",
            _require_non_empty_text(
                self.method_version,
                field_name=(
                    "differential_policy_provenance.protein_aware.method_version"
                ),
            ),
        )
        claim_status = _require_non_empty_text(
            self.claim_status,
            field_name="differential_policy_provenance.protein_aware.claim_status",
        )
        if claim_status != DIFFERENTIAL_PROTEIN_AWARE_CLAIM_STATUS_EXPERIMENTAL:
            raise PhosPyInputError(
                "differential_policy_provenance.protein_aware.claim_status must be "
                f"{DIFFERENTIAL_PROTEIN_AWARE_CLAIM_STATUS_EXPERIMENTAL!r}"
            )
        object.__setattr__(self, "claim_status", claim_status)
        model_formula = _require_non_empty_text(
            self.model_formula,
            field_name="differential_policy_provenance.protein_aware.model_formula",
        )
        if model_formula != DIFFERENTIAL_PROTEIN_AWARE_MODEL_FORMULA:
            raise PhosPyInputError(
                "differential_policy_provenance.protein_aware.model_formula must "
                "match ADR-0049"
            )
        object.__setattr__(self, "model_formula", model_formula)
        object.__setattr__(
            self,
            "logfc_interpretation",
            _require_non_empty_text(
                self.logfc_interpretation,
                field_name=(
                    "differential_policy_provenance.protein_aware.logfc_interpretation"
                ),
            ),
        )
        nuisance_coefficient_name = _require_non_empty_text(
            self.nuisance_coefficient_name,
            field_name=(
                "differential_policy_provenance.protein_aware.nuisance_coefficient_name"
            ),
        )
        if (
            nuisance_coefficient_name
            != DIFFERENTIAL_PROTEIN_AWARE_NUISANCE_COEFFICIENT_NAME
        ):
            raise PhosPyInputError(
                "differential_policy_provenance.protein_aware."
                "nuisance_coefficient_name must be "
                f"{DIFFERENTIAL_PROTEIN_AWARE_NUISANCE_COEFFICIENT_NAME!r}"
            )
        object.__setattr__(
            self,
            "nuisance_coefficient_name",
            nuisance_coefficient_name,
        )
        object.__setattr__(
            self,
            "preparation_schema_version",
            _require_positive_int(
                self.preparation_schema_version,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "preparation_schema_version"
                ),
            ),
        )
        object.__setattr__(
            self,
            "preparation_policy",
            _require_non_empty_text(
                self.preparation_policy,
                field_name=(
                    "differential_policy_provenance.protein_aware.preparation_policy"
                ),
            ),
        )
        object.__setattr__(
            self,
            "protein_mapping_policy",
            _require_non_empty_text(
                self.protein_mapping_policy,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "protein_mapping_policy"
                ),
            ),
        )
        object.__setattr__(
            self,
            "protein_mapping_policy_parameters",
            freeze_json_mapping(
                self.protein_mapping_policy_parameters,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "protein_mapping_policy_parameters"
                ),
            ),
        )
        object.__setattr__(
            self,
            "protein_reference_context",
            freeze_json_mapping(
                self.protein_reference_context,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "protein_reference_context"
                ),
            ),
        )
        object.__setattr__(
            self,
            "phosphosite_transformation_state",
            freeze_json_mapping(
                self.phosphosite_transformation_state,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "phosphosite_transformation_state"
                ),
            ),
        )
        object.__setattr__(
            self,
            "total_protein_transformation_state",
            freeze_json_mapping(
                self.total_protein_transformation_state,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "total_protein_transformation_state"
                ),
            ),
        )
        object.__setattr__(
            self,
            "prior_total_protein_correction_state",
            freeze_json_mapping(
                self.prior_total_protein_correction_state,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "prior_total_protein_correction_state"
                ),
            ),
        )
        object.__setattr__(
            self,
            "phosphosite_normalisation_state",
            _require_non_empty_text(
                self.phosphosite_normalisation_state,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "phosphosite_normalisation_state"
                ),
            ),
        )
        _require_exact_bool(
            self.protein_covariate_centered,
            expected=True,
            field_name=(
                "differential_policy_provenance.protein_aware."
                "protein_covariate_centered"
            ),
        )
        _require_exact_bool(
            self.protein_covariate_standardized,
            expected=False,
            field_name=(
                "differential_policy_provenance.protein_aware."
                "protein_covariate_standardized"
            ),
        )
        _require_exact_bool(
            self.automatic_protein_normalization,
            expected=False,
            field_name=(
                "differential_policy_provenance.protein_aware."
                "automatic_protein_normalization"
            ),
        )
        _require_exact_bool(
            self.protein_imputation,
            expected=False,
            field_name=(
                "differential_policy_provenance.protein_aware.protein_imputation"
            ),
        )
        _require_exact_bool(
            self.phosphosite_only_fallback,
            expected=False,
            field_name=(
                "differential_policy_provenance.protein_aware.phosphosite_only_fallback"
            ),
        )
        object.__setattr__(
            self,
            "execution_sample_order",
            _require_non_empty_text_tuple(
                self.execution_sample_order,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "execution_sample_order"
                ),
            ),
        )
        object.__setattr__(
            self,
            "design_subset_behavior",
            _require_non_empty_text(
                self.design_subset_behavior,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "design_subset_behavior"
                ),
            ),
        )
        object.__setattr__(
            self,
            "technical_aggregation_policy",
            _require_non_empty_text(
                self.technical_aggregation_policy,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "technical_aggregation_policy"
                ),
            ),
        )
        object.__setattr__(
            self,
            "duplicate_correlation_policy",
            _require_non_empty_text(
                self.duplicate_correlation_policy,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "duplicate_correlation_policy"
                ),
            ),
        )
        for field_name in (
            "total_site_count",
            "ordinary_eligible_site_count",
            "protein_preparation_eligible_site_count",
            "distinct_matched_protein_row_count",
            "fitted_protein_row_count",
            "tested_site_count",
            "withheld_site_count",
            "base_design_rank",
            "expected_augmented_rank",
            "common_augmented_rank",
        ):
            object.__setattr__(
                self,
                field_name,
                _require_non_negative_int(
                    getattr(self, field_name),
                    field_name=(
                        f"differential_policy_provenance.protein_aware.{field_name}"
                    ),
                ),
            )
        if self.tested_site_count > self.total_site_count:
            raise PhosPyInputError(
                "differential_policy_provenance.protein_aware.tested_site_count "
                "cannot exceed total_site_count"
            )
        if self.withheld_site_count > self.total_site_count:
            raise PhosPyInputError(
                "differential_policy_provenance.protein_aware.withheld_site_count "
                "cannot exceed total_site_count"
            )
        object.__setattr__(
            self,
            "status_counts",
            _require_count_pairs(
                self.status_counts,
                field_name="differential_policy_provenance.protein_aware.status_counts",
            ),
        )
        object.__setattr__(
            self,
            "reason_counts",
            _require_count_pairs(
                self.reason_counts,
                field_name="differential_policy_provenance.protein_aware.reason_counts",
            ),
        )
        _validate_protein_aware_count_consistency(self)
        object.__setattr__(
            self,
            "base_residual_degrees_of_freedom",
            _require_positive_finite_float(
                self.base_residual_degrees_of_freedom,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "base_residual_degrees_of_freedom"
                ),
            ),
        )
        object.__setattr__(
            self,
            "common_augmented_residual_degrees_of_freedom",
            _require_positive_finite_float(
                self.common_augmented_residual_degrees_of_freedom,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "common_augmented_residual_degrees_of_freedom"
                ),
            ),
        )
        object.__setattr__(
            self,
            "condition_number_summary_scope",
            _require_non_empty_text(
                self.condition_number_summary_scope,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "condition_number_summary_scope"
                ),
            ),
        )
        object.__setattr__(
            self,
            "condition_number_summary_statistic",
            _require_non_empty_text(
                self.condition_number_summary_statistic,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "condition_number_summary_statistic"
                ),
            ),
        )
        object.__setattr__(
            self,
            "min_augmented_condition_number",
            _require_optional_non_negative_finite_float(
                self.min_augmented_condition_number,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "min_augmented_condition_number"
                ),
            ),
        )
        object.__setattr__(
            self,
            "median_augmented_condition_number",
            _require_optional_non_negative_finite_float(
                self.median_augmented_condition_number,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "median_augmented_condition_number"
                ),
            ),
        )
        object.__setattr__(
            self,
            "max_augmented_condition_number",
            _require_optional_non_negative_finite_float(
                self.max_augmented_condition_number,
                field_name=(
                    "differential_policy_provenance.protein_aware."
                    "max_augmented_condition_number"
                ),
            ),
        )
        if not isinstance(
            cast(object, self.input_fingerprints),
            DifferentialProteinAwareInputFingerprints,
        ):
            raise PhosPyInputError(
                "differential_policy_provenance.protein_aware.input_fingerprints "
                "must be DifferentialProteinAwareInputFingerprints"
            )
        object.__setattr__(
            self,
            "limitations",
            _require_non_empty_text_tuple(
                self.limitations,
                field_name="differential_policy_provenance.protein_aware.limitations",
            ),
        )
        object.__setattr__(
            self,
            "unsupported_claims",
            _require_non_empty_text_tuple(
                self.unsupported_claims,
                field_name=(
                    "differential_policy_provenance.protein_aware.unsupported_claims"
                ),
            ),
        )


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
    protein_aware: DifferentialProteinAwarePolicyProvenance | None = field(
        default=None,
        repr=False,
    )

    def __post_init__(self) -> None:
        if not self.contrasts:
            raise PhosPyInputError(
                "differential_policy_provenance.contrasts must be non-empty"
            )
        protein_aware = self.protein_aware
        if protein_aware is not None and not isinstance(
            cast(object, protein_aware),
            DifferentialProteinAwarePolicyProvenance,
        ):
            raise PhosPyInputError(
                "differential_policy_provenance.protein_aware must be "
                "DifferentialProteinAwarePolicyProvenance or None"
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
            if protein_aware is not None:
                return
            return
        if protein_aware is not None:
            raise PhosPyInputError(
                "differential_policy_provenance.protein_aware is not valid with "
                "paired_design_policy='duplicate_correlation'"
            )
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


def _require_positive_int(value: object, *, field_name: str) -> int:
    result = _require_non_negative_int(value, field_name=field_name)
    if result < 1:
        raise PhosPyInputError(f"{field_name} must be >= 1")
    return result


def _require_non_empty_text(value: object, *, field_name: str) -> str:
    if value is None:
        raise PhosPyInputError(f"{field_name} must be non-empty")
    text = str(value).strip()
    if not text:
        raise PhosPyInputError(f"{field_name} must be non-empty")
    return text


def _require_non_empty_text_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(
        values,
        Sequence,
    ):
        raise PhosPyInputError(f"{field_name} must be a sequence of strings")
    sequence = cast(Sequence[object], values)
    result = tuple(
        _require_non_empty_text(value, field_name=field_name) for value in sequence
    )
    if not result:
        raise PhosPyInputError(f"{field_name} must be non-empty")
    return result


def _require_count_pairs(
    values: object,
    *,
    field_name: str,
) -> tuple[tuple[str, int], ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(
        values,
        Sequence,
    ):
        raise PhosPyInputError(f"{field_name} must be a sequence of pairs")
    sequence = cast(Sequence[object], values)
    pairs: list[tuple[str, int]] = []
    for position, value in enumerate(sequence):
        if isinstance(value, (str, bytes, bytearray)) or not isinstance(
            value, Sequence
        ):
            raise PhosPyInputError(
                f"{field_name}[{position}] must be a (name, count) pair"
            )
        pair = cast(Sequence[object], value)
        if len(pair) != 2:
            raise PhosPyInputError(
                f"{field_name}[{position}] must be a (name, count) pair"
            )
        key = pair[0]
        count = pair[1]
        pairs.append(
            (
                _require_non_empty_text(
                    key,
                    field_name=f"{field_name}[{position}][0]",
                ),
                _require_non_negative_int(
                    count,
                    field_name=f"{field_name}[{position}][1]",
                ),
            )
        )
    return tuple(pairs)


def _require_exact_bool(value: object, *, expected: bool, field_name: str) -> None:
    if not isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be a bool")
    if bool(value) is not expected:
        raise PhosPyInputError(f"{field_name} must be {expected!r}")


def _require_positive_finite_float(value: object, *, field_name: str) -> float:
    numeric = _require_finite_float(value, field_name=field_name)
    if numeric <= 0.0:
        raise PhosPyInputError(f"{field_name} must be > 0.0")
    return numeric


def _require_optional_non_negative_finite_float(
    value: object,
    *,
    field_name: str,
) -> float | None:
    if value is None:
        return None
    numeric = _require_finite_float(value, field_name=field_name)
    if numeric < 0.0:
        raise PhosPyInputError(f"{field_name} must be >= 0.0")
    return numeric


def _require_finite_float(value: object, *, field_name: str) -> float:
    try:
        numeric = float(cast(Any, value))
    except (TypeError, ValueError) as exc:
        raise PhosPyInputError(f"{field_name} must be a finite number") from exc
    if not math.isfinite(numeric):
        raise PhosPyInputError(f"{field_name} must be finite")
    return numeric


def _require_fingerprint(
    fingerprint: object,
    *,
    field_name: str,
    expected_name: str,
) -> None:
    if not isinstance(fingerprint, TableFingerprint):
        raise PhosPyInputError(f"{field_name} must be a TableFingerprint")
    if fingerprint.name != expected_name:
        raise PhosPyInputError(
            f"{field_name}.name must be {expected_name!r}; got {fingerprint.name!r}"
        )


def _validate_protein_aware_count_consistency(
    provenance: DifferentialProteinAwarePolicyProvenance,
) -> None:
    field_prefix = "differential_policy_provenance.protein_aware"
    if provenance.ordinary_eligible_site_count > provenance.total_site_count:
        raise PhosPyInputError(
            f"{field_prefix}.ordinary_eligible_site_count cannot exceed "
            "total_site_count"
        )
    if (
        provenance.protein_preparation_eligible_site_count
        > provenance.ordinary_eligible_site_count
    ):
        raise PhosPyInputError(
            f"{field_prefix}.protein_preparation_eligible_site_count cannot "
            "exceed ordinary_eligible_site_count"
        )
    if (
        provenance.tested_site_count
        > provenance.protein_preparation_eligible_site_count
    ):
        raise PhosPyInputError(
            f"{field_prefix}.tested_site_count cannot exceed "
            "protein_preparation_eligible_site_count"
        )
    if (
        provenance.fitted_protein_row_count
        > provenance.distinct_matched_protein_row_count
    ):
        raise PhosPyInputError(
            f"{field_prefix}.fitted_protein_row_count cannot exceed "
            "distinct_matched_protein_row_count"
        )
    if (
        provenance.tested_site_count + provenance.withheld_site_count
        != provenance.total_site_count
    ):
        raise PhosPyInputError(
            f"{field_prefix}.tested_site_count plus withheld_site_count must equal "
            "total_site_count"
        )
    status_count_total = sum(count for _, count in provenance.status_counts)
    if status_count_total != provenance.total_site_count:
        raise PhosPyInputError(
            f"{field_prefix}.status_counts must sum to total_site_count"
        )
    reason_count_total = sum(count for _, count in provenance.reason_counts)
    if reason_count_total > provenance.withheld_site_count:
        raise PhosPyInputError(
            f"{field_prefix}.reason_counts cannot exceed withheld_site_count"
        )


__all__ = [
    "DIFFERENTIAL_PROTEIN_AWARE_CLAIM_STATUS_EXPERIMENTAL",
    "DIFFERENTIAL_PROTEIN_AWARE_CONDITION_NUMBER_SUMMARY_SCOPE",
    "DIFFERENTIAL_PROTEIN_AWARE_CONDITION_NUMBER_SUMMARY_STATISTIC",
    "DIFFERENTIAL_PROTEIN_AWARE_CONTRAST_MATRIX_FINGERPRINT_NAME",
    "DIFFERENTIAL_PROTEIN_AWARE_COVARIATE_MATRIX_FINGERPRINT_NAME",
    "DIFFERENTIAL_PROTEIN_AWARE_DESIGN_MATRIX_FINGERPRINT_NAME",
    "DIFFERENTIAL_PROTEIN_AWARE_LOGFC_INTERPRETATION",
    "DIFFERENTIAL_PROTEIN_AWARE_MATCHED_PAIRS_FINGERPRINT_NAME",
    "DIFFERENTIAL_PROTEIN_AWARE_MODEL_FORMULA",
    "DIFFERENTIAL_PROTEIN_AWARE_NUISANCE_COEFFICIENT_NAME",
    "DIFFERENTIAL_PROTEIN_AWARE_PHOSPHO_MATRIX_FINGERPRINT_NAME",
    "DIFFERENTIAL_PROTEIN_AWARE_POLICY_METHOD_VERSION",
    "DIFFERENTIAL_PROTEIN_AWARE_SITE_ELIGIBILITY_FINGERPRINT_NAME",
    "DifferentialContrastDefinition",
    "DifferentialDesignMatrixSummary",
    "DifferentialEmpiricalBayesProvenance",
    "DifferentialFixedEffectCovariateProvenance",
    "DifferentialMissingValuePolicyProvenance",
    "DifferentialPolicyProvenance",
    "DifferentialProteinAwareInputFingerprints",
    "DifferentialProteinAwarePolicyProvenance",
    "DifferentialReplicatePolicyProvenance",
    "DifferentialStatisticalTestingProvenance",
    "DifferentialTechnicalReplicateGroup",
    "DifferentialUnsupportedDesignPolicyProvenance",
    "DuplicateCorrelationWorkflowProvenance",
]
