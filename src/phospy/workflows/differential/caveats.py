"""Structured result caveats for differential workflow execution."""

from __future__ import annotations

import pandas as pd

from phospy.contracts.configs import DifferentialAnalysisConfig
from phospy.contracts.configs.differential import (
    DIFFERENTIAL_PRODUCTION_MINIMUM_CONDITION_REPLICATES,
    DIFFERENTIAL_RELIABILITY_PROFILE_EXPLORATORY_SINGLE_REPLICATE,
    PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
)
from phospy.contracts.result_caveats import ResultCaveat
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DifferentialPolicyProvenance,
)
from phospy.workflows.differential.imputation_inference import (
    DifferentialImputationInferenceSummary,
    imputation_inference_summary_payload,
    summarize_differential_imputation_inference,
)
from phospy.workflows.differential.models import (
    DifferentialFeatureEligibilityInputs,
    DifferentialImputationPolicyInputs,
)
from phospy.workflows.intensity_scale_evidence import (
    build_declared_input_intensity_scale_caveat,
)
from phospy.workflows.result_caveat_helpers import (
    build_direct_trusted_dataset_construction_caveat,
)

DIFFERENTIAL_DIRECT_TRUSTED_DATASET_CAVEAT_CODE = (
    "differential_direct_trusted_dataset_construction"
)
DIFFERENTIAL_DECLARED_SCALE_OVERRIDE_CAVEAT_CODE = (
    "differential_declared_scale_override"
)
DIFFERENTIAL_IMPUTATION_WITHHOLDING_POLICY_CAVEAT_CODE = (
    "differential_imputation_withholding_policy"
)
DIFFERENTIAL_WITHHELD_FEATURES_CAVEAT_CODE = "differential_withheld_features"
DIFFERENTIAL_NARROW_PARITY_ENVELOPE_CAVEAT_CODE = "differential_narrow_parity_envelope"
DIFFERENTIAL_DUPLICATE_CORRELATION_CONSENSUS_CAVEAT_CODE = (
    "differential_duplicate_correlation_consensus"
)
DIFFERENTIAL_EXPLORATORY_SINGLE_REPLICATE_CAVEAT_CODE = (
    "differential_exploratory_single_replicate"
)


def build_differential_result_caveats(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    config: DifferentialAnalysisConfig,
    policy_provenance: DifferentialPolicyProvenance | None,
    imputation_policy_inputs: DifferentialImputationPolicyInputs | None,
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs | None,
    ruv_readiness_enabled: bool,
    ruv_readiness_ready: bool,
) -> tuple[ResultCaveat, ...]:
    """Build compact machine-readable caveats for the public result."""

    caveats: list[ResultCaveat] = []

    declared_input_scale = build_declared_input_intensity_scale_caveat(
        dataset=dataset,
        workflow_scope="differential",
    )
    if declared_input_scale is not None:
        caveats.append(declared_input_scale)

    direct_construction = build_direct_trusted_dataset_construction_caveat(
        dataset=dataset,
        code=DIFFERENTIAL_DIRECT_TRUSTED_DATASET_CAVEAT_CODE,
        workflow_scope="differential",
        workflow_label="differential analysis",
    )
    if direct_construction is not None:
        caveats.append(direct_construction)

    exploratory_details = _exploratory_single_replicate_details(
        config=config,
        policy_provenance=policy_provenance,
    )
    if exploratory_details is not None:
        caveats.append(
            ResultCaveat(
                code=DIFFERENTIAL_EXPLORATORY_SINGLE_REPLICATE_CAVEAT_CODE,
                severity="warning",
                message=(
                    "Differential analysis ran under the explicit exploratory "
                    "single-biological-replicate reliability profile. Results are "
                    "computable model outputs, not production-supported "
                    "inferential evidence."
                ),
                details=exploratory_details,
            )
        )

    declared_scale_details = _declared_scale_override_details(
        dataset=dataset,
        config=config,
    )
    if declared_scale_details is not None:
        caveats.append(
            ResultCaveat(
                code=DIFFERENTIAL_DECLARED_SCALE_OVERRIDE_CAVEAT_CODE,
                severity="warning",
                message=(
                    "Differential analysis accepted a declared input intensity "
                    "scale with diagnostic warnings because "
                    "allow_suspicious_declared_input_scale=True."
                ),
                details=declared_scale_details,
            )
        )

    if imputation_policy_inputs is not None:
        caveats.append(
            _imputation_policy_caveat(
                imputation_policy_inputs=imputation_policy_inputs,
                feature_eligibility_inputs=feature_eligibility_inputs,
            )
        )

    if config.paired_design_policy == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION:
        caveats.append(_duplicate_correlation_consensus_caveat())

    withheld_feature_details = _withheld_feature_details(feature_eligibility_inputs)
    if withheld_feature_details is not None:
        caveats.append(_withheld_feature_caveat(details=withheld_feature_details))

    parity_details = _narrow_parity_envelope_details(
        policy_provenance=policy_provenance,
        ruv_readiness_enabled=ruv_readiness_enabled,
        ruv_readiness_ready=ruv_readiness_ready,
    )
    if parity_details is not None:
        caveats.append(
            ResultCaveat(
                code=DIFFERENTIAL_NARROW_PARITY_ENVELOPE_CAVEAT_CODE,
                severity="info",
                message=(
                    _narrow_parity_envelope_message(policy_provenance=policy_provenance)
                ),
                details=parity_details,
            )
        )

    return tuple(caveats)


def finalize_differential_result_caveats(
    *,
    caveats: tuple[ResultCaveat, ...],
    imputation_policy_inputs: DifferentialImputationPolicyInputs | None,
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs | None,
) -> tuple[ResultCaveat, ...]:
    """Refresh row-count caveats after final model-fit eligibility is known."""

    finalized: list[ResultCaveat] = []
    saw_imputation_caveat = False
    saw_withheld_caveat = False
    withheld_feature_details = _withheld_feature_details(feature_eligibility_inputs)
    for caveat in caveats:
        if caveat.code == DIFFERENTIAL_IMPUTATION_WITHHOLDING_POLICY_CAVEAT_CODE:
            saw_imputation_caveat = True
            if imputation_policy_inputs is not None:
                finalized.append(
                    _imputation_policy_caveat(
                        imputation_policy_inputs=imputation_policy_inputs,
                        feature_eligibility_inputs=feature_eligibility_inputs,
                    )
                )
            continue
        if caveat.code == DIFFERENTIAL_WITHHELD_FEATURES_CAVEAT_CODE:
            saw_withheld_caveat = True
            if withheld_feature_details is not None:
                finalized.append(
                    _withheld_feature_caveat(details=withheld_feature_details)
                )
            continue
        finalized.append(caveat)
    if imputation_policy_inputs is not None and not saw_imputation_caveat:
        finalized.append(
            _imputation_policy_caveat(
                imputation_policy_inputs=imputation_policy_inputs,
                feature_eligibility_inputs=feature_eligibility_inputs,
            )
        )
    if withheld_feature_details is not None and not saw_withheld_caveat:
        finalized.append(_withheld_feature_caveat(details=withheld_feature_details))
    return tuple(finalized)


def _exploratory_single_replicate_details(
    *,
    config: DifferentialAnalysisConfig,
    policy_provenance: DifferentialPolicyProvenance | None,
) -> dict[str, object] | None:
    if (
        config.reliability_profile
        != DIFFERENTIAL_RELIABILITY_PROFILE_EXPLORATORY_SINGLE_REPLICATE
    ):
        return None
    minimum_condition_replicates = 1
    condition_replicate_counts: list[dict[str, object]] = []
    contrasted_conditions_below_production: list[str] = []
    if policy_provenance is not None:
        minimum_condition_replicates = (
            policy_provenance.replicates.minimum_condition_replicates
        )
        condition_replicate_counts = [
            {"condition": condition, "biological_replicates": int(count)}
            for condition, count in policy_provenance.replicates.condition_replicate_counts
        ]
        contrasted_condition_names = {
            condition
            for contrast in policy_provenance.contrasts
            for condition in (
                contrast.numerator_condition,
                contrast.denominator_condition,
            )
        }
        contrasted_conditions_below_production = [
            condition
            for condition, count in policy_provenance.replicates.condition_replicate_counts
            if condition in contrasted_condition_names
            and int(count) < DIFFERENTIAL_PRODUCTION_MINIMUM_CONDITION_REPLICATES
        ]
    return {
        "reliability_profile": config.reliability_profile,
        "inferential_support": "exploratory_only",
        "computable_model_output": True,
        "production_supported_inference": False,
        "minimum_condition_replicates": int(minimum_condition_replicates),
        "production_minimum_condition_replicates": (
            DIFFERENTIAL_PRODUCTION_MINIMUM_CONDITION_REPLICATES
        ),
        "condition_replicate_counts": condition_replicate_counts,
        "contrasted_conditions_below_production_minimum": (
            contrasted_conditions_below_production
        ),
        "override_config": "reliability_profile",
    }


def _declared_scale_override_details(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    config: DifferentialAnalysisConfig,
) -> dict[str, object] | None:
    if not config.allow_suspicious_declared_input_scale:
        return None
    provenance = dataset.intensity_scale_state.establishment_provenance
    if provenance is None or _enum_value(provenance.mode) != "declared":
        return None
    warnings = tuple(str(value) for value in provenance.diagnostic_warnings)
    if not warnings:
        return None
    details: dict[str, object] = {
        "scale": str(provenance.scale),
        "establishment_mode": _enum_value(provenance.mode),
        "establishment_source": _enum_value(provenance.source),
        "diagnostic_warning_count": len(warnings),
        "first_diagnostic_warning": warnings[0],
        "override_config": "allow_suspicious_declared_input_scale",
    }
    if provenance.input_declaration_source is not None:
        details["input_declaration_source"] = provenance.input_declaration_source
    return details


def _imputation_policy_details(
    inputs: DifferentialImputationPolicyInputs,
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs | None = None,
) -> dict[str, object]:
    summary = summarize_differential_imputation_inference(
        imputation_policy_inputs=inputs,
        feature_eligibility_inputs=feature_eligibility_inputs,
    )
    summary_payload = imputation_inference_summary_payload(summary)
    return {
        "policy": inputs.policy,
        "imputed_value_max_fraction": float(inputs.max_fraction),
        "testable_feature_count": int(summary.tested_feature_count),
        "result_status_column": DIFFERENTIAL_RESULT_STATUS_COLUMN,
        "result_status_reason_column": DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
        **summary_payload,
    }


def _imputation_policy_caveat(
    *,
    imputation_policy_inputs: DifferentialImputationPolicyInputs,
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs | None,
) -> ResultCaveat:
    summary = summarize_differential_imputation_inference(
        imputation_policy_inputs=imputation_policy_inputs,
        feature_eligibility_inputs=feature_eligibility_inputs,
    )
    details = _imputation_policy_details(
        imputation_policy_inputs,
        feature_eligibility_inputs,
    )
    return ResultCaveat(
        code=DIFFERENTIAL_IMPUTATION_WITHHOLDING_POLICY_CAVEAT_CODE,
        severity=(
            "warning"
            if summary.withheld_feature_count > 0
            or summary.tested_imputed_feature_count > 0
            else "info"
        ),
        message=_imputation_policy_message(summary),
        details=details,
    )


def _imputation_policy_message(
    summary: DifferentialImputationInferenceSummary,
) -> str:
    tested_imputed_feature_count = int(summary.tested_imputed_feature_count)
    withheld_feature_count = int(summary.withheld_feature_count)
    if tested_imputed_feature_count > 0:
        return (
            "Differential analysis retained tested rows that include imputed "
            "values under the imputation-aware withholding policy. The fit used "
            "analysis-ready values, not observed-only fitting, and residual "
            "degrees of freedom were not adjusted for imputation. Withheld "
            "features were excluded from model fitting and multiple-testing "
            "adjustment denominators."
        )
    if withheld_feature_count > 0:
        return (
            "Differential analysis used the imputation-aware withholding policy; "
            "no tested rows contained imputed values after withholding. Withheld "
            "features were excluded from model fitting and multiple-testing "
            "adjustment denominators."
        )
    return (
        "Differential analysis used the imputation-aware withholding policy; no "
        "analysed rows contained imputed values and no rows were withheld."
    )


def _withheld_feature_caveat(
    *,
    details: dict[str, object],
) -> ResultCaveat:
    return ResultCaveat(
        code=DIFFERENTIAL_WITHHELD_FEATURES_CAVEAT_CODE,
        severity="warning",
        message=(
            "One or more features were withheld from differential model "
            "fitting; result tables contain row-level status and reason "
            "columns."
        ),
        details=details,
    )


def _duplicate_correlation_consensus_caveat() -> ResultCaveat:
    return ResultCaveat(
        code=DIFFERENTIAL_DUPLICATE_CORRELATION_CONSENSUS_CAVEAT_CODE,
        severity="info",
        message=(
            "Differential analysis used one consensus compound-symmetry "
            "within-block correlation for GLS; it did not fit feature-specific "
            "random effects."
        ),
        details={
            "paired_design_policy": PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
            "block_treatment": "consensus_correlation",
            "covariance_structure": "compound_symmetry",
            "feature_specific_random_effects": False,
        },
    )


def _withheld_feature_details(
    inputs: DifferentialFeatureEligibilityInputs | None,
) -> dict[str, object] | None:
    if inputs is None:
        return None
    total_feature_count = int(inputs.result_status.size)
    testable_feature_count = len(inputs.testable_feature_ids)
    withheld_feature_count = total_feature_count - testable_feature_count
    if withheld_feature_count <= 0:
        return None
    return {
        "total_feature_count": total_feature_count,
        "testable_feature_count": testable_feature_count,
        "withheld_feature_count": withheld_feature_count,
        "status_counts": _status_counts(inputs.result_status),
        "result_status_column": DIFFERENTIAL_RESULT_STATUS_COLUMN,
        "result_status_reason_column": DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    }


def _narrow_parity_envelope_details(
    *,
    policy_provenance: DifferentialPolicyProvenance | None,
    ruv_readiness_enabled: bool,
    ruv_readiness_ready: bool,
) -> dict[str, object] | None:
    if policy_provenance is None:
        return None
    unsupported_count = len(
        policy_provenance.unsupported_design.intentionally_rejected_features
    )
    design_limitation_count = len(policy_provenance.design.limitations)
    if (
        unsupported_count == 0
        and design_limitation_count == 0
        and not ruv_readiness_enabled
        and not ruv_readiness_ready
    ):
        return None
    return {
        "scope": "tested_design_and_contrast_envelope",
        "model_type": (
            "moderated_gls_duplicate_correlation"
            if policy_provenance.design.paired_design_policy
            == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
            else "moderated_ols_fixed_effect"
        ),
        "full_limma_or_phosr_parity_claimed": False,
        "paired_design_policy": policy_provenance.design.paired_design_policy,
        "unsupported_design_feature_count": unsupported_count,
        "design_limitation_count": design_limitation_count,
        "fixed_block_terms_present": bool(policy_provenance.design.block_columns),
        "ruv_readiness_metadata_present": bool(
            ruv_readiness_enabled or ruv_readiness_ready
        ),
    }


def _narrow_parity_envelope_message(
    *,
    policy_provenance: DifferentialPolicyProvenance | None,
) -> str:
    if (
        policy_provenance is not None
        and policy_provenance.design.paired_design_policy
        == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
    ):
        return (
            "Differential analysis used the supported duplicate-correlation "
            "moderated GLS envelope; it is not full limma or PhosR parity, and "
            "unsupported mixed-effect, random-slope, and RUV/SPS models were not "
            "fit."
        )
    return (
        "Differential analysis used the supported scoped fixed-effect moderated "
        "OLS envelope; it is not full limma or PhosR parity, and unsupported "
        "repeated-measure, mixed-effect, and RUV/SPS models were not fit; "
        "duplicate_correlation was not selected."
    )


def _status_counts(result_status: pd.Series) -> dict[str, int]:
    counts = result_status.astype(str).value_counts(sort=False)
    return {str(status): int(count) for status, count in counts.items()}


def _enum_value(value: object) -> str:
    enum_value = getattr(value, "value", value)
    return str(enum_value)


__all__ = [
    "DIFFERENTIAL_DECLARED_SCALE_OVERRIDE_CAVEAT_CODE",
    "DIFFERENTIAL_DIRECT_TRUSTED_DATASET_CAVEAT_CODE",
    "DIFFERENTIAL_DUPLICATE_CORRELATION_CONSENSUS_CAVEAT_CODE",
    "DIFFERENTIAL_EXPLORATORY_SINGLE_REPLICATE_CAVEAT_CODE",
    "DIFFERENTIAL_IMPUTATION_WITHHOLDING_POLICY_CAVEAT_CODE",
    "DIFFERENTIAL_NARROW_PARITY_ENVELOPE_CAVEAT_CODE",
    "DIFFERENTIAL_WITHHELD_FEATURES_CAVEAT_CODE",
    "build_differential_result_caveats",
    "finalize_differential_result_caveats",
]
