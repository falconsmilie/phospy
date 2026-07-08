"""Structured result caveats for differential workflow execution."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from phospy.contracts.configs import DifferentialAnalysisConfig
from phospy.contracts.result_caveats import ResultCaveat
from phospy.provenance.models import RunProvenance
from phospy.science.datasets.direct_construction import (
    DIRECT_CONSTRUCTION_SOURCE,
    DIRECT_CONSTRUCTION_WORKFLOW_NAME,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DifferentialPolicyProvenance,
)
from phospy.workflows.differential.models import (
    DifferentialFeatureEligibilityInputs,
    DifferentialImputationPolicyInputs,
)
from phospy.workflows.intensity_scale_evidence import (
    build_declared_input_intensity_scale_caveat,
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

    direct_construction_details = _direct_construction_details(dataset.provenance)
    if direct_construction_details is not None:
        caveats.append(
            ResultCaveat(
                code=DIFFERENTIAL_DIRECT_TRUSTED_DATASET_CAVEAT_CODE,
                severity="info",
                message=(
                    "Input dataset was directly constructed as trusted "
                    "analysis-ready state; differential analysis did not rerun "
                    "dataset-building validation."
                ),
                details=direct_construction_details,
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
            ResultCaveat(
                code=DIFFERENTIAL_IMPUTATION_WITHHOLDING_POLICY_CAVEAT_CODE,
                severity="warning",
                message=(
                    "Differential analysis used the imputation-aware withholding "
                    "policy; withheld features were excluded from model fitting "
                    "and multiple-testing adjustment denominators."
                ),
                details=_imputation_policy_details(imputation_policy_inputs),
            )
        )

    withheld_feature_details = _withheld_feature_details(feature_eligibility_inputs)
    if withheld_feature_details is not None:
        caveats.append(
            ResultCaveat(
                code=DIFFERENTIAL_WITHHELD_FEATURES_CAVEAT_CODE,
                severity="warning",
                message=(
                    "One or more features were withheld from differential model "
                    "fitting; result tables contain row-level status and reason "
                    "columns."
                ),
                details=withheld_feature_details,
            )
        )

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
                    "Differential analysis used the supported fixed-effect "
                    "moderated OLS envelope; unsupported repeated-measure, "
                    "mixed-effect, RUV/SPS, and duplicateCorrelation-style "
                    "models were not fit."
                ),
                details=parity_details,
            )
        )

    return tuple(caveats)


def _direct_construction_details(
    provenance: RunProvenance | None,
) -> dict[str, object] | None:
    if provenance is None:
        return None
    construction = _construction_payload(provenance.workflow_parameters)
    source = construction.get("source")
    if (
        provenance.workflow_name != DIRECT_CONSTRUCTION_WORKFLOW_NAME
        and source != DIRECT_CONSTRUCTION_SOURCE
    ):
        return None

    details: dict[str, object] = {
        "dataset_workflow_name": provenance.workflow_name,
        "construction_source": "" if source is None else str(source),
    }
    builder_used = construction.get("builder_used")
    if isinstance(builder_used, bool):
        details["builder_used"] = builder_used
    warning = construction.get("warning")
    if isinstance(warning, str) and warning.strip():
        details["construction_warning"] = warning.strip()
    return details


def _construction_payload(parameters: Mapping[str, object]) -> Mapping[str, object]:
    construction = parameters.get("construction")
    if isinstance(construction, Mapping):
        return construction
    return {}


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
) -> dict[str, object]:
    total_feature_count = int(inputs.result_status.size)
    testable_feature_count = len(inputs.testable_feature_ids)
    return {
        "policy": inputs.policy,
        "imputed_value_max_fraction": float(inputs.max_fraction),
        "total_feature_count": total_feature_count,
        "testable_feature_count": testable_feature_count,
        "withheld_feature_count": total_feature_count - testable_feature_count,
        "status_counts": _status_counts(inputs.result_status),
        "result_status_column": DIFFERENTIAL_RESULT_STATUS_COLUMN,
    }


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
        "model_type": "moderated_ols_fixed_effect",
        "paired_design_policy": policy_provenance.design.paired_design_policy,
        "unsupported_design_feature_count": unsupported_count,
        "design_limitation_count": design_limitation_count,
        "fixed_block_terms_present": bool(policy_provenance.design.block_columns),
        "ruv_readiness_metadata_present": bool(
            ruv_readiness_enabled or ruv_readiness_ready
        ),
    }


def _status_counts(result_status: pd.Series) -> dict[str, int]:
    counts = result_status.astype(str).value_counts(sort=False)
    return {str(status): int(count) for status, count in counts.items()}


def _enum_value(value: object) -> str:
    enum_value = getattr(value, "value", value)
    return str(enum_value)


__all__ = [
    "DIFFERENTIAL_DECLARED_SCALE_OVERRIDE_CAVEAT_CODE",
    "DIFFERENTIAL_DIRECT_TRUSTED_DATASET_CAVEAT_CODE",
    "DIFFERENTIAL_IMPUTATION_WITHHOLDING_POLICY_CAVEAT_CODE",
    "DIFFERENTIAL_NARROW_PARITY_ENVELOPE_CAVEAT_CODE",
    "DIFFERENTIAL_WITHHELD_FEATURES_CAVEAT_CODE",
    "build_differential_result_caveats",
]
