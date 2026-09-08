from __future__ import annotations

import json
from dataclasses import FrozenInstanceError, replace
from typing import Any, cast

import numpy as np
import pandas as pd
import pytest

from phospy.advanced import (
    DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1,
    PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
    DifferentialAnalysisConfig,
    DifferentialProteinAwareModelConfig,
)
from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.errors.input import PhosPyInputError
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.provenance.hashing import (
    fingerprint_optional_table_strict,
    fingerprint_table_strict,
)
from phospy.provenance.models import (
    EnvironmentProvenance,
    ReferenceProvenance,
    RunProvenance,
)
from phospy.provenance.reference_context import ReferenceContext
from phospy.provenance.serialization.tables import table_fingerprint_to_payload
from phospy.science.datasets.models import DatasetPreprocessingReport
from phospy.science.datasets.preprocessing.protein_aware_alignment import (
    ProteinAwarePreparationEligibility,
    ProteinAwareSampleAlignmentDiagnostics,
    ProteinAwareTransformationStateDiagnostics,
)
from phospy.science.datasets.preprocessing.protein_aware_preparation import (
    ProteinAwarePreparationReport,
    ProteinAwarePreparationResult,
    ProteinAwareSiteEligibility,
)
from phospy.science.datasets.preprocessing.protein_mapping import ProteinMappingStatus
from phospy.science.differential.models import (
    ContrastMatrix,
    DesignMatrix,
    DifferentialAnalysisResult,
    EmpiricalBayesPriorDiagnostics,
    ProteinAwareDifferentialDiagnostics,
)
from phospy.science.differential.models.protein_aware import (
    ProteinAwareDifferentialComputationRequest,
    ProteinAwareDifferentialComputationResult,
)
from phospy.science.differential.models.provenance import (
    DIFFERENTIAL_PROTEIN_AWARE_CLAIM_STATUS_EXPERIMENTAL,
    DIFFERENTIAL_PROTEIN_AWARE_CONDITION_NUMBER_SUMMARY_SCOPE,
    DIFFERENTIAL_PROTEIN_AWARE_CONDITION_NUMBER_SUMMARY_STATISTIC,
    DIFFERENTIAL_PROTEIN_AWARE_CONTRAST_MATRIX_FINGERPRINT_NAME,
    DIFFERENTIAL_PROTEIN_AWARE_COVARIATE_MATRIX_FINGERPRINT_NAME,
    DIFFERENTIAL_PROTEIN_AWARE_DESIGN_MATRIX_FINGERPRINT_NAME,
    DIFFERENTIAL_PROTEIN_AWARE_MATCHED_PAIRS_FINGERPRINT_NAME,
    DIFFERENTIAL_PROTEIN_AWARE_MODEL_FORMULA,
    DIFFERENTIAL_PROTEIN_AWARE_NUISANCE_COEFFICIENT_NAME,
    DIFFERENTIAL_PROTEIN_AWARE_PHOSPHO_MATRIX_FINGERPRINT_NAME,
    DIFFERENTIAL_PROTEIN_AWARE_POLICY_METHOD_VERSION,
    DIFFERENTIAL_PROTEIN_AWARE_SITE_ELIGIBILITY_FINGERPRINT_NAME,
    DifferentialPolicyProvenance,
    DifferentialProteinAwarePolicyProvenance,
)
from phospy.science.differential.models.tables import (
    DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED,
    DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK,
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
)
from phospy.workflows.differential.caveats import (
    DIFFERENTIAL_PROTEIN_AWARE_CONDITIONAL_EFFECT_CAVEAT_CODE,
    DIFFERENTIAL_PROTEIN_AWARE_EXPERIMENTAL_CAVEAT_CODE,
    DIFFERENTIAL_PROTEIN_AWARE_NO_FALLBACK_CAVEAT_CODE,
    DIFFERENTIAL_PROTEIN_AWARE_NO_PREPROCESSING_CAVEAT_CODE,
    DIFFERENTIAL_PROTEIN_AWARE_UNSUPPORTED_SCOPE_CAVEAT_CODE,
)
from phospy.workflows.differential.interpreter import DifferentialAnalysisInterpreter
from phospy.workflows.differential.models import (
    DifferentialFeatureEligibilityInputs,
    InterpretedDifferentialAnalysisRequest,
    ProteinAwareDifferentialResolvedInputs,
)
from phospy.workflows.differential.protein_aware_inputs import (
    ProteinAwareDifferentialInputResolver,
)
from phospy.workflows.differential.provenance import (
    build_protein_aware_policy_provenance,
)
from phospy.workflows.differential.result_assembly import DifferentialResultAssembler
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

_SAMPLE_IDS = ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")
_PROTEINS = ("MAPK14", "GSK3B", "AKT1", "PRKACA", "RPS6", "EIF4EBP1")
_SITES = ("Y182", "S9", "T308", "S339", "S235", "T37")
_TOTAL_ROW_KEYS = (
    "protein:MAPK14",
    None,
    None,
    "protein:AKT1",
    "protein:PRKACA",
    "protein:RPS6",
)


def test_protein_aware_assembly_expands_to_full_index_with_withholding_statuses() -> (
    None
):
    result = _protein_aware_result()
    table = result.table_for("B_vs_A")

    assert table.index.tolist() == _site_index().tolist()
    assert table["site_key"].tolist() == _site_index().tolist()
    assert table[DIFFERENTIAL_RESULT_STATUS_COLUMN].tolist() == [
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE,
    ]
    assert table[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN].tolist() == [
        "",
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE,
    ]
    assert np.isfinite(
        table.loc[[table.index[0]], ["logFC", "t", "P.Value", "adj.P.Val"]].to_numpy(
            dtype=float
        )
    ).all()
    assert table.iloc[1:][["logFC", "t", "P.Value", "adj.P.Val"]].isna().all().all()


def test_protein_aware_diagnostics_are_defensive_and_site_scoped() -> None:
    result = _protein_aware_result()
    diagnostics = result.protein_aware_diagnostics
    assert diagnostics is not None

    assert diagnostics.method_id == (
        DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
    )
    assert diagnostics.claim_status == "experimental"
    assert diagnostics.preparation_policy == "prepare_model_inputs"
    assert diagnostics.protein_mapping_policy == "require_unambiguous"
    assert diagnostics.protein_covariate_centered is True
    assert diagnostics.protein_covariate_standardized is False
    assert diagnostics.minimum_condition_replicates == 2
    assert diagnostics.total_site_count == 6
    assert diagnostics.ordinary_eligible_site_count == 6
    assert diagnostics.protein_preparation_eligible_site_count == 4
    assert diagnostics.tested_site_count == 1
    assert diagnostics.withheld_site_count == 5
    assert diagnostics.execution_sample_order == _SAMPLE_IDS
    assert diagnostics.base_design_rank == 2
    assert diagnostics.expected_augmented_rank == 3
    assert diagnostics.common_augmented_rank == 3
    assert diagnostics.common_augmented_residual_degrees_of_freedom == pytest.approx(
        3.0
    )
    assert diagnostics.distinct_matched_protein_row_count == 4
    assert diagnostics.fitted_protein_row_count == 1
    assert diagnostics.median_augmented_condition_number == pytest.approx(11.0)
    assert result.diagnostics.singular_values == ()
    assert result.diagnostics.condition_number == pytest.approx(11.0)
    assert result.diagnostics.model_type == "protein_covariate_adjusted_moderated_ols"

    per_site = diagnostics.per_site_diagnostics_dataframe()
    assert per_site.index.tolist() == _site_index().tolist()
    assert per_site.loc[
        per_site.index[0], "protein_covariate_coefficient"
    ] == pytest.approx(0.25)
    assert per_site.iloc[1:]["protein_covariate_coefficient"].isna().all()
    assert (
        per_site.loc[
            per_site.index[-1],
            "protein_contrast_estimability_status",
        ]
        == "non_estimable"
    )

    per_site.loc[per_site.index[0], "protein_covariate_coefficient"] = -99.0
    exported_again = diagnostics.per_site_diagnostics_dataframe()
    assert exported_again.loc[
        exported_again.index[0],
        "protein_covariate_coefficient",
    ] == pytest.approx(0.25)

    via_result = result.protein_aware_site_diagnostics_dataframe()
    assert via_result is not None
    via_result.loc[via_result.index[0], "protein_covariate_coefficient"] = -99.0
    exported_third = result.protein_aware_site_diagnostics_dataframe()
    assert exported_third is not None
    assert exported_third.loc[
        exported_third.index[0],
        "protein_covariate_coefficient",
    ] == pytest.approx(0.25)


def test_protein_aware_scientific_equality_uses_diagnostics_content() -> None:
    left = _protein_aware_result()
    same = _protein_aware_result()
    different = _protein_aware_result(protein_coefficient=0.40)

    assert left.scientifically_equals(same)
    assert not left.scientifically_equals(different)


def test_protein_aware_payload_is_json_compatible_and_deterministic() -> None:
    result = _protein_aware_result()

    first_payload = result.to_payload()
    second_payload = result.to_payload()
    repeated_payload = _protein_aware_result().to_payload()

    assert first_payload == second_payload
    assert (
        cast(dict[str, object], first_payload["policy_provenance"])["protein_aware"]
        == cast(dict[str, object], repeated_payload["policy_provenance"])[
            "protein_aware"
        ]
    )
    assert "protein_aware_diagnostics" in first_payload
    protein_aware_diagnostics_payload = cast(
        dict[str, object],
        first_payload["protein_aware_diagnostics"],
    )
    assert protein_aware_diagnostics_payload["minimum_condition_replicates"] == 2
    json.dumps(first_payload, sort_keys=True)
    serialized_keys = json.dumps(first_payload, sort_keys=True)
    assert "protein_covariate_p_value" not in serialized_keys
    assert "protein_covariate_q_value" not in serialized_keys


def test_protein_aware_policy_provenance_records_pa007_contract() -> None:
    result = _protein_aware_result()
    assert result.policy_provenance is not None
    protein_aware = result.policy_provenance.protein_aware
    assert protein_aware is not None

    assert protein_aware.method_id == (
        DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
    )
    assert protein_aware.method_version == (
        DIFFERENTIAL_PROTEIN_AWARE_POLICY_METHOD_VERSION
    )
    assert protein_aware.claim_status == (
        DIFFERENTIAL_PROTEIN_AWARE_CLAIM_STATUS_EXPERIMENTAL
    )
    assert protein_aware.model_formula == DIFFERENTIAL_PROTEIN_AWARE_MODEL_FORMULA
    assert protein_aware.nuisance_coefficient_name == (
        DIFFERENTIAL_PROTEIN_AWARE_NUISANCE_COEFFICIENT_NAME
    )
    assert protein_aware.preparation_schema_version == 1
    assert protein_aware.preparation_policy == "prepare_model_inputs"
    assert protein_aware.protein_mapping_policy == "require_unambiguous"
    assert protein_aware.protein_reference_context["availability"] == "unavailable"
    assert protein_aware.phosphosite_transformation_state["availability"] == (
        "unavailable"
    )
    assert protein_aware.total_protein_transformation_state["availability"] == (
        "unavailable"
    )
    assert protein_aware.prior_total_protein_correction_state["applied"] is False
    assert protein_aware.phosphosite_normalisation_state == "none"
    assert protein_aware.protein_covariate_centered is True
    assert protein_aware.protein_covariate_standardized is False
    assert protein_aware.automatic_protein_normalization is False
    assert protein_aware.protein_imputation is False
    assert protein_aware.phosphosite_only_fallback is False
    assert protein_aware.execution_sample_order == _SAMPLE_IDS
    assert protein_aware.total_site_count == 6
    assert protein_aware.ordinary_eligible_site_count == 6
    assert protein_aware.protein_preparation_eligible_site_count == 4
    assert protein_aware.distinct_matched_protein_row_count == 4
    assert protein_aware.fitted_protein_row_count == 1
    assert protein_aware.tested_site_count == 1
    assert protein_aware.withheld_site_count == 5
    assert protein_aware.common_augmented_residual_degrees_of_freedom == (
        pytest.approx(3.0)
    )
    assert protein_aware.condition_number_summary_scope == (
        DIFFERENTIAL_PROTEIN_AWARE_CONDITION_NUMBER_SUMMARY_SCOPE
    )
    assert protein_aware.condition_number_summary_statistic == (
        DIFFERENTIAL_PROTEIN_AWARE_CONDITION_NUMBER_SUMMARY_STATISTIC
    )
    assert protein_aware.median_augmented_condition_number == pytest.approx(11.0)
    assert "MSstatsPTM parity" in protein_aware.unsupported_claims

    fingerprints = protein_aware.input_fingerprints
    assert fingerprints.phospho_matrix.name == (
        DIFFERENTIAL_PROTEIN_AWARE_PHOSPHO_MATRIX_FINGERPRINT_NAME
    )
    assert fingerprints.protein_matched_pairs.name == (
        DIFFERENTIAL_PROTEIN_AWARE_MATCHED_PAIRS_FINGERPRINT_NAME
    )
    assert fingerprints.protein_covariate_matrix.name == (
        DIFFERENTIAL_PROTEIN_AWARE_COVARIATE_MATRIX_FINGERPRINT_NAME
    )
    assert fingerprints.protein_site_eligibility.name == (
        DIFFERENTIAL_PROTEIN_AWARE_SITE_ELIGIBILITY_FINGERPRINT_NAME
    )
    assert fingerprints.design_matrix.name == (
        DIFFERENTIAL_PROTEIN_AWARE_DESIGN_MATRIX_FINGERPRINT_NAME
    )
    assert fingerprints.contrast_matrix.name == (
        DIFFERENTIAL_PROTEIN_AWARE_CONTRAST_MATRIX_FINGERPRINT_NAME
    )
    assert fingerprints.phospho_matrix.rows == 1
    assert fingerprints.protein_covariate_matrix.columns == len(_SAMPLE_IDS)

    payload = result.to_payload()
    policy_payload = cast(dict[str, object], payload["policy_provenance"])
    protein_aware_payload = cast(dict[str, object], policy_payload["protein_aware"])
    assert protein_aware_payload["model_formula"] == (
        DIFFERENTIAL_PROTEIN_AWARE_MODEL_FORMULA
    )
    json.dumps(payload, sort_keys=True)


def test_protein_aware_policy_provenance_records_supplied_preparation_context() -> None:
    dataset = _dataset_with_trusted_protein_aware_preparation_provenance(
        include_reference_context=True
    )
    request = _request(dataset=dataset, config=_protein_aware_config())
    validated = DifferentialAnalysisValidator().run(request)
    interpreted = DifferentialAnalysisInterpreter().run(
        DifferentialAnalysisValidator().run(_request(dataset=dataset))
    )
    resolved_inputs = ProteinAwareDifferentialInputResolver().run(
        validated,
        minimum_condition_replicates=2,
    )
    result = DifferentialResultAssembler().run_protein_aware(
        request=interpreted,
        resolved_inputs=resolved_inputs,
        computation_result=_computation_result(
            interpreted,
            resolved_inputs,
            protein_coefficient=0.25,
        ),
        workflow_provenance={"stage": "pa007-trusted-context-test"},
    )
    protein_aware = _required_protein_aware_policy(result)

    reference_context = protein_aware.protein_reference_context
    assert reference_context["availability"] == "available"
    assert reference_context["source"] == (
        "protein_aware_preparation.provenance.reference_context"
    )
    reference_context_payload = cast(
        dict[str, object],
        reference_context["reference_context"],
    )
    assert reference_context_payload["organism"] == "rat"
    assert reference_context_payload["protein_namespace"] == "test-protein-id"
    assert reference_context_payload["source_name"] == "unit-test-reference"
    assert reference_context_payload["source_version"] == "2026.1"

    assert (
        protein_aware.protein_mapping_policy_parameters["allow_reordered_samples"]
        is False
    )
    assert "dataset_binding_table_fingerprints" in (
        protein_aware.protein_mapping_policy_parameters
    )
    assert protein_aware.phosphosite_transformation_state["availability"] == (
        "available"
    )
    assert protein_aware.phosphosite_transformation_state["established_by"] == (
        "test.phospho"
    )
    assert protein_aware.total_protein_transformation_state["availability"] == (
        "available"
    )
    assert protein_aware.total_protein_transformation_state["established_by"] == (
        "test.total"
    )
    assert protein_aware.prior_total_protein_correction_state["availability"] == (
        "available"
    )
    assert protein_aware.prior_total_protein_correction_state["policy"] == "none"
    assert protein_aware.prior_total_protein_correction_state["applied"] is False


def test_protein_aware_policy_provenance_preserves_reference_payload_when_supplied() -> (
    None
):
    dataset = _dataset_with_trusted_protein_aware_preparation_provenance(
        include_reference_context=False
    )
    request = _request(dataset=dataset, config=_protein_aware_config())
    validated = DifferentialAnalysisValidator().run(request)
    interpreted = DifferentialAnalysisInterpreter().run(
        DifferentialAnalysisValidator().run(_request(dataset=dataset))
    )
    resolved_inputs = ProteinAwareDifferentialInputResolver().run(
        validated,
        minimum_condition_replicates=2,
    )
    result = DifferentialResultAssembler().run_protein_aware(
        request=interpreted,
        resolved_inputs=resolved_inputs,
        computation_result=_computation_result(
            interpreted,
            resolved_inputs,
            protein_coefficient=0.25,
        ),
        workflow_provenance={"stage": "pa007-reference-payload-test"},
    )
    protein_aware = _required_protein_aware_policy(result)

    reference_context = protein_aware.protein_reference_context
    assert reference_context["availability"] == "available"
    assert (
        reference_context["source"] == "protein_aware_preparation.provenance.reference"
    )
    reference_payload = cast(dict[str, object], reference_context["reference"])
    assert reference_payload["source_type"] == "local_fixture"
    assert reference_payload["source_name"] == "unit-test-reference"
    assert reference_payload["source_version"] == "2026.1"
    assert reference_payload["identifier_namespace"] == "test-protein-id"
    assert reference_payload["sequence_window"] == {"center_residue": "S/T/Y"}
    assert reference_payload["manifest"] == {
        "organism": "rat",
        "source_name": "unit-test-reference",
        "source_version": "2026.1",
    }
    assert reference_payload["reference_context"] is not None
    assert reference_payload["table_fingerprints"]


def test_protein_aware_policy_provenance_is_immutable_and_comparable() -> None:
    left = _protein_aware_result()
    right = _protein_aware_result()
    assert left.policy_provenance is not None
    assert right.policy_provenance is not None
    assert left.policy_provenance.protein_aware is not None
    assert right.policy_provenance.protein_aware is not None

    assert left.policy_provenance.protein_aware == (
        right.policy_provenance.protein_aware
    )
    with pytest.raises(FrozenInstanceError):
        left.policy_provenance.protein_aware.method_id = "changed"  # type: ignore[misc]
    with pytest.raises(TypeError):
        mapping = cast(
            Any,
            left.policy_provenance.protein_aware.protein_mapping_policy_parameters,
        )
        mapping["changed"] = True


@pytest.mark.parametrize(
    ("field_name", "value", "message"),
    (
        ("claim_status", "validated", "claim_status"),
        ("model_formula", "y = X beta + error", "model_formula"),
        ("nuisance_coefficient_name", "total_protein", "nuisance_coefficient_name"),
        ("protein_covariate_centered", False, "protein_covariate_centered"),
        ("protein_covariate_standardized", True, "protein_covariate_standardized"),
        (
            "automatic_protein_normalization",
            True,
            "automatic_protein_normalization",
        ),
        ("protein_imputation", True, "protein_imputation"),
        ("phosphosite_only_fallback", True, "phosphosite_only_fallback"),
    ),
)
def test_protein_aware_policy_provenance_rejects_adr_inconsistent_values(
    field_name: str,
    value: object,
    message: str,
) -> None:
    protein_aware = _required_protein_aware_policy()

    with pytest.raises(PhosPyInputError, match=message):
        _protein_aware_policy_with_override(
            protein_aware,
            field_name=field_name,
            value=value,
        )


def test_protein_aware_policy_provenance_rejects_non_json_mapping() -> None:
    protein_aware = _required_protein_aware_policy()

    with pytest.raises(PhosPyInputError, match="protein_mapping_policy_parameters"):
        replace(
            protein_aware,
            protein_mapping_policy_parameters={"not_json": object()},
        )


def test_protein_aware_policy_provenance_rejects_inconsistent_counts() -> None:
    protein_aware = _required_protein_aware_policy()

    with pytest.raises(PhosPyInputError, match="status_counts"):
        replace(
            protein_aware,
            status_counts=(
                (DIFFERENTIAL_RESULT_STATUS_TESTED, 1),
                (
                    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
                    1,
                ),
            ),
        )


def test_protein_aware_input_fingerprints_reject_unstable_names() -> None:
    protein_aware = _required_protein_aware_policy()
    fingerprints = protein_aware.input_fingerprints

    with pytest.raises(PhosPyInputError, match="phospho_matrix"):
        replace(
            fingerprints,
            phospho_matrix=replace(
                fingerprints.phospho_matrix,
                name="differential.input.not_the_phospho_matrix",
            ),
        )


def test_protein_aware_policy_rejects_duplicate_correlation_combination() -> None:
    protein_aware = _required_protein_aware_policy()
    duplicate_result = DifferentialAnalysisWorkflow().run(
        _request(
            config=DifferentialAnalysisConfig(
                paired_design_policy=PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
            ),
            include_blocks=True,
        )
    )
    duplicate_policy = duplicate_result.policy_provenance
    assert duplicate_policy is not None

    with pytest.raises(PhosPyInputError, match="duplicate_correlation"):
        replace(
            cast(DifferentialPolicyProvenance, duplicate_policy),
            protein_aware=protein_aware,
        )


def test_protein_aware_policy_fingerprints_change_for_material_input_mutations() -> (
    None
):
    interpreted = _interpreted_request()
    resolved_inputs = _resolved_inputs(interpreted, protein_coefficient=0.25)
    base_hashes = _protein_aware_policy_fingerprint_hashes(resolved_inputs)

    mutation_cases = {
        "phospho_matrix": "phospho_matrix",
        "protein_matched_pairs": "protein_matched_pairs",
        "protein_covariate_matrix": "protein_covariate_matrix",
        "protein_site_eligibility": "protein_site_eligibility",
        "design_matrix": "design_matrix",
        "contrast_matrix": "contrast_matrix",
    }
    for mutation, fingerprint_name in mutation_cases.items():
        mutated_inputs = _resolved_inputs_with_fingerprint_mutation(
            resolved_inputs,
            mutation=mutation,
        )
        mutated_hashes = _protein_aware_policy_fingerprint_hashes(mutated_inputs)
        assert mutated_hashes[fingerprint_name] != base_hashes[fingerprint_name]


def test_protein_aware_policy_fingerprints_respect_execution_sample_order() -> None:
    interpreted = _interpreted_request()
    resolved_inputs = _resolved_inputs(interpreted, protein_coefficient=0.25)
    base_hashes = _protein_aware_policy_fingerprint_hashes(resolved_inputs)
    reordered_inputs = _resolved_inputs_with_fingerprint_mutation(
        resolved_inputs,
        mutation="sample_order",
    )
    reordered_hashes = _protein_aware_policy_fingerprint_hashes(reordered_inputs)

    assert reordered_inputs.sample_order == tuple(reversed(_SAMPLE_IDS))
    assert reordered_hashes["phospho_matrix"] != base_hashes["phospho_matrix"]
    assert (
        reordered_hashes["protein_covariate_matrix"]
        != base_hashes["protein_covariate_matrix"]
    )
    assert reordered_hashes["design_matrix"] != base_hashes["design_matrix"]


def test_protein_aware_caveats_are_present_only_for_selected_lane() -> None:
    protein_aware_result = _protein_aware_result()
    ordinary_result = DifferentialAnalysisWorkflow().run(_request())
    duplicate_result = DifferentialAnalysisWorkflow().run(
        _request(
            config=DifferentialAnalysisConfig(
                paired_design_policy=PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
            ),
            include_blocks=True,
        )
    )
    protein_aware_codes = {caveat.code for caveat in protein_aware_result.caveats}
    ordinary_codes = {caveat.code for caveat in ordinary_result.caveats}
    duplicate_codes = {caveat.code for caveat in duplicate_result.caveats}
    expected_codes = {
        DIFFERENTIAL_PROTEIN_AWARE_EXPERIMENTAL_CAVEAT_CODE,
        DIFFERENTIAL_PROTEIN_AWARE_CONDITIONAL_EFFECT_CAVEAT_CODE,
        DIFFERENTIAL_PROTEIN_AWARE_NO_PREPROCESSING_CAVEAT_CODE,
        DIFFERENTIAL_PROTEIN_AWARE_NO_FALLBACK_CAVEAT_CODE,
        DIFFERENTIAL_PROTEIN_AWARE_UNSUPPORTED_SCOPE_CAVEAT_CODE,
    }

    assert expected_codes <= protein_aware_codes
    assert ordinary_codes.isdisjoint(expected_codes)
    assert duplicate_codes.isdisjoint(expected_codes)
    conditional_caveat = next(
        caveat
        for caveat in protein_aware_result.caveats
        if caveat.code == DIFFERENTIAL_PROTEIN_AWARE_CONDITIONAL_EFFECT_CAVEAT_CODE
    )
    assert "not stoichiometry or occupancy estimates" in conditional_caveat.message


def test_ordinary_payload_structure_has_no_protein_aware_diagnostics() -> None:
    result = DifferentialAnalysisWorkflow().run(_request())
    payload = result.to_payload()

    assert result.protein_aware_diagnostics is None
    assert result.protein_aware_site_diagnostics_dataframe() is None
    assert set(payload) == {
        "caveats",
        "diagnostics",
        "workflow_provenance",
        "policy_provenance",
        "empirical_bayes",
        "contrast_tables",
        "feature_eligibility",
    }
    assert "protein_aware_diagnostics" not in payload
    policy_payload = cast(dict[str, object], payload["policy_provenance"])
    assert "protein_aware" not in policy_payload


def test_duplicate_correlation_result_assembly_has_no_protein_aware_diagnostics() -> (
    None
):
    result = DifferentialAnalysisWorkflow().run(
        _request(
            config=DifferentialAnalysisConfig(
                paired_design_policy=PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
            ),
            include_blocks=True,
        )
    )

    assert result.protein_aware_diagnostics is None
    assert result.protein_aware_site_diagnostics_dataframe() is None
    assert "protein_aware_diagnostics" not in result.to_payload()
    policy_payload = cast(
        dict[str, object],
        result.to_payload()["policy_provenance"],
    )
    assert "protein_aware" not in policy_payload
    assert result.diagnostics.model_type == "moderated_gls_duplicate_correlation"


def test_protein_aware_result_rejects_misaligned_diagnostics_index() -> None:
    result = _protein_aware_result()
    diagnostics = _required_protein_aware_diagnostics(result)
    per_site = diagnostics.per_site_diagnostics_dataframe()
    original_first_site = str(per_site.index[0])
    mismatched_first_site = f"{original_first_site}|different"
    per_site = per_site.rename(index={original_first_site: mismatched_first_site})
    per_site.loc[mismatched_first_site, "site_key"] = mismatched_first_site
    misaligned_diagnostics = _protein_aware_diagnostics_copy(
        diagnostics,
        per_site_diagnostics=per_site,
    )

    with pytest.raises(PhosPyInputError, match="per-site diagnostics index"):
        _result_with_protein_aware_diagnostics(
            result,
            protein_aware_diagnostics=misaligned_diagnostics,
        )


@pytest.mark.parametrize(
    ("field_name", "value", "message"),
    (
        ("claim_status", "validated", "claim_status"),
        ("method_id", "different_method", "method_id"),
        ("protein_covariate_centered", False, "protein_covariate_centered"),
        ("protein_covariate_standardized", True, "protein_covariate_standardized"),
        (
            "protein_covariate_centering_policy",
            "z_scored",
            "protein_covariate_centering_policy",
        ),
        ("protein_covariate_imputation_policy", "mean", "imputation_policy"),
        ("fallback_policy", "fallback_to_ordinary", "fallback_policy"),
        ("minimum_condition_replicates", 0, "minimum_condition_replicates"),
    ),
)
def test_protein_aware_diagnostics_rejects_adr_inconsistent_contract_values(
    field_name: str,
    value: object,
    message: str,
) -> None:
    diagnostics = _required_protein_aware_diagnostics(_protein_aware_result())

    with pytest.raises(PhosPyInputError, match=message):
        _protein_aware_diagnostics_copy_with_override(
            diagnostics,
            field_name=field_name,
            value=value,
        )


@pytest.mark.parametrize(
    ("mismatch", "seam"),
    (
        ("sample_order", "request_sample_order"),
        ("base_design", "request_base_design"),
        ("base_contrasts", "request_base_contrasts"),
    ),
)
def test_protein_aware_assembly_rejects_resolved_inputs_from_different_request(
    mismatch: str,
    seam: str,
) -> None:
    interpreted = _interpreted_request()
    resolved_inputs = _resolved_inputs(interpreted, protein_coefficient=0.25)
    mismatched_resolved_inputs = _resolved_inputs_with_request_mismatch(
        resolved_inputs,
        mismatch=mismatch,
    )
    computation_result = _computation_result(
        interpreted,
        mismatched_resolved_inputs,
        protein_coefficient=0.25,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        DifferentialResultAssembler().run_protein_aware(
            request=interpreted,
            resolved_inputs=mismatched_resolved_inputs,
            computation_result=computation_result,
            workflow_provenance={"stage": "pa006-test"},
        )

    assert exc_info.value.seam is not None
    assert exc_info.value.seam.endswith(seam)


def test_protein_aware_assembly_carries_computation_time_failure_diagnostics() -> None:
    interpreted = _interpreted_request()
    resolved_inputs = _resolved_inputs(
        interpreted,
        protein_coefficient=0.25,
        tested_positions=(0, 5),
    )
    failed_site_id = resolved_inputs.tested_site_ids[1]
    computation_result = _computation_result(
        interpreted,
        resolved_inputs,
        protein_coefficient=0.25,
        failure_site_id=failed_site_id,
        failure_total_protein_row_key="protein:RPS6",
    )

    result = DifferentialResultAssembler().run_protein_aware(
        request=interpreted,
        resolved_inputs=resolved_inputs,
        computation_result=computation_result,
        workflow_provenance={"stage": "pa006-test"},
    )

    table = result.table_for("B_vs_A")
    assert table.loc[failed_site_id, DIFFERENTIAL_RESULT_STATUS_COLUMN] == (
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID
    )
    assert table.loc[failed_site_id, DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN] == (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT
    )
    assert (
        table.loc[[failed_site_id], ["logFC", "t", "P.Value", "adj.P.Val"]]
        .isna()
        .all()
        .all()
    )

    per_site = result.protein_aware_site_diagnostics_dataframe()
    assert per_site is not None
    assert per_site.loc[
        failed_site_id, "protein_augmented_design_rank"
    ] == pytest.approx(2.0)
    assert per_site.loc[
        failed_site_id,
        "protein_augmented_design_condition_number",
    ] == pytest.approx(1.0e11)
    assert pd.isna(per_site.loc[failed_site_id, "protein_covariate_coefficient"])
    assert (
        per_site.loc[failed_site_id, "protein_aware_failure_message"]
        == "protein augmented design is rank deficient"
    )


def _protein_aware_result(
    *,
    protein_coefficient: float = 0.25,
) -> DifferentialAnalysisResult:
    interpreted = _interpreted_request()
    resolved_inputs = _resolved_inputs(
        interpreted,
        protein_coefficient=protein_coefficient,
    )
    computation_result = _computation_result(
        interpreted,
        resolved_inputs,
        protein_coefficient=protein_coefficient,
    )
    return DifferentialResultAssembler().run_protein_aware(
        request=interpreted,
        resolved_inputs=resolved_inputs,
        computation_result=computation_result,
        workflow_provenance={"stage": "pa006-test"},
    )


def _resolved_inputs(
    interpreted: InterpretedDifferentialAnalysisRequest,
    *,
    protein_coefficient: float,
    tested_positions: tuple[int, ...] = (0,),
) -> ProteinAwareDifferentialResolvedInputs:
    site_ids = tuple(str(value) for value in interpreted.result_identity_metadata.index)
    tested_site_ids = tuple(site_ids[position] for position in tested_positions)
    base_design = cast(DesignMatrix, interpreted.computation_request.design)
    base_contrasts = cast(ContrastMatrix, interpreted.computation_request.contrasts)
    matched_pair_rows = [
        {
            "site_key": site_ids[position],
            "protein_identifier": str(_PROTEINS[position]),
            "total_protein_row_key": str(_TOTAL_ROW_KEYS[position]),
        }
        for position in tested_positions
    ]
    matched_pairs = pd.DataFrame(matched_pair_rows)
    covariate_values_by_row_key = {
        "protein:MAPK14": [10.0, 10.2, 10.4, 11.2, 11.4, 11.6],
        "protein:RPS6": [12.0, 12.3, 12.6, 13.2, 13.5, 13.8],
    }
    tested_total_row_keys = tuple(
        dict.fromkeys(str(row["total_protein_row_key"]) for row in matched_pair_rows)
    )
    covariates = pd.DataFrame(
        [covariate_values_by_row_key[row_key] for row_key in tested_total_row_keys],
        index=pd.Index(tested_total_row_keys, name="total_protein_row_key"),
        columns=pd.Index(_SAMPLE_IDS, name="sample_id"),
    )
    computation_request = ProteinAwareDifferentialComputationRequest(
        phosphosite_matrix=interpreted.computation_request.matrix.loc[
            list(tested_site_ids),
            list(_SAMPLE_IDS),
        ],
        base_design=base_design,
        base_contrasts=base_contrasts,
        sample_order=_SAMPLE_IDS,
        matched_pairs=matched_pairs,
        resolved_protein_covariates=covariates,
        empirical_bayes=interpreted.computation_request.empirical_bayes,
        multiple_testing_method=interpreted.computation_request.multiple_testing_method,
        method_id=(
            DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
        ),
    )
    feature_metadata = _protein_aware_feature_metadata(
        site_ids,
        protein_coefficient=protein_coefficient,
    )
    feature_metadata.loc[
        list(tested_site_ids),
        DIFFERENTIAL_RESULT_STATUS_COLUMN,
    ] = DIFFERENTIAL_RESULT_STATUS_TESTED
    feature_metadata.loc[
        list(tested_site_ids),
        DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    ] = ""
    feature_metadata.loc[list(tested_site_ids), "protein_aware_tested"] = True
    feature_eligibility = DifferentialFeatureEligibilityInputs(
        feature_metadata=feature_metadata,
        result_status=feature_metadata[DIFFERENTIAL_RESULT_STATUS_COLUMN],
        testable_feature_ids=tested_site_ids,
        attach_to_result_tables=True,
    )
    return ProteinAwareDifferentialResolvedInputs(
        computation_request=computation_request,
        feature_eligibility_inputs=feature_eligibility,
        matched_pairs=matched_pairs,
        candidate_matched_pairs=pd.DataFrame(
            [
                {
                    "site_key": site_ids[0],
                    "protein_identifier": "MAPK14",
                    "total_protein_row_key": "protein:MAPK14",
                },
                {
                    "site_key": site_ids[3],
                    "protein_identifier": "AKT1",
                    "total_protein_row_key": "protein:AKT1",
                },
                {
                    "site_key": site_ids[4],
                    "protein_identifier": "PRKACA",
                    "total_protein_row_key": "protein:PRKACA",
                },
                {
                    "site_key": site_ids[5],
                    "protein_identifier": "RPS6",
                    "total_protein_row_key": "protein:RPS6",
                },
            ]
        ),
        resolved_protein_covariates=covariates,
        site_eligibility_metadata=feature_metadata,
        full_site_ids=site_ids,
        tested_site_ids=tested_site_ids,
        sample_order=_SAMPLE_IDS,
        base_design=base_design,
        base_contrasts=base_contrasts,
        method_id=(
            DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
        ),
        preparation_policy="prepare_model_inputs",
        protein_mapping_policy="require_unambiguous",
        eligibility_counts=(
            ("total_site_count", 6),
            ("ordinary_testable_site_count", 6),
            ("protein_preparation_candidate_site_count", 4),
            ("protein_aware_tested_site_count", len(tested_site_ids)),
        ),
        status_counts=_status_counts(feature_metadata),
        reason_counts=_reason_counts(feature_metadata),
    )


def _protein_aware_policy_fingerprint_hashes(
    resolved_inputs: ProteinAwareDifferentialResolvedInputs,
) -> dict[str, str]:
    interpreted = _interpreted_request()
    diagnostics = _required_protein_aware_diagnostics(_protein_aware_result())
    fingerprints = build_protein_aware_policy_provenance(
        request=interpreted,
        resolved_inputs=resolved_inputs,
        protein_aware_diagnostics=diagnostics,
    ).input_fingerprints
    return {
        "phospho_matrix": fingerprints.phospho_matrix.exact_hash_value,
        "protein_matched_pairs": (fingerprints.protein_matched_pairs.exact_hash_value),
        "protein_covariate_matrix": (
            fingerprints.protein_covariate_matrix.exact_hash_value
        ),
        "protein_site_eligibility": (
            fingerprints.protein_site_eligibility.exact_hash_value
        ),
        "design_matrix": fingerprints.design_matrix.exact_hash_value,
        "contrast_matrix": fingerprints.contrast_matrix.exact_hash_value,
    }


def _resolved_inputs_with_fingerprint_mutation(
    resolved_inputs: ProteinAwareDifferentialResolvedInputs,
    *,
    mutation: str,
) -> ProteinAwareDifferentialResolvedInputs:
    sample_order = resolved_inputs.sample_order
    phosphosite_matrix = pd.DataFrame(
        resolved_inputs.computation_request.phosphosite_matrix,
        copy=True,
    )
    base_design = resolved_inputs.base_design
    base_contrasts = resolved_inputs.base_contrasts
    matched_pairs = pd.DataFrame(resolved_inputs.matched_pairs, copy=True)
    resolved_protein_covariates = pd.DataFrame(
        resolved_inputs.resolved_protein_covariates,
        copy=True,
    )
    site_eligibility_metadata = pd.DataFrame(
        resolved_inputs.site_eligibility_metadata,
        copy=True,
    )
    feature_eligibility_inputs = resolved_inputs.feature_eligibility_inputs

    if mutation == "phospho_matrix":
        phosphosite_matrix.iat[0, 0] = (
            float(cast(Any, phosphosite_matrix.iat[0, 0])) + 0.25
        )
    elif mutation == "protein_matched_pairs":
        matched_pairs.loc[0, "protein_identifier"] = "MAPK14_ALT"
    elif mutation == "protein_covariate_matrix":
        resolved_protein_covariates.iat[0, 0] = (
            float(cast(Any, resolved_protein_covariates.iat[0, 0])) + 0.25
        )
    elif mutation == "protein_site_eligibility":
        site_eligibility_metadata.loc[
            resolved_inputs.full_site_ids[1],
            "protein_aware_failure_message",
        ] = "changed protein-aware preparation failure message"
        feature_eligibility_inputs = DifferentialFeatureEligibilityInputs(
            feature_metadata=site_eligibility_metadata,
            result_status=site_eligibility_metadata[DIFFERENTIAL_RESULT_STATUS_COLUMN],
            testable_feature_ids=resolved_inputs.tested_site_ids,
            attach_to_result_tables=True,
        )
    elif mutation == "design_matrix":
        design_frame = base_design.to_dataframe()
        design_frame.iat[0, 0] = float(cast(Any, design_frame.iat[0, 0])) + 0.125
        base_design = DesignMatrix(design_frame)
    elif mutation == "contrast_matrix":
        contrast_frame = base_contrasts.to_dataframe()
        contrast_frame.iat[0, 0] = float(cast(Any, contrast_frame.iat[0, 0])) + 0.125
        base_contrasts = ContrastMatrix(contrast_frame)
    elif mutation == "sample_order":
        sample_order = tuple(reversed(sample_order))
        phosphosite_matrix = phosphosite_matrix.loc[:, list(sample_order)]
        resolved_protein_covariates = resolved_protein_covariates.loc[
            :,
            list(sample_order),
        ]
        base_design = DesignMatrix(base_design.to_dataframe().loc[list(sample_order)])
    else:
        raise AssertionError(f"unsupported fingerprint mutation: {mutation}")

    computation_request = ProteinAwareDifferentialComputationRequest(
        phosphosite_matrix=phosphosite_matrix,
        base_design=base_design,
        base_contrasts=base_contrasts,
        sample_order=sample_order,
        matched_pairs=matched_pairs,
        resolved_protein_covariates=resolved_protein_covariates,
        empirical_bayes=resolved_inputs.computation_request.empirical_bayes,
        multiple_testing_method=(
            resolved_inputs.computation_request.multiple_testing_method
        ),
        method_id=resolved_inputs.method_id,
    )
    return replace(
        resolved_inputs,
        computation_request=computation_request,
        feature_eligibility_inputs=feature_eligibility_inputs,
        matched_pairs=matched_pairs,
        resolved_protein_covariates=resolved_protein_covariates,
        site_eligibility_metadata=site_eligibility_metadata,
        sample_order=sample_order,
        base_design=base_design,
        base_contrasts=base_contrasts,
        status_counts=_status_counts(site_eligibility_metadata),
        reason_counts=_reason_counts(site_eligibility_metadata),
    )


def _computation_result(
    interpreted: InterpretedDifferentialAnalysisRequest,
    resolved_inputs: ProteinAwareDifferentialResolvedInputs,
    *,
    protein_coefficient: float,
    failure_site_id: str | None = None,
    failure_total_protein_row_key: str | None = None,
) -> ProteinAwareDifferentialComputationResult:
    failure_site_ids = () if failure_site_id is None else (failure_site_id,)
    tested_index = pd.Index(
        tuple(
            site_id
            for site_id in resolved_inputs.tested_site_ids
            if site_id not in failure_site_ids
        ),
        name="site_key",
    )
    prior_diagnostics = EmpiricalBayesPriorDiagnostics(
        method="standard",
        robust=False,
        trend=False,
        winsor_tail_p=(0.05, 0.1),
        base_prior_variance=0.5,
        base_prior_degrees_of_freedom=4.0,
        robust_outlier_count=0,
        robust_outlier_fraction=0.0,
        winsorized_low_count=0,
        winsorized_high_count=0,
        prior_variance=pd.Series([0.5], index=tested_index, name="prior_variance"),
        prior_degrees_of_freedom=pd.Series(
            [4.0],
            index=tested_index,
            name="prior_degrees_of_freedom",
        ),
    )
    return ProteinAwareDifferentialComputationResult(
        residual_variance=pd.Series([0.42], index=tested_index, name="s2"),
        posterior_residual_variance=pd.Series(
            [0.45],
            index=tested_index,
            name="s2_post",
        ),
        prior_residual_variance=pd.Series(
            [0.5],
            index=tested_index,
            name="s2_prior",
        ),
        prior_degrees_of_freedom_series_value=pd.Series(
            [4.0],
            index=tested_index,
            name="df_prior",
        ),
        prior_variance=0.5,
        prior_degrees_of_freedom=4.0,
        residual_degrees_of_freedom=3.0,
        empirical_bayes_method="standard",
        empirical_bayes_robust=False,
        empirical_bayes_trend=False,
        prior_diagnostics=prior_diagnostics,
        mean_variance_trend_diagnostics=None,
        contrast_tables={
            "B_vs_A": pd.DataFrame(
                {
                    "logFC": [1.25],
                    "t": [2.5],
                    "P.Value": [0.04],
                    "adj.P.Val": [0.04],
                },
                index=tested_index,
            )
        },
        protein_coefficient=pd.Series(
            [protein_coefficient],
            index=tested_index,
            name="protein_covariate_coefficient",
        ),
        site_diagnostics=pd.DataFrame(
            {
                "site_key": tested_index.tolist(),
                "protein_identifier": ["MAPK14"],
                "total_protein_row_key": ["protein:MAPK14"],
                "protein_coefficient": [protein_coefficient],
                "mean_intensity": [12.0],
            },
            index=tested_index,
        ),
        augmented_design_diagnostics=pd.DataFrame(
            {
                "total_protein_row_key": ["protein:MAPK14"],
                "sample_count": [6],
                "coefficient_count": [3],
                "rank": [3],
                "residual_degrees_of_freedom": [3.0],
                "condition_number": [11.0],
                "max_condition_number": [1.0e10],
                "decomposition_method": [
                    interpreted.design_decomposition.decomposition_method
                ],
                "solver": [interpreted.design_decomposition.solver],
                "column_scale_method": [
                    interpreted.design_decomposition.column_scale_method
                ],
                "rank_tolerance_policy": [
                    interpreted.design_decomposition.rank_tolerance_policy
                ],
                "rank_tolerance": [interpreted.design_decomposition.rank_tolerance],
            },
            index=pd.Index(["protein:MAPK14"], name="total_protein_row_key"),
        ),
        site_failure_diagnostics=_site_failure_diagnostics(
            failure_site_id=failure_site_id,
            total_protein_row_key=failure_total_protein_row_key,
        ),
        augmented_design_failure_diagnostics=_augmented_design_failure_diagnostics(
            total_protein_row_key=failure_total_protein_row_key,
        ),
        tested_site_ids=tested_index.tolist(),
        method_id=(
            DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
        ),
    )


def _site_failure_diagnostics(
    *,
    failure_site_id: str | None,
    total_protein_row_key: str | None,
) -> pd.DataFrame | None:
    if failure_site_id is None or total_protein_row_key is None:
        return None
    return pd.DataFrame(
        {
            "site_key": [failure_site_id],
            "total_protein_row_key": [total_protein_row_key],
            DIFFERENTIAL_RESULT_STATUS_COLUMN: [
                DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID
            ],
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN: [
                DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT
            ],
            "failure_message": ["protein augmented design is rank deficient"],
        },
        index=pd.Index([failure_site_id], name="site_key"),
    )


def _augmented_design_failure_diagnostics(
    *,
    total_protein_row_key: str | None,
) -> pd.DataFrame | None:
    if total_protein_row_key is None:
        return None
    return pd.DataFrame(
        {
            "total_protein_row_key": [total_protein_row_key],
            DIFFERENTIAL_RESULT_STATUS_COLUMN: [
                DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID
            ],
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN: [
                DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT
            ],
            "sample_count": [6],
            "coefficient_count": [3],
            "rank": [2],
            "residual_degrees_of_freedom": [4.0],
            "condition_number": [1.0e11],
            "max_condition_number": [1.0e10],
            "protein_raw_mean": [12.9],
            "protein_raw_standard_deviation": [0.72],
            "protein_centered_variance": [0.52],
            "failure_message": ["protein augmented design is rank deficient"],
        },
        index=pd.Index([total_protein_row_key], name="total_protein_row_key"),
    )


def _protein_aware_feature_metadata(
    site_ids: tuple[str, ...],
    *,
    protein_coefficient: float,
) -> pd.DataFrame:
    statuses = [
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID,
        DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE,
    ]
    reasons = [
        "",
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE,
    ]
    centered_variance = [0.53, np.nan, np.nan, 0.0, 0.70, 0.65]
    return pd.DataFrame(
        {
            "site_key": list(site_ids),
            "protein_aware_method_id": [
                DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
            ]
            * len(site_ids),
            "protein_aware_centering_policy": ["mean_centered_no_standardization"]
            * len(site_ids),
            "protein_aware_candidate": [True, False, False, True, True, True],
            "protein_aware_tested": [True, False, False, False, False, False],
            "protein_aware_preparation_eligibility": [
                "eligible_for_protein_aware_preparation",
                "fallback_to_phospho_only",
                "excluded_from_protein_aware_preparation",
                "eligible_for_protein_aware_preparation",
                "eligible_for_protein_aware_preparation",
                "eligible_for_protein_aware_preparation",
            ],
            "protein_aware_preparation_reasons": [
                ("matched_protein_available",),
                ("missing_total_protein_row",),
                ("ambiguous_protein_mapping",),
                ("matched_protein_available",),
                ("matched_protein_available",),
                ("matched_protein_available",),
            ],
            "protein_aware_mapping_status": [
                "matched",
                "missing_total_protein_row",
                "ambiguous",
                "matched",
                "matched",
                "matched",
            ],
            "protein_identifier": list(_PROTEINS),
            "total_protein_row_key": list(_TOTAL_ROW_KEYS),
            "protein_covariate_raw_mean": [10.8, np.nan, np.nan, 10.0, 11.0, 12.0],
            "protein_covariate_raw_standard_deviation": [
                0.72,
                np.nan,
                np.nan,
                0.0,
                0.84,
                0.81,
            ],
            "protein_covariate_centered_variance": centered_variance,
            "protein_augmented_design_rank": [3, np.nan, np.nan, 3, 2, 3],
            "protein_augmented_design_residual_degrees_of_freedom": [
                3.0,
                np.nan,
                np.nan,
                3.0,
                4.0,
                3.0,
            ],
            "protein_augmented_design_condition_number": [
                11.0,
                np.nan,
                np.nan,
                8.0,
                1.0e11,
                9.5,
            ],
            "protein_augmented_design_max_condition_number": [1.0e10] * len(site_ids),
            "protein_aware_failure_message": [
                "",
                "protein-aware preparation reported fallback",
                "protein-aware preparation excluded this phosphosite",
                "protein covariate has zero centered variance",
                "protein augmented design is rank deficient",
                "protein-aware contrast is not estimable",
            ],
            DIFFERENTIAL_RESULT_STATUS_COLUMN: statuses,
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN: reasons,
            "protein_covariate_coefficient": [
                protein_coefficient,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ],
        },
        index=pd.Index(site_ids, name="site_key"),
    )


def _required_protein_aware_diagnostics(
    result: DifferentialAnalysisResult,
) -> ProteinAwareDifferentialDiagnostics:
    diagnostics = result.protein_aware_diagnostics
    assert diagnostics is not None
    return diagnostics


def _required_protein_aware_policy(
    result: DifferentialAnalysisResult | None = None,
) -> DifferentialProteinAwarePolicyProvenance:
    resolved_result = _protein_aware_result() if result is None else result
    policy = resolved_result.policy_provenance
    assert policy is not None
    protein_aware = policy.protein_aware
    assert protein_aware is not None
    return protein_aware


def _protein_aware_policy_with_override(
    policy: DifferentialProteinAwarePolicyProvenance,
    *,
    field_name: str,
    value: object,
) -> DifferentialProteinAwarePolicyProvenance:
    if field_name == "claim_status":
        return replace(policy, claim_status=cast(str, value))
    if field_name == "model_formula":
        return replace(policy, model_formula=cast(str, value))
    if field_name == "nuisance_coefficient_name":
        return replace(policy, nuisance_coefficient_name=cast(str, value))
    if field_name == "protein_covariate_centered":
        return replace(policy, protein_covariate_centered=cast(bool, value))
    if field_name == "protein_covariate_standardized":
        return replace(policy, protein_covariate_standardized=cast(bool, value))
    if field_name == "automatic_protein_normalization":
        return replace(policy, automatic_protein_normalization=cast(bool, value))
    if field_name == "protein_imputation":
        return replace(policy, protein_imputation=cast(bool, value))
    if field_name == "phosphosite_only_fallback":
        return replace(policy, phosphosite_only_fallback=cast(bool, value))
    raise AssertionError(f"unsupported policy override: {field_name}")


def _protein_aware_diagnostics_copy(
    diagnostics: ProteinAwareDifferentialDiagnostics,
    *,
    method_id: str | None = None,
    claim_status: str | None = None,
    protein_covariate_centered: bool | None = None,
    protein_covariate_standardized: bool | None = None,
    protein_covariate_centering_policy: str | None = None,
    protein_covariate_imputation_policy: str | None = None,
    fallback_policy: str | None = None,
    minimum_condition_replicates: int | None = None,
    per_site_diagnostics: pd.DataFrame | None = None,
) -> ProteinAwareDifferentialDiagnostics:
    return ProteinAwareDifferentialDiagnostics(
        method_id=diagnostics.method_id if method_id is None else method_id,
        claim_status=diagnostics.claim_status if claim_status is None else claim_status,
        preparation_policy=diagnostics.preparation_policy,
        protein_mapping_policy=diagnostics.protein_mapping_policy,
        protein_covariate_centered=(
            diagnostics.protein_covariate_centered
            if protein_covariate_centered is None
            else protein_covariate_centered
        ),
        protein_covariate_standardized=(
            diagnostics.protein_covariate_standardized
            if protein_covariate_standardized is None
            else protein_covariate_standardized
        ),
        protein_covariate_centering_policy=(
            diagnostics.protein_covariate_centering_policy
            if protein_covariate_centering_policy is None
            else protein_covariate_centering_policy
        ),
        protein_covariate_imputation_policy=(
            diagnostics.protein_covariate_imputation_policy
            if protein_covariate_imputation_policy is None
            else protein_covariate_imputation_policy
        ),
        fallback_policy=diagnostics.fallback_policy
        if fallback_policy is None
        else fallback_policy,
        minimum_condition_replicates=(
            diagnostics.minimum_condition_replicates
            if minimum_condition_replicates is None
            else minimum_condition_replicates
        ),
        total_site_count=diagnostics.total_site_count,
        ordinary_eligible_site_count=diagnostics.ordinary_eligible_site_count,
        protein_preparation_eligible_site_count=(
            diagnostics.protein_preparation_eligible_site_count
        ),
        tested_site_count=diagnostics.tested_site_count,
        withheld_site_count=diagnostics.withheld_site_count,
        status_counts=diagnostics.status_counts,
        reason_counts=diagnostics.reason_counts,
        execution_sample_order=diagnostics.execution_sample_order,
        base_design_rank=diagnostics.base_design_rank,
        base_residual_degrees_of_freedom=(diagnostics.base_residual_degrees_of_freedom),
        expected_augmented_rank=diagnostics.expected_augmented_rank,
        common_augmented_rank=diagnostics.common_augmented_rank,
        common_augmented_residual_degrees_of_freedom=(
            diagnostics.common_augmented_residual_degrees_of_freedom
        ),
        distinct_matched_protein_row_count=(
            diagnostics.distinct_matched_protein_row_count
        ),
        fitted_protein_row_count=diagnostics.fitted_protein_row_count,
        condition_number_summary_scope=diagnostics.condition_number_summary_scope,
        condition_number_summary_statistic=(
            diagnostics.condition_number_summary_statistic
        ),
        min_augmented_condition_number=diagnostics.min_augmented_condition_number,
        median_augmented_condition_number=(
            diagnostics.median_augmented_condition_number
        ),
        max_augmented_condition_number=diagnostics.max_augmented_condition_number,
        per_site_diagnostics=(
            diagnostics.per_site_diagnostics_dataframe()
            if per_site_diagnostics is None
            else per_site_diagnostics
        ),
    )


def _protein_aware_diagnostics_copy_with_override(
    diagnostics: ProteinAwareDifferentialDiagnostics,
    *,
    field_name: str,
    value: object,
) -> ProteinAwareDifferentialDiagnostics:
    if field_name == "claim_status":
        return _protein_aware_diagnostics_copy(
            diagnostics,
            claim_status=cast(str, value),
        )
    if field_name == "method_id":
        return _protein_aware_diagnostics_copy(
            diagnostics,
            method_id=cast(str, value),
        )
    if field_name == "protein_covariate_centered":
        return _protein_aware_diagnostics_copy(
            diagnostics,
            protein_covariate_centered=cast(bool, value),
        )
    if field_name == "protein_covariate_standardized":
        return _protein_aware_diagnostics_copy(
            diagnostics,
            protein_covariate_standardized=cast(bool, value),
        )
    if field_name == "protein_covariate_centering_policy":
        return _protein_aware_diagnostics_copy(
            diagnostics,
            protein_covariate_centering_policy=cast(str, value),
        )
    if field_name == "protein_covariate_imputation_policy":
        return _protein_aware_diagnostics_copy(
            diagnostics,
            protein_covariate_imputation_policy=cast(str, value),
        )
    if field_name == "fallback_policy":
        return _protein_aware_diagnostics_copy(
            diagnostics,
            fallback_policy=cast(str, value),
        )
    if field_name == "minimum_condition_replicates":
        return _protein_aware_diagnostics_copy(
            diagnostics,
            minimum_condition_replicates=cast(int, value),
        )
    raise AssertionError(f"unsupported diagnostics override: {field_name}")


def _result_with_protein_aware_diagnostics(
    result: DifferentialAnalysisResult,
    *,
    protein_aware_diagnostics: ProteinAwareDifferentialDiagnostics,
) -> DifferentialAnalysisResult:
    return DifferentialAnalysisResult.from_trusted_owned(
        residual_variance=result.residual_variance_series(),
        posterior_residual_variance=result.posterior_residual_variance_series(),
        prior_residual_variance=result.prior_residual_variance_series(),
        prior_degrees_of_freedom_series_value=result.prior_degrees_of_freedom_series(),
        prior_variance=result.prior_variance,
        prior_degrees_of_freedom=result.prior_degrees_of_freedom,
        residual_degrees_of_freedom=result.residual_degrees_of_freedom,
        empirical_bayes_method=result.empirical_bayes_method,
        empirical_bayes_robust=result.empirical_bayes_robust,
        empirical_bayes_trend=result.empirical_bayes_trend,
        prior_diagnostics=result.prior_diagnostics,
        mean_variance_trend_diagnostics=result.mean_variance_trend_diagnostics,
        contrast_tables=result.contrast_tables,
        diagnostics=result.diagnostics,
        policy_provenance=result.policy_provenance,
        workflow_provenance=result.workflow_provenance,
        caveats=result.caveats,
        input_dataset_preprocessing_report=result.input_dataset_preprocessing_report,
        feature_eligibility=result.feature_eligibility,
        protein_aware_diagnostics=protein_aware_diagnostics,
    )


def _resolved_inputs_with_request_mismatch(
    resolved_inputs: ProteinAwareDifferentialResolvedInputs,
    *,
    mismatch: str,
) -> ProteinAwareDifferentialResolvedInputs:
    sample_order = resolved_inputs.sample_order
    base_design = resolved_inputs.base_design
    base_contrasts = resolved_inputs.base_contrasts
    phosphosite_matrix = resolved_inputs.computation_request.phosphosite_matrix
    resolved_protein_covariates = resolved_inputs.resolved_protein_covariates
    if mismatch == "sample_order":
        sample_order = tuple(reversed(sample_order))
        base_design = DesignMatrix(base_design.frame.loc[list(sample_order)])
        phosphosite_matrix = phosphosite_matrix.loc[:, list(sample_order)]
        resolved_protein_covariates = resolved_protein_covariates.loc[
            :,
            list(sample_order),
        ]
    elif mismatch == "base_design":
        frame = base_design.to_dataframe()
        frame.iat[0, 0] = float(cast(Any, frame.iat[0, 0])) + 0.125
        base_design = DesignMatrix(frame)
    elif mismatch == "base_contrasts":
        frame = base_contrasts.to_dataframe()
        frame.iat[0, 0] = float(cast(Any, frame.iat[0, 0])) + 0.125
        base_contrasts = ContrastMatrix(frame)
    else:
        raise AssertionError(f"unsupported mismatch case: {mismatch}")
    computation_request = ProteinAwareDifferentialComputationRequest(
        phosphosite_matrix=phosphosite_matrix,
        base_design=base_design,
        base_contrasts=base_contrasts,
        sample_order=sample_order,
        matched_pairs=resolved_inputs.matched_pairs,
        resolved_protein_covariates=resolved_protein_covariates,
        empirical_bayes=resolved_inputs.computation_request.empirical_bayes,
        multiple_testing_method=(
            resolved_inputs.computation_request.multiple_testing_method
        ),
        method_id=resolved_inputs.method_id,
    )
    return replace(
        resolved_inputs,
        computation_request=computation_request,
        resolved_protein_covariates=resolved_protein_covariates,
        sample_order=sample_order,
        base_design=base_design,
        base_contrasts=base_contrasts,
    )


def _status_counts(feature_metadata: pd.DataFrame) -> tuple[tuple[str, int], ...]:
    values = [
        str(value) for value in feature_metadata[DIFFERENTIAL_RESULT_STATUS_COLUMN]
    ]
    return tuple((value, values.count(value)) for value in dict.fromkeys(values))


def _reason_counts(feature_metadata: pd.DataFrame) -> tuple[tuple[str, int], ...]:
    values = [
        str(value).strip()
        for value in feature_metadata[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN]
        if str(value).strip()
    ]
    return tuple((value, values.count(value)) for value in dict.fromkeys(values))


def _interpreted_request(
    *,
    config: DifferentialAnalysisConfig | None = None,
    include_blocks: bool = False,
) -> InterpretedDifferentialAnalysisRequest:
    return DifferentialAnalysisInterpreter().run(
        DifferentialAnalysisValidator().run(
            _request(config=config, include_blocks=include_blocks)
        )
    )


def _request(
    *,
    dataset: object | None = None,
    config: DifferentialAnalysisConfig | None = None,
    include_blocks: bool = False,
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=_dataset() if dataset is None else dataset,
        design=ExperimentalDesign(
            samples=_sample_design(include_blocks=include_blocks)
        ),
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
        ),
        config=DifferentialAnalysisConfig() if config is None else config,
    )


def _protein_aware_config() -> DifferentialAnalysisConfig:
    return DifferentialAnalysisConfig(
        protein_aware_model=DifferentialProteinAwareModelConfig(
            method=(
                DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1
            )
        )
    )


def _sample_design(*, include_blocks: bool) -> tuple[SampleDesignRecord, ...]:
    block_ids = ("donor_1", "donor_2", "donor_3", "donor_1", "donor_2", "donor_3")
    conditions = ("A", "A", "A", "B", "B", "B")
    records: list[SampleDesignRecord] = []
    for sample_id, condition, block_id in zip(
        _SAMPLE_IDS,
        conditions,
        block_ids,
        strict=True,
    ):
        records.append(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=condition,
                biological_replicate_id=block_id,
                block_id=block_id if include_blocks else None,
            )
        )
    return tuple(records)


def _dataset_with_trusted_protein_aware_preparation_provenance(
    *,
    include_reference_context: bool,
):
    base_dataset = _dataset()
    phospho = base_dataset.phospho
    site_metadata = base_dataset.site_metadata
    total = pd.DataFrame(
        [[10.0, 10.2, 10.4, 11.2, 11.4, 11.6]],
        index=pd.Index(["protein:MAPK14"], name="protein_id"),
        columns=pd.Index(_SAMPLE_IDS, name="sample_id"),
    )
    report = _trusted_protein_aware_preparation_report(
        phospho=phospho,
        site_metadata=site_metadata,
        total=total,
        include_reference_context=include_reference_context,
    )
    preparation = ProteinAwarePreparationResult(
        matched_pairs=pd.DataFrame(
            {
                "site_key": [str(phospho.index[0])],
                "protein_identifier": ["MAPK14"],
                "total_protein_row_key": ["protein:MAPK14"],
            }
        ),
        protein_covariate_matrix=total.loc[["protein:MAPK14"], :],
        report=report,
    )
    preprocessing_report = DatasetPreprocessingReport.from_rows(
        protein_aware_preparation=report
    )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        total=total,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=True
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=True),
        preprocessing_report=preprocessing_report,
        protein_aware_preparation=preparation,
    )


def _trusted_protein_aware_preparation_report(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    total: pd.DataFrame,
    include_reference_context: bool,
) -> ProteinAwarePreparationReport:
    site_keys = tuple(str(value) for value in phospho.index.tolist())
    site_eligibility: list[ProteinAwareSiteEligibility] = []
    for position, site_key in enumerate(site_keys):
        if position == 0:
            site_eligibility.append(
                ProteinAwareSiteEligibility(
                    site_key=site_key,
                    eligibility=(
                        ProteinAwarePreparationEligibility.ELIGIBLE_FOR_PROTEIN_AWARE_PREPARATION
                    ),
                    mapping_status=ProteinMappingStatus.MATCHED,
                    protein_identifier="MAPK14",
                    total_protein_row_key="protein:MAPK14",
                    reasons=("matched_protein_available",),
                )
            )
            continue
        site_eligibility.append(
            ProteinAwareSiteEligibility(
                site_key=site_key,
                eligibility=ProteinAwarePreparationEligibility.EXCLUDED_FROM_PREPARATION,
                mapping_status=ProteinMappingStatus.MISSING_TOTAL_PROTEIN_ROW,
                protein_identifier=str(
                    site_metadata.loc[site_key, "protein_identifier"]
                ),
                total_protein_row_key=None,
                reasons=("missing_total_protein_row",),
            )
        )

    return ProteinAwarePreparationReport(
        site_eligibility=tuple(site_eligibility),
        sample_alignment=ProteinAwareSampleAlignmentDiagnostics(
            phospho_sample_columns=tuple(str(value) for value in phospho.columns),
            total_protein_sample_columns=tuple(str(value) for value in total.columns),
            exact_sample_order_match=True,
            sample_order_compatible=True,
            reordered_sample_columns=False,
            allow_reordered_samples=False,
            missing_total_protein_samples=(),
            extra_total_protein_samples=(),
        ),
        transformation_state=ProteinAwareTransformationStateDiagnostics(
            compatible=True,
            phospho_transformation_state={
                "kind": "log2",
                "transformed": True,
                "established_by": "test.phospho",
            },
            total_protein_transformation_state={
                "kind": "log2",
                "transformed": True,
                "established_by": "test.total",
            },
        ),
        preparation_policy="prepare_model_inputs",
        protein_mapping_policy="require_unambiguous",
        policy_parameters={
            "allow_reordered_samples": False,
            "dataset_binding_table_fingerprints": [
                table_fingerprint_to_payload(fingerprint)
                for fingerprint in (
                    fingerprint_optional_table_strict(
                        phospho,
                        name="dataset.phospho",
                    ),
                    fingerprint_optional_table_strict(
                        site_metadata,
                        name="dataset.site_metadata",
                    ),
                    fingerprint_optional_table_strict(total, name="dataset.total"),
                )
                if fingerprint is not None
            ],
        },
        provenance=_trusted_preparation_run_provenance(
            include_reference_context=include_reference_context
        ),
    )


def _trusted_preparation_run_provenance(
    *,
    include_reference_context: bool,
) -> RunProvenance:
    reference_context = ReferenceContext(
        organism=Organism.RAT,
        protein_namespace="test-protein-id",
        source_name="unit-test-reference",
        source_version="2026.1",
        proteome_version="unit-test-proteome",
        reference_table_sha256="a" * 64,
    )
    reference_table = pd.DataFrame(
        {"protein_identifier": ["MAPK14"]},
        index=pd.Index(["MAPK14"], name="protein_identifier"),
    )
    reference = ReferenceProvenance(
        source_type="local_fixture",
        organism=Organism.RAT,
        bundle_id="unit-test-reference-bundle",
        table_fingerprints=(
            fingerprint_table_strict(
                reference_table,
                name="protein_aware_preparation.reference_table",
            ),
        ),
        source_name="unit-test-reference",
        source_version="2026.1",
        retrieved_at="2026-08-31",
        identifier_namespace="test-protein-id",
        sequence_window={"center_residue": "S/T/Y"},
        manifest={
            "organism": "rat",
            "source_name": "unit-test-reference",
            "source_version": "2026.1",
        },
        reference_context=reference_context,
    )
    return RunProvenance(
        environment=EnvironmentProvenance(
            package_name="phospy",
            package_version="test",
            python_version="3.test",
            dependency_versions={},
        ),
        input_tables=(),
        preprocessing_stages=(),
        reference=reference,
        workflow_name="protein_aware_preparation",
        workflow_parameters={"policy": "prepare_model_inputs"},
        random_state=None,
        random_seed_policy=None,
        output_tables=(),
        reference_context=reference_context if include_reference_context else None,
    )


def _dataset():
    site_index = _site_index()
    phospho = pd.DataFrame(
        [
            [10.0, 10.4, 10.8, 11.2, 11.4, 11.8],
            [9.5, 9.7, 9.9, 10.1, 10.4, 10.5],
            [8.0, 8.2, 8.5, 9.0, 9.2, 9.5],
            [12.0, 12.1, 12.3, 11.9, 12.2, 12.4],
            [7.0, 7.2, 7.4, 8.1, 8.3, 8.6],
            [13.0, 12.8, 13.2, 13.5, 13.7, 13.9],
        ],
        index=site_index,
        columns=pd.Index(_SAMPLE_IDS, name="sample_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.tolist(),
            "display_id": [
                f"{protein};{site};"
                for protein, site in zip(_PROTEINS, _SITES, strict=True)
            ],
            **site_key_context_columns(site_index),
            "gene_symbol": list(_PROTEINS),
            "site": list(_SITES),
            "site_sequence": [("A" * 15) + site[0] + ("A" * 15) for site in _SITES],
            "protein_id": list(_PROTEINS),
        },
        index=site_index.copy(),
    )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def _site_index() -> pd.Index:
    return protein_site_key_index(
        protein_identifiers=list(_PROTEINS),
        sites=list(_SITES),
    )
