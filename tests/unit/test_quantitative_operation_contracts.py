from __future__ import annotations

from typing import cast

import pytest

from phospy.errors.build import DatasetBuildError
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.intensity_scale_state import (
    intensity_scale_state_from_payload,
    intensity_scale_state_to_payload,
)
from phospy.provenance.models import TableFingerprint
from phospy.science.configs.preprocessing import (
    DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
)
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageExecution,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.science.datasets.preprocessing.policy_models import (
    IntensityTransformPolicy,
    MissingDataPolicy,
    NormalisationPolicy,
    TotalProteinCorrectionPolicy,
)
from phospy.science.datasets.preprocessing.stage_registry import (
    get_preprocessing_stage_metadata,
    list_registered_preprocessing_stages,
)
from phospy.science.datasets.preprocessing.state_builder import (
    DatasetProcessingStateBuilder,
)
from phospy.science.transformations._authority import (
    dataset_resolver_establishment_authority,
)
from phospy.science.transformations.models import (
    QUANTITATIVE_MEANING_OPERATION_TOTAL_PROTEIN_SUBTRACT_LOG_TOTAL,
    IntensityScaleEstablishmentMode,
    IntensityScaleEvidenceLevel,
    IntensityScaleKind,
    IntensityScaleState,
    MatrixIntensityScaleState,
    QuantitativeMeaning,
    establish_intensity_scale_state,
)
from phospy.science.transformations.quantitative_contracts import (
    ALL_INTENSITY_SCALE_KINDS,
    ALL_QUANTITATIVE_MEANINGS,
    NegativeDomainPolicy,
    QuantitativeContractState,
    QuantitativeInformationLossKind,
    QuantitativeOperationContract,
    QuantitativeReversibilityKind,
    QuantitativeScaleTransitionKind,
    QuantitativeTransitionEvidence,
    preserve_quantitative_contract,
)


def test_preserve_contract_is_exhaustive_over_scale_and_meaning_enums() -> None:
    contract = preserve_quantitative_contract()

    assert ALL_INTENSITY_SCALE_KINDS == frozenset(IntensityScaleKind)
    assert ALL_QUANTITATIVE_MEANINGS == frozenset(QuantitativeMeaning)
    assert contract.accepted_input_scale_kinds == frozenset(IntensityScaleKind)
    assert contract.accepted_quantitative_meanings == frozenset(QuantitativeMeaning)
    assert frozenset(
        contract.output_scale_transition.output_scale_by_input
    ) == frozenset(IntensityScaleKind)
    assert frozenset(
        contract.output_meaning_transition.output_meaning_by_input
    ) == frozenset(QuantitativeMeaning)


def test_built_in_quantitative_contracts_cover_all_policy_branches() -> None:
    representative_plans = (
        PreprocessingPlan.default(),
        PreprocessingPlan(
            intensity_transform_policy=IntensityTransformPolicy.LOG2,
            stage_order=("intensity_transform",),
        ),
        PreprocessingPlan(
            missing_data_policy=MissingDataPolicy.IMPUTE_MINPROB,
            stage_order=("missing_data",),
        ),
        PreprocessingPlan(
            missing_data_policy=MissingDataPolicy.IMPUTE_ROW_MEDIAN,
            stage_order=("missing_data",),
        ),
        PreprocessingPlan(
            missing_data_policy=MissingDataPolicy.IMPUTE_KNN,
            stage_order=("missing_data",),
        ),
        PreprocessingPlan(
            normalisation_policy=NormalisationPolicy.MEDIAN_CENTER,
            stage_order=("normalisation",),
        ),
        PreprocessingPlan(
            normalisation_policy=NormalisationPolicy.QUANTILE,
            stage_order=("normalisation",),
        ),
        PreprocessingPlan(
            total_protein_correction_policy=(
                TotalProteinCorrectionPolicy.SUBTRACT_LOG_TOTAL
            ),
            stage_order=("total_protein_correction",),
        ),
        PreprocessingPlan(
            batch_correction_method=(
                DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH
            ),
            stage_order=("batch_correction",),
        ),
    )

    for metadata in list_registered_preprocessing_stages():
        observed_operations: set[str] = set()
        for plan in representative_plans:
            contract = metadata.resolve_quantitative_contract(plan)
            assert isinstance(contract, QuantitativeOperationContract)
            assert contract.accepted_input_scale_kinds
            assert contract.accepted_quantitative_meanings
            assert contract.required_evidence
            assert contract.negative_domain_policy in NegativeDomainPolicy
            assert contract.reversibility in QuantitativeReversibilityKind
            assert contract.information_loss in QuantitativeInformationLossKind
            assert (
                frozenset(contract.output_scale_transition.output_scale_by_input)
                >= contract.accepted_input_scale_kinds
            )
            assert (
                frozenset(contract.output_meaning_transition.output_meaning_by_input)
                >= contract.accepted_quantitative_meanings
            )
            observed_operations.add(metadata.operation_name(plan))
        assert observed_operations


def test_intensity_log2_contract_transition_matrix() -> None:
    plan = PreprocessingPlan(
        intensity_transform_policy=IntensityTransformPolicy.LOG2,
        stage_order=("intensity_transform",),
    )
    contract = get_preprocessing_stage_metadata(
        "intensity_transform"
    ).resolve_quantitative_contract(plan)

    abundance_output = contract.validate_and_transition(
        QuantitativeContractState(
            scale_kind=IntensityScaleKind.LINEAR,
            meaning=QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        ),
        stage="intensity_transform",
        operation="log2",
    )
    unknown_output = contract.validate_and_transition(
        QuantitativeContractState(
            scale_kind=IntensityScaleKind.LINEAR,
            meaning=QuantitativeMeaning.UNKNOWN,
        ),
        stage="intensity_transform",
        operation="log2",
    )

    assert contract.output_scale_transition.kind is (
        QuantitativeScaleTransitionKind.LINEAR_TO_LOG2
    )
    assert abundance_output == QuantitativeContractState(
        scale_kind=IntensityScaleKind.LOG2,
        meaning=QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
    )
    assert unknown_output == QuantitativeContractState(
        scale_kind=IntensityScaleKind.LOG2,
        meaning=QuantitativeMeaning.UNKNOWN,
    )
    with pytest.raises(PhosPyInputError, match="unsupported input scale"):
        contract.validate_and_transition(
            QuantitativeContractState(
                scale_kind=IntensityScaleKind.LOG2,
                meaning=QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
            ),
            stage="intensity_transform",
            operation="log2",
        )


def test_total_protein_contract_transition_matrix_uses_typed_evidence() -> None:
    plan = PreprocessingPlan(
        total_protein_correction_policy=TotalProteinCorrectionPolicy.SUBTRACT_LOG_TOTAL,
        stage_order=("total_protein_correction",),
    )
    contract = get_preprocessing_stage_metadata(
        "total_protein_correction"
    ).resolve_quantitative_contract(plan)
    input_state = QuantitativeContractState(
        scale_kind=IntensityScaleKind.LOG2,
        meaning=QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
    )

    fully_corrected = contract.validate_and_transition(
        input_state,
        stage="total_protein_correction",
        operation="subtract_log_total",
        evidence=QuantitativeTransitionEvidence(
            total_protein_corrected_row_count=2,
            total_protein_uncorrected_row_count=0,
        ),
    )
    mixed = contract.validate_and_transition(
        input_state,
        stage="total_protein_correction",
        operation="subtract_log_total",
        evidence=QuantitativeTransitionEvidence(
            total_protein_corrected_row_count=1,
            total_protein_uncorrected_row_count=1,
        ),
    )

    assert fully_corrected == QuantitativeContractState(
        scale_kind=IntensityScaleKind.LOG2,
        meaning=QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
    )
    assert mixed == QuantitativeContractState(
        scale_kind=IntensityScaleKind.LOG2,
        meaning=(
            QuantitativeMeaning.MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE
        ),
    )
    assert "quantitative_meaning_mixed_total_protein_correction" in (
        contract.caveat_codes(
            target=mixed.meaning,
            evidence=QuantitativeTransitionEvidence(
                total_protein_corrected_row_count=1,
                total_protein_uncorrected_row_count=1,
            ),
        )
    )


def test_pipeline_rejects_unsupported_transition_before_numerical_execution() -> None:
    class ExplodingNormalisationStage:
        stage_key = DATASET_PREPROCESSING_STAGE_NORMALISATION
        ran = False

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            self.ran = True
            raise AssertionError("stage numerical execution should not run")

    stage = ExplodingNormalisationStage()
    state = PreprocessingState(
        phospho=None,  # type: ignore[arg-type]
        site_metadata=None,  # type: ignore[arg-type]
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            normalisation_policy=NormalisationPolicy.MEDIAN_CENTER,
            stage_order=(DATASET_PREPROCESSING_STAGE_NORMALISATION,),
        ),
    )

    with pytest.raises(PhosPyInputError, match="unsupported input scale"):
        PreprocessingPipeline(stage_registry=(stage,)).run_with_trace(state)
    assert stage.ran is False


def test_pipeline_rejects_unsupported_meaning_before_numerical_execution() -> None:
    class ExplodingTotalProteinStage:
        stage_key = DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION
        ran = False

        def run(self, state: PreprocessingState) -> PreprocessingStageResult:
            self.ran = True
            raise AssertionError("stage numerical execution should not run")

    stage = ExplodingTotalProteinStage()
    state = PreprocessingState(
        phospho=None,  # type: ignore[arg-type]
        site_metadata=None,  # type: ignore[arg-type]
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(
            total_protein_correction_policy=(
                TotalProteinCorrectionPolicy.SUBTRACT_LOG_TOTAL
            ),
            stage_order=(DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,),
        ),
    )

    with pytest.raises(
        PhosPyInputError, match="unsupported input quantitative meaning"
    ):
        PreprocessingPipeline(stage_registry=(stage,)).run_with_trace(
            state,
            initial_quantitative_scale_kind=IntensityScaleKind.LOG2,
            initial_quantitative_meaning=QuantitativeMeaning.UNKNOWN,
        )
    assert stage.ran is False


def test_state_builder_generates_quantitative_meaning_provenance_from_contract() -> (
    None
):
    processing_state = DatasetProcessingStateBuilder().build(
        plan=PreprocessingPlan(
            total_protein_correction_policy=(
                TotalProteinCorrectionPolicy.SUBTRACT_LOG_TOTAL
            ),
            stage_order=(DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,),
        ),
        intensity_scale_state=_established_log2_state(
            quantity=QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        ),
        preprocessing_trace=(
            _total_protein_stage_execution(
                diagnostics_quantitative_meaning="phospho_total_log_ratio",
                evidence=QuantitativeTransitionEvidence(
                    total_protein_corrected_row_count=2,
                    total_protein_uncorrected_row_count=0,
                ),
            ),
        ),
    )

    provenance = processing_state.intensity_scale.quantitative_meaning_provenance
    assert provenance is not None
    payload = provenance.to_payload()
    assert (
        payload["operation_id"]
        == QUANTITATIVE_MEANING_OPERATION_TOTAL_PROTEIN_SUBTRACT_LOG_TOTAL
    )
    assert payload["source_quantity"] == "phosphosite_log_abundance"
    assert payload["target_quantity"] == "phospho_total_log_ratio"
    parameters = cast(dict[str, object], payload["parameters"])
    assert parameters["evidence.total_protein_corrected_row_count"] == 2
    assert parameters["evidence.total_protein_uncorrected_row_count"] == 0
    semantic_contract = cast(dict[str, object], parameters["semantic_contract"])
    assert semantic_contract["negative_domain_policy"] == "allows_negative_output"
    assert semantic_contract["information_loss"] == "ratio_transformation"
    assert semantic_contract["required_evidence"] == [
        "table_fingerprints",
        "total_protein_row_mapping",
    ]

    state_payload = intensity_scale_state_to_payload(processing_state.intensity_scale)
    restored = intensity_scale_state_from_payload(state_payload)
    assert restored.quantitative_meaning_provenance is not None
    assert restored.quantitative_meaning_provenance.to_payload() == payload


def test_state_builder_rejects_diagnostic_meaning_that_conflicts_with_contract() -> (
    None
):
    with pytest.raises(
        DatasetBuildError,
        match="diagnostics quantitative_meaning does not match the typed",
    ):
        DatasetProcessingStateBuilder().build(
            plan=PreprocessingPlan(
                total_protein_correction_policy=(
                    TotalProteinCorrectionPolicy.SUBTRACT_LOG_TOTAL
                ),
                stage_order=(DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,),
            ),
            intensity_scale_state=_established_log2_state(
                quantity=QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
            ),
            preprocessing_trace=(
                _total_protein_stage_execution(
                    diagnostics_quantitative_meaning="phospho_total_log_ratio",
                    evidence=QuantitativeTransitionEvidence(
                        total_protein_corrected_row_count=1,
                        total_protein_uncorrected_row_count=1,
                    ),
                ),
            ),
        )


def _established_log2_state(*, quantity: QuantitativeMeaning) -> IntensityScaleState:
    return establish_intensity_scale_state(
        IntensityScaleState(
            phospho=MatrixIntensityScaleState.log2(
                established_by="tests.quantitative_operation_contracts"
            ),
            total=MatrixIntensityScaleState.log2(
                established_by="tests.quantitative_operation_contracts"
            ),
            quantity=quantity,
        ),
        established_via="tests.quantitative_operation_contracts",
        establishment_mode=IntensityScaleEstablishmentMode.DERIVED,
        evidence_level=IntensityScaleEvidenceLevel.UNKNOWN,
        _authority=dataset_resolver_establishment_authority(),
    )


def _total_protein_stage_execution(
    *,
    diagnostics_quantitative_meaning: str,
    evidence: QuantitativeTransitionEvidence,
) -> PreprocessingStageExecution:
    return PreprocessingStageExecution(
        stage=DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
        operation="subtract_log_total",
        parameters={"total_protein_correction_policy": "subtract_log_total"},
        input_shape=(2, 2),
        output_shape=(2, 2),
        input_hash="input_hash",
        output_hash="output_hash",
        consumed_input_tables=(
            _fingerprint("dataset.phospho"),
            _fingerprint("dataset.total"),
            _fingerprint("dataset.site_metadata"),
        ),
        produced_output_tables=(_fingerprint("dataset.phospho"),),
        diagnostics={
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "requested_policy": "subtract_log_total",
            "resolved_policy": "subtract_log_total",
            "formula": "log2_phospho - log2_total",
            "requires_log_scale": True,
            "input_scale": "log2",
            "output_scale": "log2_ratio",
            "quantitative_meaning": diagnostics_quantitative_meaning,
            "matched_rows": 2,
            "corrected_row_count": evidence.total_protein_corrected_row_count,
            "uncorrected_row_count": evidence.total_protein_uncorrected_row_count,
        },
        quantitative_transition_evidence=evidence,
    )


def _fingerprint(name: str) -> TableFingerprint:
    return TableFingerprint(
        name=name,
        rows=2,
        columns=2,
        index_name="id",
        column_names=("sample_a", "sample_b"),
        dtypes=("float64", "float64"),
        exact_hash_algorithm="sha256",
        exact_hash_value=f"exact_{name}",
        tolerance_hash_algorithm="sha256",
        tolerance_hash_value=f"tolerance_{name}",
    )
