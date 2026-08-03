from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, cast

import pandas as pd
import pytest

from phospy.advanced.configs import SpsRuvBatchCorrectionConfig
from phospy.api.configs import (
    DatasetPreprocessingConfig,
)
from phospy.errors.build import DatasetBuildError
from phospy.errors.transformations import InvalidTransformationStateError
from phospy.provenance.models import BatchCorrectionProvenance
from phospy.science.configs import (
    DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
)
from phospy.science.configs.preprocessing import CorrectionMissingnessPolicy
from phospy.science.datasets.preprocessing.batch_correction import (
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
)
from phospy.science.datasets.preprocessing.batch_correction_provenance import (
    build_native_batch_correction_provenance,
)
from phospy.science.datasets.preprocessing.control_sites import (
    ControlSiteSet,
    ControlSiteSourceMetadata,
)
from phospy.science.datasets.preprocessing.correction_output import (
    CorrectedPreprocessingOutput,
)
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_LOCALISATION,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
    PreprocessingStateTableKey,
)
from phospy.science.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.science.datasets.preprocessing.policy_models import (
    ComparisonBuildingPolicy,
    ImputationInputScale,
    IntensityTransformPolicy,
    LocalisationEligibilityMode,
    MissingDataPolicy,
    SiteMatrixPolicy,
    TotalProteinCorrectionPolicy,
)
from phospy.science.datasets.preprocessing.quantitative_evidence import (
    QuantitativeOperationEvidence,
    QuantitativeOperationEvidenceValidator,
    RowAuditEvidence,
    TotalProteinRowMappingEvidence,
)
from phospy.science.datasets.preprocessing.stage_registry import (
    PreprocessingStageMetadata,
    list_registered_preprocessing_stages,
    resolve_registered_preprocessing_stages,
)
from phospy.science.datasets.preprocessing.state_builder import (
    DatasetProcessingStateBuilder,
)
from phospy.science.transformations._authority import (
    dataset_resolver_establishment_authority,
)
from phospy.science.transformations.models import (
    IntensityScaleEstablishmentMode,
    IntensityScaleEvidenceLevel,
    IntensityScaleKind,
    IntensityScaleState,
    MatrixIntensityScaleState,
    QuantitativeMeaning,
    establish_intensity_scale_state,
)
from phospy.science.transformations.quantitative_contracts import (
    NegativeDomainPolicy,
    QuantitativeEvidenceRequirement,
    QuantitativeInformationLossKind,
    QuantitativeOperationContract,
    QuantitativeReversibilityKind,
    QuantitativeTransitionEvidence,
    preserve_meaning_transition,
    preserve_scale_transition,
)
from tests.support.site_keys import site_key_index_from_display_ids


class _IdentityStage:
    def __init__(
        self,
        stage_key: str,
        *,
        result_evidence: QuantitativeOperationEvidence | None = None,
        transition_evidence: QuantitativeTransitionEvidence | None = None,
        diagnostics: dict[str, object] | None = None,
        add_observation_mask: bool = False,
    ) -> None:
        self.stage_key = stage_key
        self._result_evidence = result_evidence
        self._transition_evidence = transition_evidence
        self._diagnostics = diagnostics or {}
        self._add_observation_mask = add_observation_mask

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        next_state = state
        if self._add_observation_mask:
            next_state = replace(
                state,
                imputation_observation_mask=pd.DataFrame(
                    True,
                    index=state.phospho.index.copy(),
                    columns=state.phospho.columns.copy(),
                ),
            )
        return PreprocessingStageResult(
            state=next_state,
            diagnostics=self._diagnostics,
            quantitative_transition_evidence=self._transition_evidence,
            quantitative_evidence=self._result_evidence,
        )


@dataclass(frozen=True, slots=True)
class _FakeSpsRuvResult:
    corrected_preprocessing_output: CorrectedPreprocessingOutput
    diagnostics: dict[str, object]
    provenance: BatchCorrectionProvenance


class _NoopDesignValidator:
    def run(self, **_: object) -> object:
        return object()


class _NoopAdequacyValidator:
    def run(self, **_: object) -> object:
        return object()


class _FakeSpsRuvRunner:
    def run(
        self,
        *,
        phospho: pd.DataFrame,
        config: object,
        sample_metadata: pd.DataFrame | None,
        control_site_set: ControlSiteSet,
        missingness_policy: object,
        upstream_observation_mask: pd.DataFrame | None,
        site_metadata: pd.DataFrame,
    ) -> _FakeSpsRuvResult:
        del config, control_site_set, missingness_policy, site_metadata
        mask = (
            phospho.notna()
            if upstream_observation_mask is None
            else upstream_observation_mask.copy(deep=True)
        )
        corrected = phospho.copy(deep=True)
        report = BatchCorrectionReport(
            status="applied",
            policy=BatchCorrectionPolicy(
                method="sps_ruv_style",
                batch_column="batch",
                condition_column="condition",
            ),
            diagnostics=BatchCorrectionDiagnostics(
                number_of_batches=2,
                batch_levels=("batch_a", "batch_b"),
                condition_levels=("control", "treated"),
                matrix_shape_before=phospho.shape,
                matrix_shape_after=corrected.shape,
            ),
        )
        output = CorrectedPreprocessingOutput(
            corrected_matrix=corrected,
            output_observation_mask=mask.astype(bool),
            batch_correction_report=report,
            diagnostics={"status": "applied", "method": "sps_ruv_style"},
        )
        provenance = build_native_batch_correction_provenance(
            input_matrix=phospho,
            output_matrix=corrected,
            plan=_sps_ruv_plan(),
            report=report,
            metadata=None,
            diagnostics={"status": "applied", "method": "sps_ruv_style"},
            observation_mask=mask,
            control_site_source={
                "source_type": "caller_supplied",
                "organism": "rat",
                "identifier_namespace": "site_key",
                "source_name": "manual-controls",
                "source_version": "manual-v1",
                "license": "caller local use",
                "redistribution": "not redistributed",
                "selection_method": "caller_supplied",
            },
            selected_site_key_rows=(
                _site_key("MAPK14", "Y", "182"),
                _site_key("SRC", "Y", "416"),
            ),
            source="unit-test",
        )
        return _FakeSpsRuvResult(
            corrected_preprocessing_output=output,
            diagnostics={"status": "applied", "method": "sps_ruv_style"},
            provenance=provenance,
        )


@pytest.mark.parametrize(
    (
        "stage_key",
        "requirement",
        "stage_kwargs",
        "metadata_kwargs",
        "match",
    ),
    [
        (
            "custom_random_seed",
            QuantitativeEvidenceRequirement.RANDOM_SEED,
            {},
            {},
            "missing_requirement='random_seed'",
        ),
        (
            "custom_row_audit",
            QuantitativeEvidenceRequirement.ROW_AUDIT,
            {},
            {},
            "missing_requirement='row_audit'",
        ),
        (
            "custom_table_fingerprints",
            QuantitativeEvidenceRequirement.TABLE_FINGERPRINTS,
            {},
            {
                "produced_output_tables": (),
                "provenance_input_tables": ("dataset.phospho",),
                "provenance_output_table": "dataset.phospho",
            },
            "missing_requirement='table_fingerprints'",
        ),
        (
            "custom_intensity_event",
            QuantitativeEvidenceRequirement.TYPED_INTENSITY_TRANSFORMATION_EVENT,
            {},
            {},
            "missing_requirement='typed_intensity_transformation_event'",
        ),
        (
            "custom_total_mapping",
            QuantitativeEvidenceRequirement.TOTAL_PROTEIN_ROW_MAPPING,
            {
                "transition_evidence": QuantitativeTransitionEvidence(
                    total_protein_corrected_row_count=1,
                    total_protein_uncorrected_row_count=0,
                ),
                "result_evidence": QuantitativeOperationEvidence(
                    total_protein_row_mapping=TotalProteinRowMappingEvidence(
                        input_phosphosite_row_count=1,
                        corrected_phosphosite_row_ids=("row_a",),
                        uncorrected_phosphosite_row_ids=(),
                        corrected_phosphosite_to_total_row_ids=(("row_a", "total_a"),),
                        total_protein_row_count=1,
                    )
                ),
            },
            {},
            "missing_requirement='total_protein_row_mapping'",
        ),
        (
            "custom_missingness_mask",
            QuantitativeEvidenceRequirement.MISSINGNESS_MASK,
            {},
            {},
            "missing_requirement='missingness_mask'",
        ),
    ],
)
def test_adversarial_custom_stage_missing_declared_evidence_fails_before_trace_acceptance(
    stage_key: str,
    requirement: QuantitativeEvidenceRequirement,
    stage_kwargs: dict[str, Any],
    metadata_kwargs: dict[str, Any],
    match: str,
) -> None:
    assert requirement is not QuantitativeEvidenceRequirement.NONE
    stage = _IdentityStage(stage_key, **stage_kwargs)
    metadata = _custom_metadata(
        stage_key=stage_key,
        required_evidence=(requirement,),
        **metadata_kwargs,
    )

    with pytest.raises(DatasetBuildError, match=match):
        PreprocessingPipeline(
            stage_registry=(stage,),
            stage_contract_registry=(metadata,),
        ).run_with_trace(_custom_state(stage.stage_key))


def test_positive_custom_stage_supplies_evidence_and_reconstruction_accepts_trace() -> (
    None
):
    stage_key = "custom_positive_evidence"
    stage = _IdentityStage(
        stage_key,
        diagnostics={"diagnostics": {"random_seed": 17}},
        add_observation_mask=True,
        result_evidence=QuantitativeOperationEvidence(
            row_audit=RowAuditEvidence(record_count=0)
        ),
    )
    metadata = _custom_metadata(
        stage_key=stage_key,
        required_evidence=(
            QuantitativeEvidenceRequirement.TABLE_FINGERPRINTS,
            QuantitativeEvidenceRequirement.MISSINGNESS_MASK,
            QuantitativeEvidenceRequirement.ROW_AUDIT,
            QuantitativeEvidenceRequirement.RANDOM_SEED,
        ),
        produced_output_tables=(
            PreprocessingStateTableKey.DATASET_PHOSPHO,
            PreprocessingStateTableKey.DATASET_IMPUTATION_OBSERVATION_MASK,
        ),
        provenance_input_tables=("dataset.phospho",),
        provenance_output_table="dataset.phospho",
    )

    _, trace = PreprocessingPipeline(
        stage_registry=(stage,),
        stage_contract_registry=(metadata,),
    ).run_with_trace(_custom_state(stage_key))

    record = trace[0]
    assert record.random_seed == 17
    assert record.quantitative_evidence is not None
    assert record.quantitative_evidence.row_audit is not None
    assert record.quantitative_evidence.observation_mask is not None

    processing_state = DatasetProcessingStateBuilder().build(
        plan=PreprocessingPlan(stage_order=(stage_key,)),
        intensity_scale_state=_established_linear_state(),
        preprocessing_trace=trace,
    )
    assert processing_state.intensity_scale.kind is IntensityScaleKind.LINEAR


def test_processing_state_rejects_tampered_trace_missing_required_seed() -> None:
    stage_key = "custom_tampered_seed"
    stage = _IdentityStage(
        stage_key,
        diagnostics={"diagnostics": {"random_seed": 23}},
    )
    metadata = _custom_metadata(
        stage_key=stage_key,
        required_evidence=(QuantitativeEvidenceRequirement.RANDOM_SEED,),
    )
    _, trace = PreprocessingPipeline(
        stage_registry=(stage,),
        stage_contract_registry=(metadata,),
    ).run_with_trace(_custom_state(stage_key))
    tampered = replace(trace[0], random_seed=None)

    with pytest.raises(DatasetBuildError, match="missing_requirement='random_seed'"):
        DatasetProcessingStateBuilder().build(
            plan=PreprocessingPlan(stage_order=(stage_key,)),
            intensity_scale_state=_established_linear_state(),
            preprocessing_trace=(tampered,),
        )


def test_processing_state_rejects_tampered_total_protein_trace_missing_mapping() -> (
    None
):
    plan = PreprocessingPlan(
        total_protein_correction_policy=TotalProteinCorrectionPolicy.SUBTRACT_LOG_TOTAL,
        intensity_transform_policy=IntensityTransformPolicy.LOG2,
        stage_order=(DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,),
    )
    state = PreprocessingState(
        phospho=_log2_phospho(),
        site_metadata=_site_metadata(_log2_phospho().index),
        sample_metadata=None,
        total=_total(_log2_phospho().columns),
        plan=plan,
    )
    _, trace = PreprocessingPipeline().run_with_trace(
        state,
        initial_quantitative_scale_kind=IntensityScaleKind.LOG2,
        initial_quantitative_meaning=QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
    )
    tampered = replace(trace[0], quantitative_evidence=None)

    with pytest.raises(
        DatasetBuildError,
        match="missing_requirement='total_protein_row_mapping'",
    ):
        DatasetProcessingStateBuilder().build(
            plan=plan,
            intensity_scale_state=_established_log2_state(),
            preprocessing_trace=(tampered,),
        )


def test_built_in_registry_evidence_requirements_are_exhaustive() -> None:
    validator_requirements = (
        QuantitativeOperationEvidenceValidator.supported_requirements()
    )
    assert validator_requirements == frozenset(QuantitativeEvidenceRequirement)

    representative_plans = (
        PreprocessingPlan.default(),
        PreprocessingPlan(
            intensity_transform_policy=IntensityTransformPolicy.LOG2,
            stage_order=(DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,),
        ),
        PreprocessingPlan(
            missing_data_policy=MissingDataPolicy.IMPUTE_MINPROB,
            intensity_transform_policy=IntensityTransformPolicy.LOG2,
            stage_order=(
                DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
                DATASET_PREPROCESSING_STAGE_MISSING_DATA,
            ),
        ),
        PreprocessingPlan(
            total_protein_correction_policy=(
                TotalProteinCorrectionPolicy.SUBTRACT_LOG_TOTAL
            ),
            intensity_transform_policy=IntensityTransformPolicy.LOG2,
            stage_order=(DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,),
        ),
        _sps_ruv_plan(),
    )
    for metadata in list_registered_preprocessing_stages():
        for plan in representative_plans:
            contract = metadata.resolve_quantitative_contract(plan)
            assert contract.required_evidence
            assert (
                QuantitativeEvidenceRequirement.NONE not in contract.required_evidence
                or len(contract.required_evidence) == 1
            )
            assert contract.required_evidence <= validator_requirements


@pytest.mark.parametrize(
    "case_name",
    [
        "log2_intensity_transformation",
        "minprob_imputation",
        "row_median_imputation",
        "localisation_filtering",
        "group_coverage_filtering",
        "site_matrix_construction",
        "fixed_effect_batch_correction",
        "sps_ruv_style_batch_correction",
        "total_protein_correction",
        "comparison_derivation",
    ],
)
def test_built_in_stage_normal_execution_supplies_required_evidence(
    case_name: str,
) -> None:
    state, pipeline, initial_scale, initial_meaning = _built_in_case(case_name)
    _current, trace = pipeline.run_with_trace(
        state,
        initial_quantitative_scale_kind=initial_scale,
        initial_quantitative_meaning=initial_meaning,
    )

    assert trace, case_name
    for record in trace:
        assert record.quantitative_contract is not None
        requirements = record.quantitative_contract.required_evidence
        assert requirements
        if requirements != frozenset({QuantitativeEvidenceRequirement.NONE}):
            assert record.quantitative_evidence is not None or any(
                requirement
                in {
                    QuantitativeEvidenceRequirement.RANDOM_SEED,
                    QuantitativeEvidenceRequirement.TYPED_INTENSITY_TRANSFORMATION_EVENT,
                    QuantitativeEvidenceRequirement.TABLE_FINGERPRINTS,
                }
                for requirement in requirements
            )


def test_registration_rejects_empty_and_mixed_none_evidence() -> None:
    with pytest.raises(
        InvalidTransformationStateError,
        match="requires required_evidence",
    ):
        _contract(required_evidence=())

    with pytest.raises(
        InvalidTransformationStateError,
        match="NONE may not be combined",
    ):
        _contract(
            required_evidence=(
                QuantitativeEvidenceRequirement.NONE,
                QuantitativeEvidenceRequirement.RANDOM_SEED,
            )
        )


def test_validator_exhaustiveness_failure_is_deterministic(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import phospy.science.datasets.preprocessing.quantitative_evidence as module

    reduced = {
        key: value
        for key, value in module._VALIDATOR_BY_REQUIREMENT.items()
        if key is not QuantitativeEvidenceRequirement.RANDOM_SEED
    }
    monkeypatch.setattr(module, "_VALIDATOR_BY_REQUIREMENT", reduced)

    with pytest.raises(DatasetBuildError, match="missing validators.*random_seed"):
        QuantitativeOperationEvidenceValidator.require_all_enum_members_handled()


def test_stage_registration_with_unknown_evidence_requirement_fails() -> None:
    with pytest.raises(InvalidTransformationStateError, match="unsupported"):
        _contract(required_evidence=("future_requirement",))  # type: ignore[arg-type]

    metadata = _custom_metadata(
        stage_key="custom_supported_requirement",
        required_evidence=(QuantitativeEvidenceRequirement.RANDOM_SEED,),
    )
    resolved = resolve_registered_preprocessing_stages((metadata,))
    assert any(item.stage_key == "custom_supported_requirement" for item in resolved)


def _built_in_case(
    case_name: str,
) -> tuple[
    PreprocessingState,
    PreprocessingPipeline,
    IntensityScaleKind,
    QuantitativeMeaning,
]:
    if case_name == "log2_intensity_transformation":
        phospho = _linear_phospho()
        return (
            PreprocessingState(
                phospho=phospho,
                site_metadata=_site_metadata(phospho.index),
                sample_metadata=None,
                total=_total(phospho.columns),
                plan=PreprocessingPlan(
                    intensity_transform_policy=IntensityTransformPolicy.LOG2,
                    stage_order=(DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,),
                ),
            ),
            PreprocessingPipeline(),
            IntensityScaleKind.LINEAR,
            QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        )
    if case_name == "minprob_imputation":
        phospho = _minprob_linear_phospho_with_missing()
        return (
            PreprocessingState(
                phospho=phospho,
                site_metadata=_site_metadata(phospho.index),
                sample_metadata=None,
                total=None,
                plan=PreprocessingPlan(
                    intensity_transform_policy=IntensityTransformPolicy.LOG2,
                    missing_data_policy=MissingDataPolicy.IMPUTE_MINPROB,
                    missing_data_seed=101,
                    missing_data_q=0.01,
                    missing_data_width=0.3,
                    missing_data_max_missing_fraction_per_row=0.5,
                    stage_order=(
                        DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
                        DATASET_PREPROCESSING_STAGE_MISSING_DATA,
                    ),
                ),
            ),
            PreprocessingPipeline(),
            IntensityScaleKind.LINEAR,
            QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        )
    if case_name == "row_median_imputation":
        phospho = _linear_phospho_with_missing()
        return (
            PreprocessingState(
                phospho=phospho,
                site_metadata=_site_metadata(phospho.index),
                sample_metadata=None,
                total=None,
                plan=PreprocessingPlan(
                    missing_data_policy=MissingDataPolicy.IMPUTE_ROW_MEDIAN,
                    missing_data_min_observed_values=1,
                    missing_data_input_scale=ImputationInputScale.LINEAR,
                    stage_order=(DATASET_PREPROCESSING_STAGE_MISSING_DATA,),
                ),
            ),
            PreprocessingPipeline(),
            IntensityScaleKind.LINEAR,
            QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        )
    if case_name == "localisation_filtering":
        phospho = _linear_phospho()
        return (
            PreprocessingState(
                phospho=phospho,
                site_metadata=_site_metadata_with_missing_localisation(phospho.index),
                sample_metadata=None,
                total=None,
                plan=PreprocessingPlan(
                    localisation_mode=(
                        LocalisationEligibilityMode.ALLOW_MISSING_WITH_WAIVER
                    ),
                    localisation_waiver_reason="unit-test waiver",
                    stage_order=(DATASET_PREPROCESSING_STAGE_LOCALISATION,),
                ),
            ),
            PreprocessingPipeline(),
            IntensityScaleKind.LINEAR,
            QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        )
    if case_name == "group_coverage_filtering":
        phospho = _group_coverage_phospho()
        return (
            PreprocessingState(
                phospho=phospho,
                site_metadata=_site_metadata(phospho.index),
                sample_metadata=_group_sample_metadata(phospho.columns),
                total=None,
                plan=PreprocessingPlan(
                    group_coverage_filter_enabled=True,
                    group_coverage_filter_group_column="group",
                    group_coverage_filter_min_finite_observations_per_group=1,
                    group_coverage_filter_min_groups_passing_threshold=2,
                    stage_order=(DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER,),
                ),
            ),
            PreprocessingPipeline(),
            IntensityScaleKind.LINEAR,
            QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        )
    if case_name == "site_matrix_construction":
        phospho = _linear_phospho()
        return (
            PreprocessingState(
                phospho=phospho,
                site_metadata=_site_matrix_metadata(phospho.index),
                sample_metadata=None,
                total=None,
                plan=PreprocessingPlan(
                    site_matrix_policy=SiteMatrixPolicy.BUILD_FROM_METADATA,
                    stage_order=(DATASET_PREPROCESSING_STAGE_SITE_MATRIX,),
                ),
            ),
            PreprocessingPipeline(),
            IntensityScaleKind.LINEAR,
            QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        )
    if case_name == "fixed_effect_batch_correction":
        phospho = _batch_phospho()
        return (
            PreprocessingState(
                phospho=phospho,
                site_metadata=_site_metadata(phospho.index),
                sample_metadata=_batch_sample_metadata(phospho.columns),
                total=None,
                plan=PreprocessingPlan(
                    batch_correction_method=(
                        DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH
                    ),
                    stage_order=(DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,),
                ),
            ),
            PreprocessingPipeline(
                batch_correction_metadata_validator=cast(
                    Any,
                    _NoopDesignValidator(),
                ),
                batch_correction_adequacy_validator=cast(
                    Any,
                    _NoopAdequacyValidator(),
                ),
            ),
            IntensityScaleKind.LOG2,
            QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        )
    if case_name == "sps_ruv_style_batch_correction":
        phospho = _batch_phospho()
        return (
            PreprocessingState(
                phospho=phospho,
                site_metadata=_site_matrix_metadata(phospho.index),
                sample_metadata=_batch_sample_metadata(phospho.columns),
                total=None,
                plan=_sps_ruv_plan(),
            ),
            PreprocessingPipeline(
                batch_correction_runner=cast(Any, _FakeSpsRuvRunner()),
            ),
            IntensityScaleKind.LOG2,
            QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        )
    if case_name == "total_protein_correction":
        phospho = _log2_phospho()
        return (
            PreprocessingState(
                phospho=phospho,
                site_metadata=_site_metadata(phospho.index),
                sample_metadata=None,
                total=_total(phospho.columns),
                plan=PreprocessingPlan(
                    total_protein_correction_policy=(
                        TotalProteinCorrectionPolicy.SUBTRACT_LOG_TOTAL
                    ),
                    intensity_transform_policy=IntensityTransformPolicy.LOG2,
                    stage_order=(DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,),
                ),
            ),
            PreprocessingPipeline(),
            IntensityScaleKind.LOG2,
            QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        )
    if case_name == "comparison_derivation":
        phospho = _linear_phospho()
        return (
            PreprocessingState(
                phospho=phospho,
                site_metadata=_site_metadata(phospho.index),
                sample_metadata=_comparison_sample_metadata(phospho.columns),
                total=None,
                plan=PreprocessingPlan(
                    comparison_building_policy=(
                        ComparisonBuildingPolicy.SAMPLE_METADATA_PAIRS
                    ),
                    stage_order=(DATASET_PREPROCESSING_STAGE_COMPARISONS,),
                ),
            ),
            PreprocessingPipeline(),
            IntensityScaleKind.LINEAR,
            QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        )
    raise AssertionError(f"unsupported built-in evidence case {case_name!r}")


def _custom_metadata(
    *,
    stage_key: str,
    required_evidence: tuple[QuantitativeEvidenceRequirement, ...],
    consumed_input_tables: tuple[PreprocessingStateTableKey, ...] = (
        PreprocessingStateTableKey.DATASET_PHOSPHO,
    ),
    produced_output_tables: tuple[PreprocessingStateTableKey, ...] = (
        PreprocessingStateTableKey.DATASET_PHOSPHO,
    ),
    provenance_input_tables: tuple[str, ...] = (),
    provenance_output_table: str | None = None,
) -> PreprocessingStageMetadata:
    return PreprocessingStageMetadata(
        stage_key=stage_key,
        display_label=stage_key,
        provenance_stage=stage_key,
        operation_name=lambda _plan: "custom",
        serialize_parameters=lambda _plan: {},
        consumed_input_tables=consumed_input_tables,
        produced_output_tables=produced_output_tables,
        quantitative_contract=_contract(
            required_evidence=required_evidence,
            provenance_input_tables=provenance_input_tables,
            provenance_output_table=provenance_output_table,
        ),
        diagnostics_metadata={"known_diagnostics_fields": ("random_seed",)},
    )


def _contract(
    *,
    required_evidence: tuple[QuantitativeEvidenceRequirement | str, ...],
    provenance_input_tables: tuple[str, ...] = (),
    provenance_output_table: str | None = None,
) -> QuantitativeOperationContract:
    meanings = frozenset({QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE})
    scales = frozenset({IntensityScaleKind.LINEAR})
    return QuantitativeOperationContract(
        accepted_input_scale_kinds=scales,
        accepted_quantitative_meanings=meanings,
        output_scale_transition=preserve_scale_transition(scales),
        output_meaning_transition=preserve_meaning_transition(meanings),
        preserves_abundance=True,
        negative_domain_policy=NegativeDomainPolicy.PRESERVES_INPUT_DOMAIN,
        required_evidence=cast(
            frozenset[QuantitativeEvidenceRequirement],
            frozenset(required_evidence),
        ),
        reversibility=QuantitativeReversibilityKind.REVERSIBLE,
        information_loss=QuantitativeInformationLossKind.NONE,
        provenance_input_tables=provenance_input_tables,
        provenance_output_table=provenance_output_table,
    )


def _custom_state(stage_key: str) -> PreprocessingState:
    phospho = _linear_phospho()
    return PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan(stage_order=(stage_key,)),
    )


def _linear_phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [4.0, 8.0],
            "sample_b": [6.0, 12.0],
        },
        index=pd.Index(["row_a", "row_b"], name="source_row"),
    )


def _log2_phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [2.0, 3.0],
            "sample_b": [2.5, 3.5],
        },
        index=pd.Index(["MAPK14", "AKT1"], name="gene_symbol"),
    )


def _linear_phospho_with_missing() -> pd.DataFrame:
    frame = _linear_phospho()
    frame.loc["row_b", "sample_b"] = float("nan")
    return frame


def _log2_phospho_with_missing() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [2.0, float("nan"), float("nan")],
            "sample_b": [2.5, 3.0, float("nan")],
            "sample_c": [3.0, 3.5, 4.0],
        },
        index=pd.Index(["row_keep", "row_impute", "row_drop"]),
    )


def _minprob_linear_phospho_with_missing() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [4.0, float("nan"), float("nan")],
            "sample_b": [5.0, 8.0, float("nan")],
            "sample_c": [8.0, 9.0, 16.0],
        },
        index=pd.Index(["row_keep", "row_impute", "row_drop"]),
    )


def _group_coverage_phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 1.0],
            "sample_b": [2.0, float("nan")],
            "sample_c": [3.0, float("nan")],
            "sample_d": [4.0, 4.0],
        },
        index=pd.Index(["row_keep", "row_drop"]),
    )


def _batch_phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [2.0, 3.0],
            "sample_b": [2.2, 3.1],
            "sample_c": [4.0, 5.0],
            "sample_d": [4.2, 5.1],
        },
        index=pd.Index(
            [
                _site_key("MAPK14", "Y", "182"),
                _site_key("SRC", "Y", "416"),
            ],
            name="site_key",
        ),
    )


def _total(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            str(columns[0]): [1.0, 1.5],
            str(columns[1]): [1.2, 1.7],
        },
        index=pd.Index(["MAPK14", "AKT1"], name="gene_symbol"),
    )


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    labels = [str(item) for item in index.tolist()]
    gene_symbols = [f"GENE_{position}" for position, _label in enumerate(labels)]
    if len(labels) >= 1:
        gene_symbols[0] = "MAPK14"
    if len(labels) >= 2:
        gene_symbols[1] = "AKT1"
    sites = [f"S{position + 1}" for position, _label in enumerate(labels)]
    if len(labels) >= 1:
        sites[0] = "Y182"
    if len(labels) >= 2:
        sites[1] = "T308"
    return pd.DataFrame(
        {
            "gene_symbol": gene_symbols,
            "site": sites,
            "site_sequence": [f"SEQ_{position}" for position in range(len(labels))],
            "localisation_confidence": [
                0.95 - (0.01 * position) for position in range(len(labels))
            ],
        },
        index=index.copy(),
    )


def _site_metadata_with_missing_localisation(index: pd.Index) -> pd.DataFrame:
    metadata = _site_metadata(index)
    metadata.loc[metadata.index[0], "localisation_confidence"] = pd.NA
    return metadata


def _site_matrix_metadata(index: pd.Index) -> pd.DataFrame:
    display_ids = ["MAPK14;Y182;", "SRC;Y416;"][: len(index)]
    site_keys = site_key_index_from_display_ids(display_ids)
    return pd.DataFrame(
        {
            "site_key": site_keys.astype(str).tolist(),
            "display_id": display_ids,
            "gene_symbol": ["MAPK14", "SRC"][: len(index)],
            "site": ["Y182", "Y416"][: len(index)],
            "site_sequence": ["SEQ_A", "SEQ_S"][: len(index)],
            "localisation_confidence": [0.95, 0.9][: len(index)],
        },
        index=index.copy(),
    )


def _group_sample_metadata(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {"group": ["left", "left", "right", "right"][: len(columns)]},
        index=columns.copy(),
    )


def _batch_sample_metadata(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "batch": ["batch_a", "batch_b", "batch_a", "batch_b"][: len(columns)],
            "condition": ["control", "control", "treated", "treated"][: len(columns)],
            "replicate": ["r1", "r2", "r1", "r2"][: len(columns)],
        },
        index=columns.copy(),
    )


def _comparison_sample_metadata(columns: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {"comparison_group": ["control", "treated"][: len(columns)]},
        index=columns.copy(),
    )


def _sps_ruv_plan() -> PreprocessingPlan:
    return PreprocessingPlan.from_config(
        DatasetPreprocessingConfig(
            batch_correction=SpsRuvBatchCorrectionConfig(
                control_site_set=ControlSiteSet.from_site_keys(
                    (
                        _site_key("MAPK14", "Y", "182"),
                        _site_key("SRC", "Y", "416"),
                    ),
                    source_metadata=ControlSiteSourceMetadata(
                        organism="rat",
                        identifier_namespace="site_key",
                        source_name="manual-controls",
                        source_version="manual-v1",
                        license="caller local use",
                        redistribution="not redistributed",
                    ),
                ),
                batch_column="batch",
                condition_columns=("condition",),
                replicate_column="replicate",
                missingness_policy=CorrectionMissingnessPolicy(),
                n_unwanted_factors=1,
                diagnostics_enabled=True,
                provenance_enabled=True,
            )
        )
    )


def _site_key(gene_symbol: str, residue: str, position: str) -> str:
    return str(
        site_key_index_from_display_ids(
            [f"{gene_symbol};{residue}{position};"],
            protein_namespace="gene_symbol",
        )[0]
    )


def _established_linear_state() -> IntensityScaleState:
    return establish_intensity_scale_state(
        IntensityScaleState(
            phospho=MatrixIntensityScaleState.linear(
                established_by="tests.quantitative_evidence"
            ),
            total=None,
            quantity=QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        ),
        established_via="tests.quantitative_evidence",
        establishment_mode=IntensityScaleEstablishmentMode.DECLARED,
        evidence_level=IntensityScaleEvidenceLevel.DECLARED_BY_USER,
        _authority=dataset_resolver_establishment_authority(),
    )


def _established_log2_state() -> IntensityScaleState:
    return establish_intensity_scale_state(
        IntensityScaleState(
            phospho=MatrixIntensityScaleState.log2(
                established_by="tests.quantitative_evidence"
            ),
            total=MatrixIntensityScaleState.log2(
                established_by="tests.quantitative_evidence"
            ),
            quantity=QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        ),
        established_via="tests.quantitative_evidence",
        establishment_mode=IntensityScaleEstablishmentMode.DECLARED,
        evidence_level=IntensityScaleEvidenceLevel.DECLARED_BY_USER,
        _authority=dataset_resolver_establishment_authority(),
    )
