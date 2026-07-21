from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceContextCompatibilityPolicy,
    SampleDesignRecord,
    SignalomeWorkflowRequest,
)
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
)
from phospy.errors import WorkflowValidationError
from phospy.errors.transformations import InvalidTransformationStateError
from phospy.science.references.models import ReferenceBundle
from phospy.science.transformations._authority import (
    bundle_quantitative_meaning_restoration_authority,
)
from phospy.science.transformations.models import (
    IntensityScaleState,
    MatrixIntensityScaleState,
    QuantitativeMeaning,
    QuantitativeMeaningEvidenceMode,
    QuantitativeMeaningTransitionProvenance,
)
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator
from tests.support.intensity_scale_states import (
    supported_log2_intensity_scale_state,
    supported_log2_processing_state,
)
from tests.support.signalome_config import build_signalome_config
from tests.support.site_keys import protein_site_key_index, site_key_context_columns
from tests.support.unsafe_dataset_states import (
    unsafe_replace_dataset_intensity_scale_state,
)


def _dataset() -> AnalysisReadyPhosphoDataset:
    genes = ["MAPK14", "AKT1"]
    sites = ["Y182", "T308"]
    site_index = protein_site_key_index(protein_identifiers=genes, sites=sites)
    display_ids = ["MAPK14;Y182;", "AKT1;T308;"]
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {
                "A_1": [1.0, 2.0],
                "A_2": [1.1, 2.1],
                "B_1": [2.0, 2.2],
                "B_2": [2.1, 2.3],
            },
            index=site_index.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_index.astype(str).tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_index),
                "gene_symbol": genes,
                "site": sites,
                "site_sequence": [("A" * 15) + site[0] + ("A" * 15) for site in sites],
                "protein_id": genes,
            },
            index=site_index.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def _dataset_with_meaning(
    meaning: QuantitativeMeaning,
) -> AnalysisReadyPhosphoDataset:
    dataset = _dataset()
    provenance = QuantitativeMeaningTransitionProvenance(
        source_quantity=dataset.intensity_scale_state.quantity,
        target_quantity=meaning,
        operation_id=("tests.workflow_quantitative_input_contracts.restore_test_state"),
        producer_id="tests.unit.test_workflow_quantitative_input_contracts",
        evidence_mode=(
            QuantitativeMeaningEvidenceMode.RESTORED_FROM_TRUSTED_SERIALIZED_PROVENANCE
        ),
    )
    unsafe_replace_dataset_intensity_scale_state(
        dataset,
        dataset.intensity_scale_state.restore_quantitative_meaning_provenance(
            provenance=provenance,
            authority=bundle_quantitative_meaning_restoration_authority(),
        ),
    )
    return dataset


def _reference_bundle() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "AKT1"],
                "substrate_site": ["MAPK14;Y182;", "AKT1;T308;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                    "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
                ]
            },
            index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
        ),
    )


def _kinase_request(dataset: AnalysisReadyPhosphoDataset) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=dataset,
        references=_reference_bundle(),
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        activity_config=None,
    )


def _signalome_request(
    dataset: AnalysisReadyPhosphoDataset,
) -> SignalomeWorkflowRequest:
    site_index = dataset._borrow_phospho_frame().index.copy()
    score_matrix = pd.DataFrame(
        {"MAP2K6": [0.5, 0.8], "AKT1": [0.7, 0.4]},
        index=site_index.copy(),
    )
    prediction_matrix = pd.DataFrame(
        {"MAP2K6": [0.9, 0.2], "AKT1": [0.2, 0.9]},
        index=site_index.copy(),
    )
    return SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=dataset,
            references=_reference_bundle(),
            scoring_result=KinaseScoringResult(profile_scores=score_matrix),
            prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
        ),
        config=build_signalome_config(
            substrate_support_cutoff=0.5,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
    )


def _differential_request(
    dataset: AnalysisReadyPhosphoDataset,
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=dataset,
        design=ExperimentalDesign(
            samples=(
                SampleDesignRecord(sample_id="A_1", condition="A"),
                SampleDesignRecord(sample_id="A_2", condition="A"),
                SampleDesignRecord(sample_id="B_1", condition="B"),
                SampleDesignRecord(sample_id="B_2", condition="B"),
            )
        ),
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
        ),
    )


def _assert_expected_actual_meaning_error(
    exc: pytest.ExceptionInfo[WorkflowValidationError],
    *,
    context: str,
    expected_meaning: str,
    actual_meaning: QuantitativeMeaning,
) -> None:
    message = str(exc.value)
    assert f"{context} requires quantitative meaning in" in message
    assert repr(expected_meaning) in message
    assert f"got {actual_meaning.value!r}" in message


def test_differential_validator_accepts_log_abundance_input_meaning() -> None:
    dataset = _dataset()

    validated = DifferentialAnalysisValidator().run(_differential_request(dataset))

    assert validated.dataset is dataset


@pytest.mark.parametrize(
    "meaning",
    [
        pytest.param(
            QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE,
            id="contrast-logfc",
        ),
        pytest.param(QuantitativeMeaning.ACTIVITY_SCORE, id="activity-score"),
    ],
)
def test_differential_validator_rejects_incompatible_quantitative_meanings(
    meaning: QuantitativeMeaning,
) -> None:
    dataset = _dataset_with_meaning(meaning)

    with pytest.raises(WorkflowValidationError) as exc_info:
        DifferentialAnalysisValidator().run(_differential_request(dataset))

    _assert_expected_actual_meaning_error(
        exc_info,
        context="differential workflow request dataset",
        expected_meaning=QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE.value,
        actual_meaning=meaning,
    )


def test_kinase_validator_accepts_log_abundance_input_meaning() -> None:
    dataset = _dataset()
    request = _kinase_request(dataset)

    validated = KinaseWorkflowValidator().run(request)

    assert validated is request


@pytest.mark.parametrize(
    "meaning",
    [
        pytest.param(
            QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE,
            id="contrast-logfc",
        ),
        pytest.param(
            QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE,
            id="differential-effect",
        ),
    ],
)
def test_kinase_validator_rejects_incompatible_quantitative_meanings(
    meaning: QuantitativeMeaning,
) -> None:
    dataset = _dataset_with_meaning(meaning)

    with pytest.raises(WorkflowValidationError) as exc_info:
        KinaseWorkflowValidator().run(_kinase_request(dataset))

    _assert_expected_actual_meaning_error(
        exc_info,
        context="kinase workflow request dataset",
        expected_meaning=QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE.value,
        actual_meaning=meaning,
    )


def test_signalome_validator_accepts_log_abundance_input_meaning() -> None:
    dataset = _dataset()
    request = _signalome_request(dataset)

    validated = SignalomeWorkflowValidator().run(request)

    assert validated is request


@pytest.mark.parametrize(
    "meaning",
    [
        pytest.param(
            QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE,
            id="contrast-logfc",
        ),
        pytest.param(QuantitativeMeaning.ACTIVITY_SCORE, id="activity-score"),
    ],
)
def test_signalome_validator_rejects_incompatible_quantitative_meanings(
    meaning: QuantitativeMeaning,
) -> None:
    dataset = _dataset_with_meaning(meaning)

    with pytest.raises(WorkflowValidationError) as exc_info:
        SignalomeWorkflowValidator().run(_signalome_request(dataset))

    _assert_expected_actual_meaning_error(
        exc_info,
        context="signalome workflow request kinase_result.dataset",
        expected_meaning=QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE.value,
        actual_meaning=meaning,
    )


def test_log_abundance_meaning_requires_log2_scale_at_state_construction() -> None:
    with pytest.raises(
        InvalidTransformationStateError,
        match="phosphosite_log_abundance.*log2 intensity scale",
    ):
        IntensityScaleState(
            phospho=MatrixIntensityScaleState.linear(established_by="test.linear"),
            total=None,
            quantity=QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        )


def test_raw_abundance_meaning_requires_linear_scale_at_state_construction() -> None:
    with pytest.raises(
        InvalidTransformationStateError,
        match="phosphosite_abundance.*linear intensity scale",
    ):
        IntensityScaleState(
            phospho=MatrixIntensityScaleState.log2(established_by="test.log2"),
            total=None,
            quantity=QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        )


def test_contrast_logfc_meaning_requires_log2_scale_at_state_construction() -> None:
    with pytest.raises(
        InvalidTransformationStateError,
        match="contrast_log2_fold_change.*log2 intensity scale",
    ):
        IntensityScaleState(
            phospho=MatrixIntensityScaleState.linear(established_by="test.linear"),
            total=None,
            quantity=QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE,
        )


def test_contrast_matrix_cannot_enter_kinase_abundance_workflow() -> None:
    dataset = _dataset_with_meaning(QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE)

    with pytest.raises(WorkflowValidationError, match="requires quantitative meaning"):
        KinaseWorkflowValidator().run(_kinase_request(dataset))


def test_activity_score_matrix_cannot_enter_kinase_abundance_workflow() -> None:
    dataset = _dataset_with_meaning(QuantitativeMeaning.ACTIVITY_SCORE)

    with pytest.raises(WorkflowValidationError, match="activity_score"):
        KinaseWorkflowValidator().run(_kinase_request(dataset))


def test_differential_effect_matrix_cannot_enter_differential_input() -> None:
    dataset = _dataset_with_meaning(QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE)

    with pytest.raises(WorkflowValidationError, match="differential_effect_size"):
        DifferentialAnalysisValidator().run(_differential_request(dataset))


@pytest.mark.parametrize(
    "run_validator",
    [
        pytest.param(
            lambda dataset: DifferentialAnalysisValidator().run(
                _differential_request(dataset)
            ),
            id="differential",
        ),
        pytest.param(
            lambda dataset: KinaseWorkflowValidator().run(_kinase_request(dataset)),
            id="kinase",
        ),
        pytest.param(
            lambda dataset: SignalomeWorkflowValidator().run(
                _signalome_request(dataset)
            ),
            id="signalome",
        ),
    ],
)
def test_mixed_quantitative_meaning_is_rejected_by_default(
    run_validator,
) -> None:
    dataset = _dataset_with_meaning(
        QuantitativeMeaning.MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE
    )

    with pytest.raises(WorkflowValidationError, match="mixed"):
        run_validator(dataset)


@pytest.mark.parametrize(
    "run_validator",
    [
        pytest.param(
            lambda dataset: DifferentialAnalysisValidator().run(
                _differential_request(dataset)
            ),
            id="differential",
        ),
        pytest.param(
            lambda dataset: KinaseWorkflowValidator().run(_kinase_request(dataset)),
            id="kinase",
        ),
        pytest.param(
            lambda dataset: SignalomeWorkflowValidator().run(
                _signalome_request(dataset)
            ),
            id="signalome",
        ),
    ],
)
def test_unknown_quantitative_meaning_is_rejected_by_default(
    run_validator,
) -> None:
    dataset = _dataset_with_meaning(QuantitativeMeaning.UNKNOWN)

    with pytest.raises(WorkflowValidationError, match="unknown quantitative meaning"):
        run_validator(dataset)
