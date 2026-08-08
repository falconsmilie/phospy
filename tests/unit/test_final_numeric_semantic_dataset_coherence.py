from __future__ import annotations

from collections.abc import Mapping
from dataclasses import replace

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder, AnalysisReadyPhosphoDataset
from phospy.advanced import (
    DatasetIntensityTransformConfig,
    DatasetNormalisationConfig,
)
from phospy.api import (
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    Organism,
)
from phospy.errors.transformations import InvalidTransformationStateError
from phospy.errors.validation import (
    DatasetValidationError,
    TransformationValidationError,
)
from phospy.provenance import from_payload, to_payload
from phospy.provenance.models import TrustedDatasetConstructionEvidence
from phospy.science.datasets.builders.preprocessing import (
    build_dataset_processing_state,
)
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.science.transformations.models import (
    IntensityScaleKind,
    IntensityScaleState,
    MatrixIntensityScaleState,
    QuantitativeMeaning,
)
from phospy.science.transformations.state_coherence import (
    ObservedNumericDomain,
    require_quantitative_numeric_domain_coherence,
)
from tests.support.analysis_ready_dataset_factories import (
    complete_trusted_dataset_construction_assertions_for_tests,
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_log2_intensity_scale_state,
    with_restored_quantitative_meaning_for_tests,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

_CENTRED_Y_SEQUENCE = "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA"
_CENTRED_T_SEQUENCE = "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA"


def _site_index() -> pd.Index:
    return protein_site_key_index(
        protein_identifiers=["MAPK14", "AKT1"],
        sites=["Y182", "T308"],
    )


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            "display_id": ["MAPK14;Y182;", "AKT1;T308;"],
            **site_key_context_columns(index),
            "gene_symbol": ["MAPK14", "AKT1"],
            "protein_id": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [_CENTRED_Y_SEQUENCE, _CENTRED_T_SEQUENCE],
            "localisation_confidence": [0.95, 0.96],
        },
        index=index.copy(),
    )


def _phospho(values: list[list[float]]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [row[0] for row in values],
            "sample_b": [row[1] for row in values],
        },
        index=_site_index(),
    )


def _processing_state_for(
    intensity_scale_state: IntensityScaleState,
):
    return build_dataset_processing_state(
        plan=PreprocessingPlan.default(),
        intensity_scale_state=intensity_scale_state,
    )


def _state(
    *,
    scale_kind: IntensityScaleKind,
    meaning: QuantitativeMeaning,
) -> IntensityScaleState:
    matrix_state = (
        MatrixIntensityScaleState.linear(established_by="test")
        if scale_kind is IntensityScaleKind.LINEAR
        else MatrixIntensityScaleState.log2(established_by="test")
    )
    return IntensityScaleState(phospho=matrix_state, quantity=meaning)


_DOMAIN_VALUES = {
    ObservedNumericDomain.ZERO_ONLY: [[0.0, 0.0], [0.0, 0.0]],
    ObservedNumericDomain.NON_NEGATIVE: [[0.0, 2.0], [3.0, 4.0]],
    ObservedNumericDomain.SIGNED: [[-1.0, 2.0], [3.0, -0.25]],
    ObservedNumericDomain.NEGATIVE_ONLY: [[-3.0, -2.0], [-1.5, -0.25]],
}

_COMPATIBLE_MEANINGS_BY_SCALE = {
    IntensityScaleKind.LINEAR: {
        QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        QuantitativeMeaning.ACTIVITY_SCORE,
        QuantitativeMeaning.UNKNOWN,
    },
    IntensityScaleKind.LOG2: {
        QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
        QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE,
        QuantitativeMeaning.DIFFERENTIAL_EFFECT_SIZE,
        QuantitativeMeaning.ACTIVITY_SCORE,
        QuantitativeMeaning.MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE,
        QuantitativeMeaning.UNKNOWN,
    },
}


@pytest.mark.parametrize("scale_kind", list(IntensityScaleKind))
@pytest.mark.parametrize("meaning", list(QuantitativeMeaning))
@pytest.mark.parametrize("domain", list(ObservedNumericDomain))
def test_scale_meaning_numeric_domain_matrix(
    scale_kind: IntensityScaleKind,
    meaning: QuantitativeMeaning,
    domain: ObservedNumericDomain,
) -> None:
    compatible_meanings = _COMPATIBLE_MEANINGS_BY_SCALE[scale_kind]
    if meaning not in compatible_meanings:
        with pytest.raises(InvalidTransformationStateError):
            _state(scale_kind=scale_kind, meaning=meaning)
        return

    state = _state(scale_kind=scale_kind, meaning=meaning)
    matrix = _phospho(_DOMAIN_VALUES[domain])
    should_accept = meaning is not QuantitativeMeaning.UNKNOWN and (
        meaning is not QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE
        or domain
        in {ObservedNumericDomain.ZERO_ONLY, ObservedNumericDomain.NON_NEGATIVE}
    )
    if should_accept:
        require_quantitative_numeric_domain_coherence(
            phospho=matrix,
            total=None,
            intensity_scale_state=state,
        )
    else:
        with pytest.raises(TransformationValidationError):
            require_quantitative_numeric_domain_coherence(
                phospho=matrix,
                total=None,
                intensity_scale_state=state,
            )


def test_linear_phosphosite_abundance_with_signed_values_is_rejected() -> None:
    state = supported_linear_intensity_scale_state(has_total_matrix=False)
    phospho = _phospho([[-1.0, 2.0], [3.0, 4.0]])

    with pytest.raises(DatasetValidationError) as exc_info:
        trusted_analysis_ready_dataset_from_tables(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            organism=Organism.RAT,
            intensity_scale_state=state,
            processing_state=_processing_state_for(state),
        )

    message = str(exc_info.value)
    assert "scale='linear'" in message
    assert "meaning='phosphosite_abundance'" in message
    assert "observed_numeric_domain='signed'" in message
    assert "linear phosphosite_abundance must be non-negative" in message


def test_public_builder_rejects_previously_warning_only_false_linear_state() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [-1.0, 3.0], "sample_b": [2.0, 4.0]},
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "protein_id": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [_CENTRED_Y_SEQUENCE, _CENTRED_T_SEQUENCE],
            "localisation_confidence": [0.95, 0.96],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(
        DatasetValidationError,
        match="linear phosphosite_abundance must be non-negative",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                organism=Organism.RAT,
                input_intensity_scale="linear",
                preprocessing_config=DatasetPreprocessingConfig(
                    intensity_transform=DatasetIntensityTransformConfig(
                        policy="identity"
                    )
                ),
            )
        )


def test_signed_centered_log_quantity_remains_valid_when_labelled_log_abundance() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.0, 3.0], "sample_b": [4.0, 0.0]},
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "protein_id": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [_CENTRED_Y_SEQUENCE, _CENTRED_T_SEQUENCE],
            "localisation_confidence": [0.95, 0.96],
        },
        index=phospho.index.copy(),
    )

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="log2",
            preprocessing_config=DatasetPreprocessingConfig(
                normalisation=DatasetNormalisationConfig(policy="median_center")
            ),
        )
    )

    assert dataset.intensity_scale_state.quantity is (
        QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE
    )
    assert float(dataset.phospho.min().min()) < 0.0


def test_signed_fold_change_quantity_serializes_and_fingerprints() -> None:
    base_state = supported_log2_intensity_scale_state(has_total_matrix=False)
    signed_state = with_restored_quantitative_meaning_for_tests(
        base_state,
        QuantitativeMeaning.CONTRAST_LOG2_FOLD_CHANGE,
    )
    signed_phospho = _phospho([[-1.25, 2.0], [0.0, -0.5]])
    positive_phospho = _phospho([[1.25, 2.0], [0.0, 0.5]])

    signed_dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=signed_phospho,
        site_metadata=_site_metadata(signed_phospho.index),
        organism=Organism.RAT,
        intensity_scale_state=signed_state,
        processing_state=_processing_state_for(signed_state),
    )
    positive_dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=positive_phospho,
        site_metadata=_site_metadata(positive_phospho.index),
        organism=Organism.RAT,
        intensity_scale_state=signed_state,
        processing_state=_processing_state_for(signed_state),
    )

    assert signed_dataset.provenance is not None
    payload = to_payload(signed_dataset.provenance)
    assert to_payload(from_payload(payload)) == payload
    signed_fingerprint = _phospho_exact_fingerprint(signed_dataset)
    positive_fingerprint = _phospho_exact_fingerprint(positive_dataset)
    assert signed_fingerprint != positive_fingerprint


def test_unknown_quantitative_meaning_is_not_established_by_numeric_domain() -> None:
    base_state = supported_linear_intensity_scale_state(has_total_matrix=False)
    unknown_state = with_restored_quantitative_meaning_for_tests(
        base_state,
        QuantitativeMeaning.UNKNOWN,
    )
    phospho = _phospho([[1.0, 2.0], [3.0, 4.0]])

    with pytest.raises(
        DatasetValidationError,
        match="unknown quantitative meaning has no numeric-domain contract",
    ):
        trusted_analysis_ready_dataset_from_tables(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            organism=Organism.RAT,
            intensity_scale_state=unknown_state,
            processing_state=_processing_state_for(unknown_state),
        )


def test_quantitative_meaning_waiver_alone_cannot_bypass_numeric_domain_rule() -> None:
    state = supported_linear_intensity_scale_state(has_total_matrix=False)
    phospho = _phospho([[-1.0, 2.0], [3.0, 4.0]])
    assertions = replace(
        complete_trusted_dataset_construction_assertions_for_tests(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            intensity_scale_state=state,
            processing_state=_processing_state_for(state),
        ),
        quantitative_meaning=TrustedDatasetConstructionEvidence.waiver(
            reason="legacy export lacks independent quantitative meaning evidence",
            policy="legacy_quantitative_meaning_waiver",
        ),
        numeric_semantic_domain=None,
    )

    with pytest.raises(
        DatasetValidationError,
        match="linear phosphosite_abundance must be non-negative",
    ):
        trusted_analysis_ready_dataset_from_tables(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            organism=Organism.RAT,
            intensity_scale_state=state,
            processing_state=_processing_state_for(state),
            trusted_construction_assertions=assertions,
        )


def test_typed_numeric_semantic_waiver_can_bypass_and_remains_visible() -> None:
    state = supported_linear_intensity_scale_state(has_total_matrix=False)
    phospho = _phospho([[-1.0, 2.0], [3.0, 4.0]])
    assertions = replace(
        complete_trusted_dataset_construction_assertions_for_tests(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            intensity_scale_state=state,
            processing_state=_processing_state_for(state),
        ),
        numeric_semantic_domain=TrustedDatasetConstructionEvidence.waiver(
            reason="legacy centred export intentionally replayed for audit",
            policy="explicit_numeric_semantic_domain_waiver",
        ),
    )

    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        organism=Organism.RAT,
        intensity_scale_state=state,
        processing_state=_processing_state_for(state),
        trusted_construction_assertions=assertions,
    )

    assert dataset.trusted_construction_assertions is not None
    assert "numeric_semantic_domain" in (
        dataset.trusted_construction_assertions.waived_assertions
    )
    assert dataset.provenance is not None
    construction = dataset.provenance.workflow_parameters["construction"]
    assert isinstance(construction, Mapping)
    trusted_assertions = construction["trusted_construction_assertions"]
    assert isinstance(trusted_assertions, Mapping)
    waiver = trusted_assertions["numeric_semantic_domain"]
    assert isinstance(waiver, Mapping)
    assert waiver["kind"] == "waiver"
    assert "numeric_semantic_domain" in trusted_assertions["waived_assertions"]
    assert construction["trusted_construction_assertion_fingerprint"] == (
        assertions.assertion_fingerprint
    )


def test_linear_total_protein_abundance_rejects_negative_total_values() -> None:
    state = supported_linear_intensity_scale_state(has_total_matrix=True)
    phospho = _phospho([[1.0, 2.0], [3.0, 4.0]])
    total = pd.DataFrame(
        {"sample_a": [-10.0, 20.0], "sample_b": [11.0, 21.0]},
        index=pd.Index(["MAPK14", "AKT1"], name="protein_id"),
    )

    with pytest.raises(DatasetValidationError) as exc_info:
        trusted_analysis_ready_dataset_from_tables(
            phospho=phospho,
            site_metadata=_site_metadata(phospho.index),
            total=total,
            organism=Organism.RAT,
            intensity_scale_state=state,
            processing_state=_processing_state_for(state),
        )

    message = str(exc_info.value)
    assert "table='dataset.total'" in message
    assert "meaning='total_protein_abundance'" in message
    assert "linear total_protein_abundance must be non-negative" in message


def _phospho_exact_fingerprint(dataset: AnalysisReadyPhosphoDataset) -> str:
    assert dataset.provenance is not None
    for fingerprint in dataset.provenance.output_tables:
        if fingerprint.name == "dataset.phospho":
            return fingerprint.exact_hash_value
    raise AssertionError("dataset.phospho fingerprint not found")
