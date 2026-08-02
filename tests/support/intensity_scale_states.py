from __future__ import annotations

import pandas as pd

from phospy.science.datasets.builders.preprocessing import (
    build_dataset_processing_state,
)
from phospy.science.datasets.builders.transformation_resolver import (
    DatasetIntensityScaleResolver,
)
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.science.transformations._authority import (
    bundle_quantitative_meaning_restoration_authority,
)
from phospy.science.transformations.models import (
    IntensityScaleKind,
    IntensityScaleState,
    MatrixIntensityScaleState,
    QuantitativeMeaning,
    QuantitativeMeaningEvidenceMode,
    QuantitativeMeaningTransitionProvenance,
)
from phospy.science.transformations.transformers import IdentityTransformer


def supported_linear_intensity_scale_state(
    *,
    has_total_matrix: bool,
) -> IntensityScaleState:
    phospho = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index(["GENEA;S1;"], name="site_id"),
    )
    total = None
    if has_total_matrix:
        total = pd.DataFrame(
            {"sample_a": [1.0]},
            index=pd.Index(["GENEA"], name="protein_id"),
        )
    declared_state = IntensityScaleState(
        phospho=MatrixIntensityScaleState.linear(established_by="trusted.input"),
        total=(
            MatrixIntensityScaleState.linear(established_by="trusted.input")
            if has_total_matrix
            else None
        ),
    )
    established_state = (
        DatasetIntensityScaleResolver(transformer=IdentityTransformer())
        .run(
            phospho=phospho,
            total=total,
            expected_scale_kind=IntensityScaleKind.LINEAR,
            declared_input_scale_state=declared_state,
            input_declaration_source="tests.support.intensity_scale_states",
        )
        .intensity_scale_state
    )
    return build_dataset_processing_state(
        plan=PreprocessingPlan.default(),
        intensity_scale_state=established_state,
    ).intensity_scale


def supported_linear_processing_state(*, has_total_matrix: bool):
    intensity_scale_state = supported_linear_intensity_scale_state(
        has_total_matrix=has_total_matrix
    )
    return build_dataset_processing_state(
        plan=PreprocessingPlan.default(),
        intensity_scale_state=intensity_scale_state,
    )


def with_restored_quantitative_meaning_for_tests(
    intensity_scale_state: IntensityScaleState,
    meaning: QuantitativeMeaning,
) -> IntensityScaleState:
    provenance = QuantitativeMeaningTransitionProvenance(
        source_quantity=intensity_scale_state.quantity,
        target_quantity=meaning,
        operation_id="tests.support.intensity_scale_states.restore_meaning",
        producer_id="tests.support.intensity_scale_states",
        evidence_mode=(
            QuantitativeMeaningEvidenceMode.RESTORED_FROM_TRUSTED_SERIALIZED_PROVENANCE
        ),
    )
    return intensity_scale_state.restore_quantitative_meaning_provenance(
        provenance=provenance,
        authority=bundle_quantitative_meaning_restoration_authority(),
    )


def supported_log2_intensity_scale_state(
    *,
    has_total_matrix: bool,
) -> IntensityScaleState:
    phospho = pd.DataFrame(
        {"sample_a": [1.0]},
        index=pd.Index(["GENEA;S1;"], name="site_id"),
    )
    total = None
    if has_total_matrix:
        total = pd.DataFrame(
            {"sample_a": [1.0]},
            index=pd.Index(["GENEA"], name="protein_id"),
        )
    declared_state = IntensityScaleState(
        phospho=MatrixIntensityScaleState.log2(established_by="trusted.input"),
        total=(
            MatrixIntensityScaleState.log2(established_by="trusted.input")
            if has_total_matrix
            else None
        ),
    )
    established_state = (
        DatasetIntensityScaleResolver(transformer=IdentityTransformer())
        .run(
            phospho=phospho,
            total=total,
            expected_scale_kind=IntensityScaleKind.LOG2,
            declared_input_scale_state=declared_state,
            input_declaration_source="tests.support.intensity_scale_states",
        )
        .intensity_scale_state
    )
    return build_dataset_processing_state(
        plan=PreprocessingPlan.default(),
        intensity_scale_state=established_state,
    ).intensity_scale


def supported_log2_processing_state(*, has_total_matrix: bool):
    intensity_scale_state = supported_log2_intensity_scale_state(
        has_total_matrix=has_total_matrix
    )
    return build_dataset_processing_state(
        plan=PreprocessingPlan.default(),
        intensity_scale_state=intensity_scale_state,
    )


def supported_log2_intensity_scale_state_with_meaning(
    *,
    has_total_matrix: bool,
    meaning: QuantitativeMeaning,
) -> IntensityScaleState:
    return with_restored_quantitative_meaning_for_tests(
        supported_log2_intensity_scale_state(has_total_matrix=has_total_matrix),
        meaning,
    )


def supported_log2_processing_state_with_meaning(
    *,
    has_total_matrix: bool,
    meaning: QuantitativeMeaning,
):
    intensity_scale_state = supported_log2_intensity_scale_state_with_meaning(
        has_total_matrix=has_total_matrix,
        meaning=meaning,
    )
    return build_dataset_processing_state(
        plan=PreprocessingPlan.default(),
        intensity_scale_state=intensity_scale_state,
    )
