from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.datasets.builders.preprocessing import build_dataset_processing_state
from phospy.datasets.builders.transformation_resolver import (
    DatasetIntensityScaleResolver,
)
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.datasets.preprocessing.models import PreprocessingPlan
from phospy.errors.input import PhosPyInputError
from phospy.errors.transformations import (
    InvalidTransformationStateError,
    TransformationStateEstablishmentError,
    TransformerExecutionError,
)
from phospy.errors.validation import TransformationValidationError
from phospy.io.bundles._shared.intensity_scale_state import (
    intensity_scale_state_from_payload,
)
from phospy.io.bundles._shared.processing_state import (
    processing_state_from_payload,
)
from phospy.references.models import Organism
from phospy.transformations.contracts import TransformationResult
from phospy.transformations.models import (
    IntensityScaleKind,
    IntensityScaleState,
    MatrixIntensityScaleState,
    establish_intensity_scale_state,
)
from phospy.transformations.transformers import IdentityTransformer
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
)


def _phospho() -> pd.DataFrame:
    return pd.DataFrame({"sample_a": [1.0]}, index=["MAPK14;Y182;"])


def _total() -> pd.DataFrame:
    return pd.DataFrame({"sample_a": [2.0]}, index=["MAPK14"])


def _site_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["LDFGLARHTDDEMTGYVATRWYRAPEIMLNW"],
        },
        index=["MAPK14;Y182;"],
    )


def _processing_state_for(intensity_scale_state: IntensityScaleState):
    return build_dataset_processing_state(
        plan=PreprocessingPlan.default(),
        intensity_scale_state=intensity_scale_state,
    )


def test_resolver_fails_when_state_is_unknown_and_no_transformer_is_configured() -> (
    None
):
    resolver = DatasetIntensityScaleResolver()

    with pytest.raises(
        TransformationStateEstablishmentError, match="unable to establish"
    ):
        resolver.run(
            phospho=_phospho(),
            total=None,
        )


def test_resolver_uses_configured_transformer_to_establish_state() -> None:
    class DeclaredLog2Transformer:
        def run(
            self,
            phospho: pd.DataFrame,
            total: pd.DataFrame | None = None,
        ) -> TransformationResult:
            return TransformationResult(
                phospho=phospho,
                total=total,
                state=IntensityScaleState(
                    phospho=MatrixIntensityScaleState.log2(
                        established_by="test.transformer"
                    ),
                    total=(
                        MatrixIntensityScaleState.log2(
                            established_by="test.transformer"
                        )
                        if total is not None
                        else None
                    ),
                ),
            )

    resolver = DatasetIntensityScaleResolver(transformer=DeclaredLog2Transformer())

    resolved = resolver.run(
        phospho=_phospho(),
        total=_total(),
    )

    assert resolved.intensity_scale_state.phospho.kind.value == "log2"
    assert resolved.intensity_scale_state.total is not None
    assert resolved.intensity_scale_state.total.kind.value == "log2"
    assert resolved.intensity_scale_state.is_established
    assert resolved.intensity_scale_state.established_via is not None
    assert "DeclaredLog2Transformer" in resolved.intensity_scale_state.established_via


def test_resolver_rejects_invalid_transformer_state_contract() -> None:
    class MismatchedKindTransformer:
        def run(
            self,
            phospho: pd.DataFrame,
            total: pd.DataFrame | None = None,
        ) -> TransformationResult:
            return TransformationResult(
                phospho=phospho,
                total=total,
                state=IntensityScaleState(
                    phospho=MatrixIntensityScaleState.linear(
                        established_by="test.transformer"
                    ),
                    total=(
                        MatrixIntensityScaleState.log2(
                            established_by="test.transformer"
                        )
                        if total is not None
                        else None
                    ),
                ),
            )

    resolver = DatasetIntensityScaleResolver(transformer=MismatchedKindTransformer())
    with pytest.raises(
        TransformationStateEstablishmentError,
        match="configured transformer produced an invalid intensity scale state",
    ):
        resolver.run(
            phospho=_phospho(),
            total=_total(),
        )


def test_resolver_rejects_transformer_that_changes_total_matrix_presence() -> None:
    class DropsTotalTransformer:
        def run(
            self,
            phospho: pd.DataFrame,
            total: pd.DataFrame | None = None,
        ) -> TransformationResult:
            return TransformationResult(
                phospho=phospho,
                total=None,
                state=IntensityScaleState.raw(has_total_matrix=False),
            )

    resolver = DatasetIntensityScaleResolver(transformer=DropsTotalTransformer())

    with pytest.raises(
        TransformationStateEstablishmentError,
        match="changed total-matrix presence",
    ):
        resolver.run(
            phospho=_phospho(),
            total=_total(),
        )


def test_resolver_translates_unexpected_transformer_errors() -> None:
    class FailingTransformer:
        def run(
            self,
            phospho: pd.DataFrame,
            total: pd.DataFrame | None = None,
        ) -> TransformationResult:
            raise RuntimeError("boom")

    resolver = DatasetIntensityScaleResolver(transformer=FailingTransformer())

    with pytest.raises(
        TransformerExecutionError, match="configured transformer failed"
    ):
        resolver.run(
            phospho=_phospho(),
            total=None,
        )


def test_dataset_boundary_rejects_declared_intensity_scale_state_bypass() -> None:
    intensity_scale_state = IntensityScaleState.raw(has_total_matrix=False)
    with pytest.raises(
        TransformationValidationError,
        match="must be established through a supported PhosPy path",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            organism=Organism.RAT,
            intensity_scale_state=intensity_scale_state,
            processing_state=_processing_state_for(intensity_scale_state),
        )


def test_dataset_boundary_accepts_supported_established_state() -> None:
    supported_state = (
        DatasetIntensityScaleResolver(transformer=IdentityTransformer())
        .run(
            phospho=_phospho(),
            total=None,
        )
        .intensity_scale_state
    )
    dataset = AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        intensity_scale_state=supported_state,
        processing_state=_processing_state_for(supported_state),
    )
    assert dataset.intensity_scale_state.is_established


def test_direct_mint_established_raw_is_rejected() -> None:
    with pytest.raises(
        InvalidTransformationStateError,
        match="can be established only through supported PhosPy",
    ):
        IntensityScaleState.established_raw(has_total_matrix=False)


def test_direct_establishment_function_call_is_rejected() -> None:
    with pytest.raises(
        InvalidTransformationStateError,
        match="can be established only through supported PhosPy",
    ):
        establish_intensity_scale_state(
            IntensityScaleState.raw(has_total_matrix=False),
            established_via="phospy.datasets.builders.transformation_resolver",
        )


def test_fake_authority_object_is_rejected_even_with_supported_source() -> None:
    with pytest.raises(
        InvalidTransformationStateError,
        match="can be established only through supported PhosPy",
    ):
        establish_intensity_scale_state(
            IntensityScaleState.raw(has_total_matrix=False),
            established_via="phospy.datasets.builders.transformation_resolver",
            _authority=object(),
        )


def test_dataset_boundary_distinguishes_declared_from_supported_state() -> None:
    declared = IntensityScaleState.raw(has_total_matrix=False)
    supported = supported_linear_intensity_scale_state(has_total_matrix=False)
    assert not declared.is_established
    assert supported.is_established


def test_identity_transformer_is_strict_passthrough_establisher() -> None:
    phospho = _phospho().copy(deep=True)
    total = _total().copy(deep=True)

    result = IdentityTransformer().run(phospho=phospho, total=total)

    pdt.assert_frame_equal(result.phospho, phospho)
    pdt.assert_frame_equal(result.total, total)
    assert result.state.is_established
    assert result.state.kind.value == "linear"


def test_resolver_can_establish_log2_state_for_identity_transformer_when_expected() -> (
    None
):
    resolver = DatasetIntensityScaleResolver(transformer=IdentityTransformer())
    resolved = resolver.run(
        phospho=_phospho(),
        total=_total(),
        expected_scale_kind=IntensityScaleKind.LOG2,
    )
    assert resolved.intensity_scale_state.phospho.kind.value == "log2"
    assert resolved.intensity_scale_state.total is not None
    assert resolved.intensity_scale_state.total.kind.value == "log2"


def test_bundle_reconstruction_lane_establishes_state() -> None:
    state = intensity_scale_state_from_payload(
        {
            "phospho": {
                "kind": "linear",
                "transformed": False,
                "established_by": "bundle.fixture",
            },
            "total": None,
            "quantity": "phosphosite_abundance",
        }
    )

    assert state.is_established
    assert state.established_via == "phospy.io.bundles._shared.intensity_scale_state"


def test_processing_state_payload_rejects_minimal_total_correction_without_explicit_quantitative_meaning() -> (
    None
):
    payload = {
        "intensity_scale": {
            "phospho": {
                "kind": "log2",
                "transformed": True,
                "established_by": "bundle.fixture",
            },
            "total": {
                "kind": "log2",
                "transformed": True,
                "established_by": "bundle.fixture",
            },
            "quantity": "phospho_total_log_ratio",
        },
        "missing_data": {
            "policy": "forbid",
            "min_observed_values": None,
            "complete_matrix": True,
            "imputed": False,
        },
        "normalisation": {"policy": "none"},
        "total_protein_correction": {
            "policy": "subtract_log_total",
            "applied": True,
            "requires_log_scale": True,
            "diagnostics": {
                "diagnostics_schema_version": 1,
                "policy": "subtract_log_total",
                "requested_policy": "subtract_log_total",
                "resolved_policy": "subtract_log_total",
                "quantitative_meaning": "phospho_total_log_ratio",
            },
        },
        "site_matrix": {
            "policy": "as_input",
            "constructed": False,
            "missing_data_policy": "drop_any_missing",
            "minimum_observed_values": None,
            "duplicate_site_policy": "max_mean_signal",
        },
        "comparisons": {
            "policy": "none",
            "sample_group_column": "comparison_group",
            "pairs": None,
        },
    }

    with pytest.raises(
        PhosPyInputError,
        match=(
            "dataset.metadata.processing_state.total_protein_correction."
            "quantitative_meaning is required"
        ),
    ):
        processing_state_from_payload(payload)


def test_processing_state_payload_rejects_applied_total_correction_without_versioned_diagnostics() -> (
    None
):
    payload = {
        "intensity_scale": {
            "phospho": {
                "kind": "log2",
                "transformed": True,
                "established_by": "bundle.fixture",
            },
            "total": {
                "kind": "log2",
                "transformed": True,
                "established_by": "bundle.fixture",
            },
            "quantity": "phospho_total_log_ratio",
        },
        "missing_data": {
            "policy": "forbid",
            "min_observed_values": None,
            "complete_matrix": True,
            "imputed": False,
        },
        "normalisation": {"policy": "none"},
        "total_protein_correction": {
            "policy": "subtract_log_total",
            "applied": True,
            "requires_log_scale": True,
            "quantitative_meaning": "phospho_total_log_ratio",
            "diagnostics": None,
        },
        "site_matrix": {
            "policy": "as_input",
            "constructed": False,
            "missing_data_policy": "drop_any_missing",
            "minimum_observed_values": None,
            "duplicate_site_policy": "max_mean_signal",
        },
        "comparisons": {
            "policy": "none",
            "sample_group_column": "comparison_group",
            "pairs": None,
        },
    }

    with pytest.raises(
        PhosPyInputError,
        match=(
            "dataset.metadata.processing_state.total_protein_correction.diagnostics "
            "must be an object with"
        ),
    ):
        processing_state_from_payload(payload)
