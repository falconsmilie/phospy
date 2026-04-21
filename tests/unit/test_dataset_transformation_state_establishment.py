from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.datasets.builders.transformation_resolver import (
    DatasetTransformationResolver,
)
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.transformations import (
    InvalidTransformationStateError,
    TransformationStateEstablishmentError,
    TransformerExecutionError,
)
from phospy.errors.validation import TransformationValidationError
from phospy.io.bundles._shared.transformation_state import (
    transformation_state_from_payload,
)
from phospy.references.models import Organism
from phospy.transformations.contracts import TransformationResult
from phospy.transformations.models import (
    MatrixTransformationState,
    TransformationState,
    establish_transformation_state,
)
from phospy.transformations.transformers import IdentityTransformer
from tests.support.transformation_states import supported_linear_state


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


def test_resolver_fails_when_state_is_unknown_and_no_transformer_is_configured() -> (
    None
):
    resolver = DatasetTransformationResolver()

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
                state=TransformationState(
                    phospho=MatrixTransformationState.log2(
                        established_by="test.transformer"
                    ),
                    total=(
                        MatrixTransformationState.log2(
                            established_by="test.transformer"
                        )
                        if total is not None
                        else None
                    ),
                ),
            )

    resolver = DatasetTransformationResolver(transformer=DeclaredLog2Transformer())

    resolved = resolver.run(
        phospho=_phospho(),
        total=_total(),
    )

    assert resolved.transformation_state.phospho.kind.value == "log2"
    assert resolved.transformation_state.total is not None
    assert resolved.transformation_state.total.kind.value == "log2"
    assert resolved.transformation_state.is_established
    assert resolved.transformation_state.established_via is not None
    assert "DeclaredLog2Transformer" in resolved.transformation_state.established_via


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
                state=TransformationState(
                    phospho=MatrixTransformationState.linear(
                        established_by="test.transformer"
                    ),
                    total=(
                        MatrixTransformationState.log2(
                            established_by="test.transformer"
                        )
                        if total is not None
                        else None
                    ),
                ),
            )

    resolver = DatasetTransformationResolver(transformer=MismatchedKindTransformer())
    with pytest.raises(
        TransformationStateEstablishmentError,
        match="configured transformer produced an invalid transformation state",
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
                state=TransformationState.raw(has_total_matrix=False),
            )

    resolver = DatasetTransformationResolver(transformer=DropsTotalTransformer())

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

    resolver = DatasetTransformationResolver(transformer=FailingTransformer())

    with pytest.raises(
        TransformerExecutionError, match="configured transformer failed"
    ):
        resolver.run(
            phospho=_phospho(),
            total=None,
        )


def test_dataset_boundary_rejects_declared_transformation_state_bypass() -> None:
    with pytest.raises(
        TransformationValidationError,
        match="must be established through a supported PhosPy path",
    ):
        AnalysisReadyPhosphoDataset(
            phospho=_phospho(),
            site_metadata=_site_metadata(),
            organism=Organism.RAT,
            transformation_state=TransformationState.raw(has_total_matrix=False),
        )


def test_dataset_boundary_accepts_supported_established_state() -> None:
    supported_state = (
        DatasetTransformationResolver(transformer=IdentityTransformer())
        .run(
            phospho=_phospho(),
            total=None,
        )
        .transformation_state
    )
    dataset = AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        transformation_state=supported_state,
    )
    assert dataset.transformation_state.is_established


def test_direct_mint_established_raw_is_rejected() -> None:
    with pytest.raises(
        InvalidTransformationStateError,
        match="can be established only through supported PhosPy",
    ):
        TransformationState.established_raw(has_total_matrix=False)


def test_direct_establishment_function_call_is_rejected() -> None:
    with pytest.raises(
        InvalidTransformationStateError,
        match="can be established only through supported PhosPy",
    ):
        establish_transformation_state(
            TransformationState.raw(has_total_matrix=False),
            established_via="phospy.datasets.builders.transformation_resolver",
        )


def test_fake_authority_object_is_rejected_even_with_supported_source() -> None:
    with pytest.raises(
        InvalidTransformationStateError,
        match="can be established only through supported PhosPy",
    ):
        establish_transformation_state(
            TransformationState.raw(has_total_matrix=False),
            established_via="phospy.datasets.builders.transformation_resolver",
            _authority=object(),
        )


def test_dataset_boundary_distinguishes_declared_from_supported_state() -> None:
    declared = TransformationState.raw(has_total_matrix=False)
    supported = supported_linear_state(has_total_matrix=False)
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


def test_bundle_reconstruction_lane_establishes_state() -> None:
    state = transformation_state_from_payload(
        {
            "phospho": {
                "kind": "linear",
                "transformed": False,
                "established_by": "bundle.fixture",
            },
            "total": None,
        }
    )

    assert state.is_established
    assert state.established_via == "phospy.io.bundles._shared.transformation_state"
