from __future__ import annotations

import pandas as pd
import pytest

from phospy.datasets.builders.transformation_resolver import (
    DatasetTransformationResolver,
)
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.transformations import (
    TransformationStateEstablishmentError,
    TransformerExecutionError,
)
from phospy.errors.validation import TransformationValidationError
from phospy.references.models import Organism
from phospy.transformations.contracts import TransformationResult
from phospy.transformations.models import MatrixTransformationState, TransformationState


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
    dataset = AnalysisReadyPhosphoDataset(
        phospho=_phospho(),
        site_metadata=_site_metadata(),
        organism=Organism.RAT,
        transformation_state=TransformationState.established_raw(
            has_total_matrix=False
        ),
    )
    assert dataset.transformation_state.is_established
