from __future__ import annotations

from dataclasses import replace

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.stages.intensity_transform import (
    IntensityTransformStage,
)
from phospy.science.transformations.contracts import TransformationResult
from phospy.science.transformations.models import (
    IntensityScaleKind,
    IntensityScaleState,
)
from phospy.science.transformations.transformers import Log2Transformer


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 3.0],
            "sample_b": [2.0, 7.0],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
            "localisation_confidence": [0.95, 0.92],
        },
        index=index.copy(),
    )


def test_log2_transformer_transforms_expected_numeric_values() -> None:
    transformed = Log2Transformer(pseudocount=1.0).run(phospho=_phospho(), total=None)
    expected = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.584962500721156, 3.0],
        },
        index=_phospho().index.copy(),
    )

    pdt.assert_frame_equal(transformed.phospho, expected)
    assert transformed.total is None


@pytest.mark.parametrize(
    ("pseudocount", "expected"),
    [
        (float("nan"), "pseudocount must be finite"),
        (-1.0, "pseudocount must be greater than or equal to 0"),
    ],
)
def test_log2_transformer_rejects_invalid_pseudocount(
    pseudocount: float,
    expected: str,
) -> None:
    with pytest.raises(PhosPyInputError, match=expected):
        Log2Transformer(pseudocount=pseudocount)


def test_log2_transformer_returns_log2_intensity_scale_state() -> None:
    total = pd.DataFrame(
        {"sample_a": [3.0], "sample_b": [7.0]},
        index=pd.Index(["MAPK14"], name="protein_id"),
    )
    transformed = Log2Transformer(pseudocount=1.0).run(
        phospho=_phospho(),
        total=total,
    )

    assert transformed.state.phospho.kind is IntensityScaleKind.LOG2
    assert transformed.state.phospho.transformed is True
    assert transformed.state.total is not None
    assert transformed.state.total.kind is IntensityScaleKind.LOG2
    assert transformed.state.total.transformed is True
    assert transformed.provenance["pseudocount"] == 1.0
    assert transformed.provenance["affected_matrices"] == ["phospho", "total"]
    assert transformed.provenance["output_intensity_scale_kind"] == "log2"
    assert isinstance(transformed.provenance["transformer_state"], dict)


def test_log2_transformer_exposes_scale_capabilities() -> None:
    transformer = Log2Transformer(pseudocount=1.0)
    assert transformer.preserves_input_scale_state is False
    assert transformer.changes_numeric_values is True
    assert transformer.requires_established_input_state is False


def test_intensity_transform_stage_delegates_to_transformer() -> None:
    class SpyTransformer:
        def __init__(self) -> None:
            self.calls = 0
            self.last_phospho: pd.DataFrame | None = None
            self.last_total: pd.DataFrame | None = None

        def run(
            self,
            phospho: pd.DataFrame,
            total: pd.DataFrame | None = None,
        ) -> TransformationResult:
            self.calls += 1
            self.last_phospho = phospho
            self.last_total = total
            return TransformationResult(
                phospho=phospho + 5.0,
                total=total,
                state=IntensityScaleState.raw(has_total_matrix=total is not None),
                provenance={"policy": "identity", "transformer_name": "test.spy"},
            )

    class StageUnderTest(IntensityTransformStage):
        def __init__(self, transformer: SpyTransformer) -> None:
            self._transformer = transformer

        def _resolve_transformer(self, state: PreprocessingState):
            return self._transformer

    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=replace(
            PreprocessingPlan.default(),
            stage_order=("intensity_transform",),
        ),
    )
    spy = SpyTransformer()
    result = StageUnderTest(spy).run(state)

    assert spy.calls == 1
    assert spy.last_phospho is phospho
    assert spy.last_total is None
    pdt.assert_frame_equal(result.state.phospho, phospho + 5.0)
    diagnostics = result.diagnostics["diagnostics"]
    assert diagnostics["policy"] == "identity"
    assert diagnostics["transformer_name"] == "test.spy"
    assert isinstance(diagnostics["input_phospho_hash"], str)
    assert isinstance(diagnostics["output_phospho_hash"], str)


def test_intensity_transform_stage_uses_log2_transformer_for_log2_policy() -> None:
    phospho = _phospho()
    state = PreprocessingState(
        phospho=phospho,
        site_metadata=_site_metadata(phospho.index),
        sample_metadata=None,
        total=None,
        plan=replace(
            PreprocessingPlan.default(),
            intensity_transform_policy="log2",
            intensity_transform_pseudocount=2.0,
            stage_order=("intensity_transform",),
        ),
    )

    transformer = IntensityTransformStage()._resolve_transformer(state)

    assert isinstance(transformer, Log2Transformer)
    assert transformer.pseudocount == 2.0
