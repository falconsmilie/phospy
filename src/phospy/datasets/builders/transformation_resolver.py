"""Internal intensity-scale-state establishment for dataset builder execution."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.transformations import (
    PhosPyTransformationError,
    TransformationStateEstablishmentError,
    TransformerExecutionError,
)
from phospy.errors.validation import TransformationValidationError
from phospy.transformations._authority import (
    _dataset_resolver_establishment_authority,
)
from phospy.transformations.contracts import Transformer
from phospy.transformations.models import (
    IntensityScaleKind,
    IntensityScaleState,
    MatrixIntensityScaleState,
    establish_intensity_scale_state,
)
from phospy.transformations.transformers import IdentityTransformer
from phospy.validation.transformations.state import IntensityScaleStateValidator


@dataclass(frozen=True, slots=True)
class ResolvedIntensityScale:
    """Quantitative matrices paired with an established intensity scale state."""

    phospho: pd.DataFrame
    total: pd.DataFrame | None
    intensity_scale_state: IntensityScaleState


class DatasetIntensityScaleResolver:
    """Resolve intensity scale state for a dataset build request.

    Supported establishment path:
    A configured transformer that establishes state from matrices.

    Public builder policy wires this resolver with identity pass-through
    establishment by default, keeping transformation behavior narrow and honest.
    """

    def __init__(self, *, transformer: Transformer | None = None) -> None:
        self._transformer = transformer
        self._state_validator = IntensityScaleStateValidator()

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        total: pd.DataFrame | None,
        expected_scale_kind: IntensityScaleKind | None = None,
    ) -> ResolvedIntensityScale:
        if self._transformer is None:
            raise TransformationStateEstablishmentError(
                "unable to establish intensity scale state with confidence: "
                "no supported intensity-scale establisher is configured. "
                "Configure AnalysisReadyDatasetBuilder("
                "executor=DatasetBuildExecutor(transformer=...))."
            )

        try:
            transformed = self._transformer.run(phospho=phospho, total=total)
        except PhosPyTransformationError:
            raise
        except (TypeError, ValueError) as exc:
            raise TransformerExecutionError(
                "configured transformer failed while establishing intensity scale "
                "state from dataset matrices"
            ) from exc

        if (total is None) is not (transformed.total is None):
            raise TransformationStateEstablishmentError(
                "configured transformer changed total-matrix presence while "
                "establishing intensity scale state; this is unsupported. "
                "Use a transformer that preserves phospho/total matrix presence."
            )

        state = transformed.state
        if (
            expected_scale_kind is not None
            and state.phospho.kind is not expected_scale_kind
            and isinstance(self._transformer, IdentityTransformer)
        ):
            state = _build_identity_state_for_scale_kind(
                expected_scale_kind=expected_scale_kind,
                has_total_matrix=transformed.total is not None,
            )
        if (
            expected_scale_kind is not None
            and state.phospho.kind is not expected_scale_kind
        ):
            raise TransformationStateEstablishmentError(
                "configured transformer produced an intensity scale that is "
                "incompatible with the configured preprocessing intensity "
                f"transform policy; expected '{expected_scale_kind.value}' but "
                f"received '{state.phospho.kind.value}'"
            )

        self._validate_state(
            state,
            has_total_matrix=transformed.total is not None,
            source="configured transformer",
        )
        transformer_source = (
            f"{self._transformer.__class__.__module__}."
            f"{self._transformer.__class__.__qualname__}"
        )
        established_state = establish_intensity_scale_state(
            state,
            established_via=transformer_source,
            _authority=_dataset_resolver_establishment_authority(),
        )
        return ResolvedIntensityScale(
            phospho=transformed.phospho,
            total=transformed.total,
            intensity_scale_state=established_state,
        )

    def _validate_state(
        self,
        state: IntensityScaleState,
        *,
        has_total_matrix: bool,
        source: str,
    ) -> None:
        try:
            self._state_validator.run(
                intensity_scale_state=state,
                has_total_matrix=has_total_matrix,
            )
        except TransformationValidationError as exc:
            raise TransformationStateEstablishmentError(
                f"{source} produced an invalid intensity scale state: {exc}"
            ) from exc


def _build_identity_state_for_scale_kind(
    *,
    expected_scale_kind: IntensityScaleKind,
    has_total_matrix: bool,
) -> IntensityScaleState:
    if expected_scale_kind is IntensityScaleKind.LINEAR:
        return IntensityScaleState.raw(has_total_matrix=has_total_matrix)
    if expected_scale_kind is IntensityScaleKind.LOG2:
        phospho_state = MatrixIntensityScaleState.log2(
            established_by="phospy.transformations.transformers.identity"
        )
        if has_total_matrix:
            return IntensityScaleState(
                phospho=phospho_state,
                total=MatrixIntensityScaleState.log2(
                    established_by="phospy.transformations.transformers.identity"
                ),
            )
        return IntensityScaleState(phospho=phospho_state, total=None)
    raise TransformationStateEstablishmentError(
        f"unsupported expected intensity scale kind for resolver: {expected_scale_kind}"
    )
