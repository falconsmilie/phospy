"""Internal transformation-state establishment for dataset builder execution."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.transformations import (
    PhosPyTransformationError,
    TransformationStateEstablishmentError,
    TransformerExecutionError,
)
from phospy.transformations.contracts import Transformer
from phospy.transformations.models import (
    TransformationState,
    establish_transformation_state,
)
from phospy.validation.transformations.state import TransformationStateValidator


@dataclass(frozen=True, slots=True)
class ResolvedTransformation:
    """Quantitative matrices paired with an established transformation state."""

    phospho: pd.DataFrame
    total: pd.DataFrame | None
    transformation_state: TransformationState


class DatasetTransformationResolver:
    """Resolve transformation state for a dataset build request.

    Supported establishment path:
    A configured transformer that establishes state from matrices.

    Public builder policy wires this resolver with identity pass-through
    establishment by default, keeping transformation behavior narrow and honest.
    """

    def __init__(self, *, transformer: Transformer | None = None) -> None:
        self._transformer = transformer
        self._state_validator = TransformationStateValidator()

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        total: pd.DataFrame | None,
    ) -> ResolvedTransformation:
        if self._transformer is None:
            raise TransformationStateEstablishmentError(
                "unable to establish transformation state with confidence: "
                "no supported transformation establisher is configured. "
                "Configure AnalysisReadyDatasetBuilder("
                "executor=DatasetBuildExecutor(transformer=...))."
            )

        try:
            transformed = self._transformer.run(phospho=phospho, total=total)
        except PhosPyTransformationError:
            raise
        except Exception as exc:
            raise TransformerExecutionError(
                "configured transformer failed while establishing transformation "
                "state from dataset matrices"
            ) from exc

        if (total is None) is not (transformed.total is None):
            raise TransformationStateEstablishmentError(
                "configured transformer changed total-matrix presence while "
                "establishing transformation state; this is unsupported. "
                "Use a transformer that preserves phospho/total matrix presence."
            )

        self._validate_state(
            transformed.state,
            has_total_matrix=transformed.total is not None,
            source="configured transformer",
        )
        transformer_source = (
            f"{self._transformer.__class__.__module__}."
            f"{self._transformer.__class__.__qualname__}"
        )
        established_state = establish_transformation_state(
            transformed.state,
            established_via=transformer_source,
        )
        return ResolvedTransformation(
            phospho=transformed.phospho,
            total=transformed.total,
            transformation_state=established_state,
        )

    def _validate_state(
        self,
        state: TransformationState,
        *,
        has_total_matrix: bool,
        source: str,
    ) -> None:
        try:
            self._state_validator.run(
                transformation_state=state,
                has_total_matrix=has_total_matrix,
            )
        except Exception as exc:
            raise TransformationStateEstablishmentError(
                f"{source} produced an invalid transformation state: {exc}"
            ) from exc
