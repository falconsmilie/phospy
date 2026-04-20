"""Transformation-state validator."""

from __future__ import annotations

from phospy.errors.validation import TransformationValidationError
from phospy.transformations.models import TransformationState


class TransformationStateValidator:
    """Validate explicit transformation state coherence."""

    def run(
        self,
        transformation_state: TransformationState,
        *,
        has_total_matrix: bool,
        require_established: bool = False,
    ) -> TransformationState:
        if not isinstance(transformation_state, TransformationState):
            raise TransformationValidationError(
                "dataset.transformation_state must be a TransformationState instance"
            )
        if require_established and not transformation_state.is_established:
            raise TransformationValidationError(
                "dataset.transformation_state must be established through a "
                "supported PhosPy path; use AnalysisReadyDatasetBuilder or a "
                "supported transformer/bundle reconstruction path"
            )
        if has_total_matrix and transformation_state.total is None:
            raise TransformationValidationError(
                "transformation_state.total is required when dataset.total is provided"
            )
        if not has_total_matrix and transformation_state.total is not None:
            raise TransformationValidationError(
                "transformation_state.total must be None when dataset.total is absent"
            )
        if (
            transformation_state.total is not None
            and transformation_state.total.kind is not transformation_state.phospho.kind
        ):
            raise TransformationValidationError(
                "phospho and total matrices must share one transformation kind"
            )
        return transformation_state
