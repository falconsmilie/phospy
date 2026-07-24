"""Intensity-transform stage for dataset preprocessing."""

from __future__ import annotations

from dataclasses import replace

from phospy.errors.input import PhosPyInputError
from phospy.provenance.hashing import hash_table_tolerance
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
    PreprocessingStateTableKey,
)
from phospy.science.datasets.preprocessing.policy_models import IntensityTransformPolicy
from phospy.science.datasets.preprocessing.stage_contract import (
    DeterminismKind,
    PreprocessingStageContract,
    PreprocessingStageFactoryContext,
)
from phospy.science.transformations.contracts import Transformer
from phospy.science.transformations.transformers import (
    IdentityTransformer,
    Log2Transformer,
)


class IntensityTransformStage:
    """Apply configured quantitative intensity transform to phospho values."""

    stage_key = DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM

    def _resolve_transformer(self, state: PreprocessingState) -> Transformer:
        policy = state.plan.intensity_transform_policy
        if policy is IntensityTransformPolicy.IDENTITY:
            return IdentityTransformer()
        if policy is IntensityTransformPolicy.LOG2:
            return Log2Transformer(
                pseudocount=float(state.plan.intensity_transform_pseudocount)
            )
        raise PhosPyInputError(
            "dataset build request preprocessing_config contains an unsupported "
            "intensity_transform.policy"
        )

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        transformed = self._resolve_transformer(state).run(
            phospho=state.phospho,
            total=state.total,
        )
        transformed_phospho = transformed.phospho
        transformed_total = transformed.total
        next_state = replace(
            state,
            phospho=transformed_phospho,
            total=transformed_total,
        )
        input_phospho_hash = hash_table_tolerance(
            state.phospho,
            name="intensity_transform.input.phospho",
        )
        output_phospho_hash = hash_table_tolerance(
            transformed_phospho,
            name="intensity_transform.output.phospho",
        )
        intensity_transformation_event = transformed.intensity_transformation_event
        if intensity_transformation_event is not None:
            intensity_transformation_event = replace(
                intensity_transformation_event,
                input_fingerprint=input_phospho_hash,
                output_fingerprint=output_phospho_hash,
            )
        diagnostics = {
            **dict(transformed.provenance),
            "input_phospho_hash": input_phospho_hash,
            "output_phospho_hash": output_phospho_hash,
        }
        diagnostics.setdefault(
            "output_intensity_scale_kind",
            transformed.state.kind.value,
        )
        diagnostics.setdefault(
            "policy",
            state.plan.intensity_transform_policy.value,
        )
        diagnostics.setdefault(
            "pseudocount",
            float(state.plan.intensity_transform_pseudocount),
        )
        diagnostics.setdefault(
            "affected_matrices",
            ["phospho"] if transformed_total is None else ["phospho", "total"],
        )
        if state.total is not None:
            diagnostics["input_total_hash"] = hash_table_tolerance(
                state.total,
                name="intensity_transform.input.total",
            )
        if transformed_total is not None:
            diagnostics["output_total_hash"] = hash_table_tolerance(
                transformed_total,
                name="intensity_transform.output.total",
            )
        return PreprocessingStageResult(
            state=next_state,
            diagnostics={
                "dropped_row_ids": (),
                "dropped_row_count": 0,
                "imputed_cell_count": 0,
                "imputed_row_ids": (),
                "notes": "stage executed",
                "diagnostics": diagnostics,
            },
            intensity_transformation_event=intensity_transformation_event,
        )


def _resolve_operation(plan: PreprocessingPlan) -> str:
    return plan.intensity_transform_policy.value


def _resolve_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    return {"pseudocount": float(plan.intensity_transform_pseudocount)}


def _build_intensity_transform_stage(
    _context: PreprocessingStageFactoryContext,
) -> IntensityTransformStage:
    return IntensityTransformStage()


INTENSITY_TRANSFORM_STAGE_CONTRACT = PreprocessingStageContract(
    stage_key=DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    display_label=DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    provenance_stage=DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    operation_name=_resolve_operation,
    serialize_parameters=_resolve_parameters,
    consumed_input_tables=(
        PreprocessingStateTableKey.DATASET_PHOSPHO,
        PreprocessingStateTableKey.DATASET_TOTAL,
    ),
    produced_output_tables=(
        PreprocessingStateTableKey.DATASET_PHOSPHO,
        PreprocessingStateTableKey.DATASET_TOTAL,
    ),
    stage_factory=_build_intensity_transform_stage,
    backend="numpy",
    determinism_kind=DeterminismKind.DETERMINISTIC,
    diagnostics_metadata={
        "known_diagnostics_fields": (
            "policy",
            "pseudocount",
            "output_intensity_scale_kind",
            "affected_matrices",
            "input_phospho_hash",
            "output_phospho_hash",
            "input_total_hash",
            "output_total_hash",
        )
    },
)


__all__ = ["INTENSITY_TRANSFORM_STAGE_CONTRACT", "IntensityTransformStage"]
