"""Shared primitives for method-owned quantitative input contracts."""

from __future__ import annotations

from typing import cast

from phospy.errors.validation import (
    TransformationValidationError,
    WorkflowValidationError,
)
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.quantitative_method_contracts import (
    MethodQuantitativeInputContract,
    ResolvedMethodQuantitativeInputContract,
)
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
)
from phospy.validation.transformations.state import IntensityScaleStateValidator
from phospy.validation.workflows.quantitative import is_mixed_quantitative_meaning


class MethodQuantitativeInputValidator:
    """Validate a dataset against a method-supplied quantitative contract."""

    def __init__(
        self,
        *,
        intensity_scale_state_validator: IntensityScaleStateValidator | None = None,
    ) -> None:
        self._intensity_scale_state_validator = (
            intensity_scale_state_validator or IntensityScaleStateValidator()
        )

    def run(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        contract: MethodQuantitativeInputContract,
        context: str,
    ) -> ResolvedMethodQuantitativeInputContract:
        if not isinstance(cast(object, dataset), AnalysisReadyPhosphoDataset):
            raise WorkflowValidationError(
                f"{context} dataset must be AnalysisReadyPhosphoDataset"
            )
        if not isinstance(contract, MethodQuantitativeInputContract):
            raise WorkflowValidationError(
                f"{context} method quantitative contract must be "
                "MethodQuantitativeInputContract"
            )
        state = dataset.intensity_scale_state
        has_total_matrix = DatasetInternalView(dataset).total is not None
        try:
            self._intensity_scale_state_validator.run(
                intensity_scale_state=state,
                has_total_matrix=has_total_matrix,
                require_established=True,
            )
        except TransformationValidationError as exc:
            raise WorkflowValidationError(
                f"{context} quantitative input state is invalid: {exc}"
            ) from exc

        if not contract.quantitative_input_required:
            return ResolvedMethodQuantitativeInputContract(
                contract=contract,
                resolved_scale=state.phospho.kind,
                resolved_meaning=state.quantity,
                enforcement_context=context,
            )
        quantity = state.quantity
        if quantity is None:
            raise WorkflowValidationError(
                f"{context} requires quantitative meaning for method "
                f"{contract.method_id!r}; got None"
            )
        accepted_meanings = cast(
            tuple[QuantitativeMeaning, ...],
            contract.accepted_meanings,
        )
        if (
            quantity is QuantitativeMeaning.UNKNOWN
            and quantity not in accepted_meanings
        ):
            raise WorkflowValidationError(
                f"{context} requires quantitative meaning in "
                f"{{{_format_meanings(accepted_meanings)}}} for method "
                f"{contract.method_id!r}; got {quantity.value!r}. unknown "
                "quantitative meaning is rejected "
                "unless the method contract explicitly accepts it"
            )
        if (
            is_mixed_quantitative_meaning(quantity)
            and quantity not in accepted_meanings
        ):
            raise WorkflowValidationError(
                f"{context} method {contract.method_id!r} requires "
                "quantitative meaning in "
                f"{{{_format_meanings(accepted_meanings)}}}; got "
                f"{quantity.value!r}. mixed "
                "quantitative matrix semantics are "
                "not accepted by this method contract"
            )
        if quantity not in accepted_meanings:
            raise WorkflowValidationError(
                f"{context} requires quantitative meaning in "
                f"{{{_format_meanings(accepted_meanings)}}} for method "
                f"{contract.method_id!r}; got {quantity.value!r}"
            )
        accepted_scales = cast(
            tuple[IntensityScaleKind, ...],
            contract.accepted_scales,
        )
        scale = state.phospho.kind
        if scale not in accepted_scales:
            raise WorkflowValidationError(
                f"{context} requires phospho intensity scale in "
                f"{{{_format_scales(accepted_scales)}}} for method "
                f"{contract.method_id!r}; got {scale.value!r}. Linear and log2 "
                "inputs are not treated as equivalent and no implicit "
                "transformation is applied"
            )
        return ResolvedMethodQuantitativeInputContract(
            contract=contract,
            resolved_scale=scale,
            resolved_meaning=quantity,
            enforcement_context=context,
        )


def _format_meanings(values: tuple[QuantitativeMeaning, ...]) -> str:
    return ", ".join(
        repr(value.value) for value in sorted(values, key=lambda x: x.value)
    )


def _format_scales(values: tuple[IntensityScaleKind, ...]) -> str:
    return ", ".join(
        repr(value.value) for value in sorted(values, key=lambda x: x.value)
    )


__all__ = ["MethodQuantitativeInputValidator"]
