"""Shared quantitative-meaning contracts for workflow inputs."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import cast

from phospy.errors.validation import (
    TransformationValidationError,
    WorkflowValidationError,
)
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.transformations.models import (
    IntensityScaleKind,
    IntensityScaleState,
    QuantitativeMeaning,
)
from phospy.validation.transformations.state import IntensityScaleStateValidator

_MIXED_TOTAL_PROTEIN_QUANTITATIVE_MEANING = (
    QuantitativeMeaning.MIXED_PHOSPHO_TOTAL_LOG_RATIO_AND_PHOSPHOSITE_LOG_ABUNDANCE
)

PHOSPHOSITE_ABUNDANCE_INPUT_MEANINGS = frozenset(
    {
        QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
        QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO,
    }
)
MIXED_TOTAL_PROTEIN_INPUT_MEANINGS = frozenset(
    {_MIXED_TOTAL_PROTEIN_QUANTITATIVE_MEANING}
)


@dataclass(frozen=True, slots=True)
class WorkflowQuantitativeInputContract:
    """Allowed quantitative state for a workflow input matrix."""

    allowed_meanings: frozenset[QuantitativeMeaning]
    allowed_scales: frozenset[IntensityScaleKind]
    require_established_scale: bool = True
    allow_unknown_meaning: bool = False
    allow_mixed_meaning: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "allowed_meanings",
            frozenset(_normalize_quantitative_meanings(self.allowed_meanings)),
        )
        object.__setattr__(
            self,
            "allowed_scales",
            frozenset(_normalize_intensity_scale_kinds(self.allowed_scales)),
        )
        if not self.allowed_meanings:
            raise ValueError(
                "WorkflowQuantitativeInputContract.allowed_meanings must not be empty"
            )
        if not self.allowed_scales:
            raise ValueError(
                "WorkflowQuantitativeInputContract.allowed_scales must not be empty"
            )
        if not isinstance(cast(object, self.require_established_scale), bool):
            raise ValueError(
                "WorkflowQuantitativeInputContract.require_established_scale must be a bool"
            )
        if not isinstance(cast(object, self.allow_unknown_meaning), bool):
            raise ValueError(
                "WorkflowQuantitativeInputContract.allow_unknown_meaning must be a bool"
            )
        if not isinstance(cast(object, self.allow_mixed_meaning), bool):
            raise ValueError(
                "WorkflowQuantitativeInputContract.allow_mixed_meaning must be a bool"
            )


class WorkflowQuantitativeInputValidator:
    """Validate workflow input quantitative semantics against a shared contract."""

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
        contract: WorkflowQuantitativeInputContract,
        context: str,
        dataset_view: DatasetInternalView | None = None,
    ) -> IntensityScaleState:
        if not isinstance(cast(object, dataset), AnalysisReadyPhosphoDataset):
            raise WorkflowValidationError(
                f"{context} dataset must be AnalysisReadyPhosphoDataset"
            )
        state = dataset.intensity_scale_state
        resolved_dataset_view = dataset_view or DatasetInternalView(dataset)
        has_total_matrix = resolved_dataset_view.total is not None
        try:
            self._intensity_scale_state_validator.run(
                intensity_scale_state=state,
                has_total_matrix=has_total_matrix,
                require_established=contract.require_established_scale,
            )
        except TransformationValidationError as exc:
            raise WorkflowValidationError(
                f"{context} quantitative input state is invalid: {exc}"
            ) from exc

        quantity = state.quantity
        if quantity is None:
            raise WorkflowValidationError(
                f"{context} quantitative input state is missing quantitative meaning"
            )
        allowed = _format_quantitative_meanings(contract.allowed_meanings)
        if (
            quantity is QuantitativeMeaning.UNKNOWN
            and not contract.allow_unknown_meaning
        ):
            raise WorkflowValidationError(
                f"{context} requires quantitative meaning in {{{allowed}}}; "
                f"got {quantity.value!r}. "
                "unknown quantitative meaning is rejected by default; "
                "provide an explicit workflow override before using unknown matrix "
                "semantics"
            )
        if is_mixed_quantitative_meaning(quantity) and not contract.allow_mixed_meaning:
            raise WorkflowValidationError(
                f"{context} requires quantitative meaning in {{{allowed}}}; "
                f"got {quantity.value!r}. "
                "Mixed quantitative matrix semantics are not allowed unless the "
                "workflow explicitly enables them"
            )
        if quantity not in contract.allowed_meanings:
            raise WorkflowValidationError(
                f"{context} requires quantitative meaning in {{{allowed}}}; "
                f"got {quantity.value!r}"
            )
        scale = state.phospho.kind
        if scale not in contract.allowed_scales:
            allowed_scales = _format_intensity_scales(contract.allowed_scales)
            raise WorkflowValidationError(
                f"{context} requires phospho intensity scale in {{{allowed_scales}}}; "
                f"got {scale.value!r}"
            )
        return state


def phosphosite_abundance_workflow_input_contract(
    *,
    allow_mixed_total_protein_quantitative_meaning: bool = False,
) -> WorkflowQuantitativeInputContract:
    """Return the shared phosphosite-abundance workflow input contract."""

    allowed_meanings = set(PHOSPHOSITE_ABUNDANCE_INPUT_MEANINGS)
    if allow_mixed_total_protein_quantitative_meaning:
        allowed_meanings.update(MIXED_TOTAL_PROTEIN_INPUT_MEANINGS)
    return WorkflowQuantitativeInputContract(
        allowed_meanings=frozenset(allowed_meanings),
        allowed_scales=frozenset({IntensityScaleKind.LINEAR, IntensityScaleKind.LOG2}),
        require_established_scale=True,
        allow_mixed_meaning=allow_mixed_total_protein_quantitative_meaning,
    )


def kinase_profile_scoring_workflow_input_contract(
    *,
    allow_mixed_total_protein_quantitative_meaning: bool = False,
) -> WorkflowQuantitativeInputContract:
    """Return the quantitative contract for kinase profile scoring inputs."""

    return phosphosite_abundance_workflow_input_contract(
        allow_mixed_total_protein_quantitative_meaning=(
            allow_mixed_total_protein_quantitative_meaning
        )
    )


def signalome_workflow_input_contract(
    *,
    allow_mixed_total_protein_quantitative_meaning: bool = False,
) -> WorkflowQuantitativeInputContract:
    """Return the quantitative contract for signalome upstream dataset inputs."""

    return phosphosite_abundance_workflow_input_contract(
        allow_mixed_total_protein_quantitative_meaning=(
            allow_mixed_total_protein_quantitative_meaning
        )
    )


def is_mixed_quantitative_meaning(quantity: QuantitativeMeaning | None) -> bool:
    return quantity is _MIXED_TOTAL_PROTEIN_QUANTITATIVE_MEANING


def _normalize_quantitative_meanings(
    values: Iterable[QuantitativeMeaning | str],
) -> tuple[QuantitativeMeaning, ...]:
    resolved: list[QuantitativeMeaning] = []
    for value in values:
        if isinstance(value, QuantitativeMeaning):
            resolved.append(value)
            continue
        resolved.append(QuantitativeMeaning(str(value)))
    return tuple(resolved)


def _normalize_intensity_scale_kinds(
    values: Iterable[IntensityScaleKind | str],
) -> tuple[IntensityScaleKind, ...]:
    resolved: list[IntensityScaleKind] = []
    for value in values:
        if isinstance(value, IntensityScaleKind):
            resolved.append(value)
            continue
        resolved.append(IntensityScaleKind(str(value)))
    return tuple(resolved)


def _format_quantitative_meanings(values: frozenset[QuantitativeMeaning]) -> str:
    return ", ".join(
        repr(value.value) for value in sorted(values, key=lambda item: item.value)
    )


def _format_intensity_scales(values: frozenset[IntensityScaleKind]) -> str:
    return ", ".join(
        repr(value.value) for value in sorted(values, key=lambda item: item.value)
    )


DIFFERENTIAL_LOG_ABUNDANCE_INPUT_CONTRACT = WorkflowQuantitativeInputContract(
    allowed_meanings=frozenset({QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE}),
    allowed_scales=frozenset({IntensityScaleKind.LOG2}),
    require_established_scale=True,
)


__all__ = [
    "DIFFERENTIAL_LOG_ABUNDANCE_INPUT_CONTRACT",
    "PHOSPHOSITE_ABUNDANCE_INPUT_MEANINGS",
    "WorkflowQuantitativeInputContract",
    "WorkflowQuantitativeInputValidator",
    "is_mixed_quantitative_meaning",
    "kinase_profile_scoring_workflow_input_contract",
    "phosphosite_abundance_workflow_input_contract",
    "signalome_workflow_input_contract",
]
