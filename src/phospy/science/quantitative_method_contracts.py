"""Method-owned quantitative input contract models.

These models describe the scientific quantitative input assumptions for scoring
or activity methods. They intentionally contain no workflow-specific policy
lookup tables; individual methods own those declarations.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import cast

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.activities.semantics import (
    ActivityInputMatrix,
    ActivityProfileAxis,
    ActivityQuantitativeSemantics,
)
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
)


@dataclass(frozen=True, slots=True)
class MethodQuantitativeInputContract:
    """Scientific quantitative input assumptions declared by one method."""

    method_id: str
    accepted_scales: tuple[IntensityScaleKind | str, ...]
    accepted_meanings: tuple[QuantitativeMeaning | str, ...]
    required_centring: str
    required_standardisation: str
    missing_value_treatment: str
    profile_axis_requirements: str
    statistical_interpretation: str
    accepted_activity_profile_axes: tuple[ActivityProfileAxis | str, ...] = ()
    accepted_activity_quantitative_semantics: tuple[
        ActivityQuantitativeSemantics | str,
        ...,
    ] = ()
    p_value_interpretation: str | None = None
    quantitative_input_required: bool = True
    scale_sensitivity: str = (
        "Scale-sensitive method; accepted scales are not interchangeable and "
        "the method does not transform input values."
    )
    no_implicit_transformation: bool = True
    contract_version: str = "1"

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "method_id",
            _require_non_empty_text(
                self.method_id,
                field_name="method_quantitative_input_contract.method_id",
            ),
        )
        object.__setattr__(
            self,
            "accepted_scales",
            _normalize_intensity_scales(self.accepted_scales),
        )
        object.__setattr__(
            self,
            "accepted_meanings",
            _normalize_quantitative_meanings(self.accepted_meanings),
        )
        object.__setattr__(
            self,
            "accepted_activity_profile_axes",
            _normalize_activity_profile_axes(self.accepted_activity_profile_axes),
        )
        object.__setattr__(
            self,
            "accepted_activity_quantitative_semantics",
            _normalize_activity_quantitative_semantics_values(
                self.accepted_activity_quantitative_semantics
            ),
        )
        for field_name in (
            "required_centring",
            "required_standardisation",
            "missing_value_treatment",
            "profile_axis_requirements",
            "statistical_interpretation",
            "scale_sensitivity",
            "contract_version",
        ):
            object.__setattr__(
                self,
                field_name,
                _require_non_empty_text(
                    getattr(self, field_name),
                    field_name=f"method_quantitative_input_contract.{field_name}",
                ),
            )
        if self.p_value_interpretation is not None:
            object.__setattr__(
                self,
                "p_value_interpretation",
                _require_non_empty_text(
                    self.p_value_interpretation,
                    field_name=(
                        "method_quantitative_input_contract.p_value_interpretation"
                    ),
                ),
            )
        if not isinstance(cast(object, self.quantitative_input_required), bool):
            raise ValueError(
                "method_quantitative_input_contract.quantitative_input_required "
                "must be a bool"
            )
        if not isinstance(cast(object, self.no_implicit_transformation), bool):
            raise ValueError(
                "method_quantitative_input_contract.no_implicit_transformation "
                "must be a bool"
            )
        if self.quantitative_input_required:
            if not self.accepted_scales:
                raise ValueError(
                    "method_quantitative_input_contract.accepted_scales must not "
                    "be empty when quantitative input is required"
                )
            if not self.accepted_meanings:
                raise ValueError(
                    "method_quantitative_input_contract.accepted_meanings must not "
                    "be empty when quantitative input is required"
                )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-safe method contract payload."""

        return {
            "contract_version": self.contract_version,
            "method_id": self.method_id,
            "accepted_scales": [
                scale.value
                for scale in cast(tuple[IntensityScaleKind, ...], self.accepted_scales)
            ],
            "accepted_meanings": [
                meaning.value
                for meaning in cast(
                    tuple[QuantitativeMeaning, ...], self.accepted_meanings
                )
            ],
            "accepted_activity_profile_axes": [
                axis.value
                for axis in cast(
                    tuple[ActivityProfileAxis, ...],
                    self.accepted_activity_profile_axes,
                )
            ],
            "accepted_activity_quantitative_semantics": [
                quantity.value
                for quantity in cast(
                    tuple[ActivityQuantitativeSemantics, ...],
                    self.accepted_activity_quantitative_semantics,
                )
            ],
            "required_centring": self.required_centring,
            "required_standardisation": self.required_standardisation,
            "missing_value_treatment": self.missing_value_treatment,
            "profile_axis_requirements": self.profile_axis_requirements,
            "statistical_interpretation": self.statistical_interpretation,
            "p_value_interpretation": self.p_value_interpretation,
            "quantitative_input_required": bool(self.quantitative_input_required),
            "scale_sensitivity": self.scale_sensitivity,
            "no_implicit_transformation": bool(self.no_implicit_transformation),
        }


@dataclass(frozen=True, slots=True)
class ResolvedMethodQuantitativeInputContract:
    """Method input contract plus the observed input semantics used in a run."""

    contract: MethodQuantitativeInputContract
    resolved_scale: IntensityScaleKind | str | None = None
    resolved_meaning: QuantitativeMeaning | str | None = None
    resolved_activity_profile_axis: ActivityProfileAxis | str | None = None
    resolved_activity_quantitative_semantics: (
        ActivityQuantitativeSemantics | str | None
    ) = None
    enforcement_context: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.contract, MethodQuantitativeInputContract):
            raise ValueError(
                "resolved_method_quantitative_input_contract.contract must be "
                "MethodQuantitativeInputContract"
            )
        object.__setattr__(
            self,
            "resolved_scale",
            _normalize_optional_intensity_scale(self.resolved_scale),
        )
        object.__setattr__(
            self,
            "resolved_meaning",
            _normalize_optional_quantitative_meaning(self.resolved_meaning),
        )
        object.__setattr__(
            self,
            "resolved_activity_profile_axis",
            _normalize_optional_activity_profile_axis(
                self.resolved_activity_profile_axis
            ),
        )
        object.__setattr__(
            self,
            "resolved_activity_quantitative_semantics",
            _normalize_optional_activity_quantitative_semantics(
                self.resolved_activity_quantitative_semantics
            ),
        )
        object.__setattr__(
            self,
            "enforcement_context",
            (
                ""
                if self.enforcement_context is None
                else str(self.enforcement_context).strip()
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-safe resolved method input contract payload."""

        payload = self.contract.to_payload()
        resolved_scale = cast(IntensityScaleKind | None, self.resolved_scale)
        resolved_meaning = cast(QuantitativeMeaning | None, self.resolved_meaning)
        resolved_axis = cast(
            ActivityProfileAxis | None,
            self.resolved_activity_profile_axis,
        )
        resolved_activity_quantity = cast(
            ActivityQuantitativeSemantics | None,
            self.resolved_activity_quantitative_semantics,
        )
        payload.update(
            {
                "resolved_scale": (
                    None if resolved_scale is None else resolved_scale.value
                ),
                "resolved_meaning": (
                    None if resolved_meaning is None else resolved_meaning.value
                ),
                "resolved_activity_profile_axis": (
                    None if resolved_axis is None else resolved_axis.value
                ),
                "resolved_activity_quantitative_semantics": (
                    None
                    if resolved_activity_quantity is None
                    else resolved_activity_quantity.value
                ),
                "enforcement_context": self.enforcement_context or None,
            }
        )
        return payload


def resolve_activity_input_contract(
    *,
    activity_input: ActivityInputMatrix,
    contract: MethodQuantitativeInputContract,
    context: str,
) -> ResolvedMethodQuantitativeInputContract:
    """Validate typed activity-input semantics against a method contract."""

    if not isinstance(activity_input, ActivityInputMatrix):
        raise WorkflowBoundaryError(f"{context} must be ActivityInputMatrix")
    accepted_axes = cast(
        tuple[ActivityProfileAxis, ...],
        contract.accepted_activity_profile_axes,
    )
    accepted_quantities = cast(
        tuple[ActivityQuantitativeSemantics, ...],
        contract.accepted_activity_quantitative_semantics,
    )
    observed_axis = cast(ActivityProfileAxis, activity_input.semantics.profile_axis)
    observed_quantity = cast(
        ActivityQuantitativeSemantics,
        activity_input.semantics.quantitative_semantics,
    )
    if accepted_axes and observed_axis not in accepted_axes:
        allowed = _format_values(axis.value for axis in accepted_axes)
        raise WorkflowBoundaryError(
            f"{context} method {contract.method_id!r} requires activity profile "
            f"axis in {{{allowed}}}; got {observed_axis.value!r}; no implicit "
            "transformation is applied"
        )
    if accepted_quantities and observed_quantity not in accepted_quantities:
        allowed = _format_values(quantity.value for quantity in accepted_quantities)
        raise WorkflowBoundaryError(
            f"{context} method {contract.method_id!r} requires activity "
            f"quantitative semantics in {{{allowed}}}; got "
            f"{observed_quantity.value!r}; no implicit transformation is applied"
        )
    return ResolvedMethodQuantitativeInputContract(
        contract=contract,
        resolved_activity_profile_axis=observed_axis,
        resolved_activity_quantitative_semantics=observed_quantity,
        enforcement_context=context,
    )


def _normalize_intensity_scales(
    values: Iterable[IntensityScaleKind | str],
) -> tuple[IntensityScaleKind, ...]:
    return tuple(_normalize_intensity_scale(value) for value in values)


def _normalize_quantitative_meanings(
    values: Iterable[QuantitativeMeaning | str],
) -> tuple[QuantitativeMeaning, ...]:
    return tuple(_normalize_quantitative_meaning(value) for value in values)


def _normalize_activity_profile_axes(
    values: Iterable[ActivityProfileAxis | str],
) -> tuple[ActivityProfileAxis, ...]:
    return tuple(_normalize_activity_profile_axis(value) for value in values)


def _normalize_activity_quantitative_semantics_values(
    values: Iterable[ActivityQuantitativeSemantics | str],
) -> tuple[ActivityQuantitativeSemantics, ...]:
    return tuple(_normalize_activity_quantitative_semantics(value) for value in values)


def _normalize_intensity_scale(value: IntensityScaleKind | str) -> IntensityScaleKind:
    if isinstance(value, IntensityScaleKind):
        return value
    return IntensityScaleKind(str(value))


def _normalize_quantitative_meaning(
    value: QuantitativeMeaning | str,
) -> QuantitativeMeaning:
    if isinstance(value, QuantitativeMeaning):
        return value
    return QuantitativeMeaning(str(value))


def _normalize_activity_profile_axis(
    value: ActivityProfileAxis | str,
) -> ActivityProfileAxis:
    if isinstance(value, ActivityProfileAxis):
        return value
    return ActivityProfileAxis(str(value))


def _normalize_activity_quantitative_semantics(
    value: ActivityQuantitativeSemantics | str,
) -> ActivityQuantitativeSemantics:
    if isinstance(value, ActivityQuantitativeSemantics):
        return value
    return ActivityQuantitativeSemantics(str(value))


def _normalize_optional_intensity_scale(
    value: IntensityScaleKind | str | None,
) -> IntensityScaleKind | None:
    if value is None:
        return None
    return _normalize_intensity_scale(value)


def _normalize_optional_quantitative_meaning(
    value: QuantitativeMeaning | str | None,
) -> QuantitativeMeaning | None:
    if value is None:
        return None
    return _normalize_quantitative_meaning(value)


def _normalize_optional_activity_profile_axis(
    value: ActivityProfileAxis | str | None,
) -> ActivityProfileAxis | None:
    if value is None:
        return None
    return _normalize_activity_profile_axis(value)


def _normalize_optional_activity_quantitative_semantics(
    value: ActivityQuantitativeSemantics | str | None,
) -> ActivityQuantitativeSemantics | None:
    if value is None:
        return None
    return _normalize_activity_quantitative_semantics(value)


def _require_non_empty_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{field_name} must be a non-empty string")
    normalized = value.strip()
    if not normalized:
        raise ValueError(f"{field_name} must be a non-empty string")
    return normalized


def _format_values(values: Iterable[str]) -> str:
    return ", ".join(repr(str(value)) for value in sorted(values))


def method_contracts_to_markdown_table(
    contracts: Iterable[MethodQuantitativeInputContract],
) -> str:
    """Render a compact documentation table from method contracts."""

    rows = [
        "| Method | Accepted scale | Accepted meaning | Required centring/standardisation | Missing values | Profile axis | Statistical interpretation | P-value interpretation |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for contract in contracts:
        payload = contract.to_payload()
        centring_standardisation = (
            f"{payload['required_centring']}; {payload['required_standardisation']}"
        )
        rows.append(
            "| "
            + " | ".join(
                _escape_markdown_cell(value)
                for value in (
                    payload["method_id"],
                    ", ".join(cast(list[str], payload["accepted_scales"])) or "n/a",
                    ", ".join(cast(list[str], payload["accepted_meanings"])) or "n/a",
                    centring_standardisation,
                    payload["missing_value_treatment"],
                    payload["profile_axis_requirements"],
                    payload["statistical_interpretation"],
                    payload["p_value_interpretation"] or "none",
                )
            )
            + " |"
        )
    return "\n".join(rows)


def _escape_markdown_cell(value: object) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def contract_payload_from_mapping(
    value: Mapping[str, object],
) -> Mapping[str, object]:
    """Typed identity helper used by tests and provenance readers."""

    return dict(value)


__all__ = [
    "MethodQuantitativeInputContract",
    "ResolvedMethodQuantitativeInputContract",
    "contract_payload_from_mapping",
    "method_contracts_to_markdown_table",
    "resolve_activity_input_contract",
]
