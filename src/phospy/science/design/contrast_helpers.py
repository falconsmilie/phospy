"""Helpers for constructing explicit differential contrasts."""

from __future__ import annotations

from collections import Counter
from typing import cast

from phospy.errors.validation import WorkflowValidationError
from phospy.science.design.models import Contrast, ExperimentalDesign


def all_pairwise_contrasts(design: ExperimentalDesign) -> tuple[Contrast, ...]:
    """Return one explicit contrast for every condition pair.

    Conditions follow their first appearance in ``design``. For each pair, the
    later condition is the numerator and the earlier condition is the
    denominator, so names follow ``<numerator>_vs_<denominator>``.
    """

    condition_labels = _condition_labels(design)
    contrasts = tuple(
        Contrast(
            name=_contrast_name(
                numerator_condition=numerator_condition,
                denominator_condition=denominator_condition,
            ),
            numerator_condition=numerator_condition,
            denominator_condition=denominator_condition,
        )
        for numerator_index, numerator_condition in enumerate(condition_labels)
        for denominator_condition in condition_labels[:numerator_index]
    )
    _validate_unique_contrast_names(contrasts)
    return contrasts


def contrasts_vs_control(
    design: ExperimentalDesign,
    control_condition: str,
) -> tuple[Contrast, ...]:
    """Return each non-control condition as numerator versus control."""

    condition_labels = _condition_labels(design)
    control = _normalize_control_condition(control_condition)
    if control not in condition_labels:
        known = ", ".join(repr(condition) for condition in condition_labels)
        raise WorkflowValidationError(
            f"control condition {control!r} is not present in experimental "
            f"design conditions: {known}"
        )
    contrasts = tuple(
        Contrast(
            name=_contrast_name(
                numerator_condition=condition,
                denominator_condition=control,
            ),
            numerator_condition=condition,
            denominator_condition=control,
        )
        for condition in condition_labels
        if condition != control
    )
    _validate_unique_contrast_names(contrasts)
    return contrasts


def _condition_labels(design: ExperimentalDesign) -> tuple[str, ...]:
    if not isinstance(cast(object, design), ExperimentalDesign):
        raise WorkflowValidationError("contrast helpers require an ExperimentalDesign")
    return design.condition_labels()


def _normalize_control_condition(control_condition: str) -> str:
    if not isinstance(control_condition, str):
        raise WorkflowValidationError("control_condition must be a string")
    normalized = control_condition.strip()
    if normalized == "":
        raise WorkflowValidationError("control_condition must be a non-empty string")
    return normalized


def _contrast_name(
    *,
    numerator_condition: str,
    denominator_condition: str,
) -> str:
    return f"{numerator_condition}_vs_{denominator_condition}"


def _validate_unique_contrast_names(contrasts: tuple[Contrast, ...]) -> None:
    duplicate_names = sorted(
        name
        for name, count in Counter(contrast.name for contrast in contrasts).items()
        if count > 1
    )
    if duplicate_names:
        raise WorkflowValidationError(
            "contrast helper generated duplicate contrast names: "
            + ", ".join(duplicate_names)
            + ". Rename conditions or supply explicit Contrast objects."
        )


__all__ = [
    "all_pairwise_contrasts",
    "contrasts_vs_control",
]
