"""Narrow parsing helpers for CLI/runtime adapter boundaries."""

from __future__ import annotations

from phospy.errors.input import PhosPyInputError
from phospy.references.models import Organism, ReferencePreset


def organism_from_value(value: Organism | str | None) -> Organism | None:
    """Parse an optional organism token into ``Organism``."""

    if value is None:
        return None
    if isinstance(value, Organism):
        return value
    if not isinstance(value, str):
        raise PhosPyInputError(
            "unsupported organism value type. expected Organism, str, or None"
        )
    normalized = value.strip().lower()
    for organism in Organism:
        if organism.value == normalized:
            return organism
    supported = ", ".join(member.value for member in Organism)
    raise PhosPyInputError(
        f"unsupported organism '{value}'. supported organisms: {supported}"
    )


def reference_preset_from_value(value: ReferencePreset | str) -> ReferencePreset:
    """Parse a reference preset token into ``ReferencePreset``."""

    if isinstance(value, ReferencePreset):
        return value
    if not isinstance(value, str):
        raise PhosPyInputError(
            "unsupported reference preset value type. expected ReferencePreset or str"
        )
    normalized = value.strip().lower()
    for preset in ReferencePreset:
        if preset.value == normalized:
            return preset
    supported = ", ".join(member.value for member in ReferencePreset)
    raise PhosPyInputError(
        f"unsupported reference preset '{value}'. supported presets: {supported}"
    )
