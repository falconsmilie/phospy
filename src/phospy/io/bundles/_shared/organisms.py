"""Organism parsing helpers for bundle metadata payloads."""

from __future__ import annotations

from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.primitives import require_str
from phospy.science.references.models import Organism


def parse_optional_organism(value: object, *, field_name: str) -> Organism | None:
    """Parse an optional organism enum token."""

    if value is None:
        return None
    return parse_required_organism(value, field_name=field_name)


def parse_required_organism(value: object, *, field_name: str) -> Organism:
    """Parse a required organism enum token."""

    token = require_str(value, field_name=field_name)
    try:
        return Organism(token)
    except ValueError as exc:
        supported = ", ".join(member.value for member in Organism)
        raise PhosPyInputError(
            f"unsupported organism '{token}' in {field_name}; supported: {supported}"
        ) from exc
