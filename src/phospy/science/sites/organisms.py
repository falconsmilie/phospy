"""Organism parsing for phosphosite identity boundaries."""

from __future__ import annotations

import math
import re
from enum import Enum
from typing import TypeVar

import pandas as pd

from phospy.science.references.models import Organism

ErrorType = TypeVar("ErrorType", bound=Exception)

_ALIAS_SEPARATOR = re.compile(r"[\s_-]+")
_ORGANISM_ALIASES: dict[str, Organism] = {
    "human": Organism.HUMAN,
    "homo sapiens": Organism.HUMAN,
    "h sapiens": Organism.HUMAN,
    "hsapiens": Organism.HUMAN,
    "hsa": Organism.HUMAN,
    "9606": Organism.HUMAN,
    "mouse": Organism.MOUSE,
    "mus musculus": Organism.MOUSE,
    "m musculus": Organism.MOUSE,
    "mmusculus": Organism.MOUSE,
    "mmu": Organism.MOUSE,
    "10090": Organism.MOUSE,
    "rat": Organism.RAT,
    "rattus norvegicus": Organism.RAT,
    "r norvegicus": Organism.RAT,
    "rnorvegicus": Organism.RAT,
    "rno": Organism.RAT,
    "10116": Organism.RAT,
}


def normalize_organism(
    value: object,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> Organism:
    """Resolve one required organism value to the shared ``Organism`` enum."""

    if isinstance(value, Organism):
        return value
    if isinstance(value, Enum):
        value = value.value
    if _is_missing(value):
        raise error_type(f"{field_name} must be a supported organism")
    if not isinstance(value, str):
        raise error_type(f"{field_name} must be a supported organism string")
    token = value.strip()
    if token == "":
        raise error_type(f"{field_name} must be a supported organism")
    alias_key = _normalise_alias_token(token)
    organism = _ORGANISM_ALIASES.get(alias_key)
    if organism is None:
        supported = ", ".join(member.value for member in Organism)
        raise error_type(
            f"{field_name} has unsupported organism {token!r}; "
            f"supported organisms: {supported}"
        )
    return organism


def normalize_optional_organism(
    value: object | None,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> Organism | None:
    """Resolve one optional organism value to ``Organism`` or ``None``."""

    if _is_missing(value):
        return None
    if isinstance(value, str) and value.strip() == "":
        return None
    return normalize_organism(
        value,
        field_name=field_name,
        error_type=error_type,
    )


def _normalise_alias_token(value: str) -> str:
    return _ALIAS_SEPARATOR.sub(" ", value.strip().casefold()).strip()


def _is_missing(value: object | None) -> bool:
    if value is None or value is pd.NA or value is pd.NaT:
        return True
    if isinstance(value, float):
        return math.isnan(value)
    return False


__all__ = [
    "normalize_optional_organism",
    "normalize_organism",
]
