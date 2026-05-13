"""Reference domain package."""

from __future__ import annotations

from typing import TYPE_CHECKING

__all__ = [
    "Organism",
    "ReferenceBundle",
    "ReferenceManifest",
    "ReferencePreset",
    "SequenceWindowDefinition",
]

if TYPE_CHECKING:
    from phospy.science.references.models import (
        Organism,
        ReferenceBundle,
        ReferenceManifest,
        ReferencePreset,
        SequenceWindowDefinition,
    )


def __getattr__(name: str) -> object:
    if name in __all__:
        from phospy.science.references import models as _models

        return getattr(_models, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
