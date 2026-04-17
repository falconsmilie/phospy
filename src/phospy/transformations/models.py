"""Transformation domain models."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class TransformationKind(str, Enum):
    """Supported transformation-state kinds."""

    RAW = "raw"
    LOG2 = "log2"
    NORMALIZED = "normalized"


@dataclass(frozen=True, slots=True)
class TransformationState:
    """Validated transformation metadata for dataset inputs."""

    kind: TransformationKind = TransformationKind.RAW
    label: str = "raw"
