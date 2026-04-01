from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from .constants import (
    DEFAULT_CORRECTED_COLS,
    DEFAULT_PHOSPHO_COLS,
    DEFAULT_TOTAL_COLS,
)
from .validation.compatibility import validate_core_column_alignment


@dataclass(frozen=True, slots=True)
class DatasetSchema:
    """Immutable aligned column groups used by core dataset preprocessing."""

    total_cols: tuple[str, ...] = DEFAULT_TOTAL_COLS
    phospho_cols: tuple[str, ...] = DEFAULT_PHOSPHO_COLS
    corrected_cols: tuple[str, ...] = DEFAULT_CORRECTED_COLS

    def __post_init__(self) -> None:
        total_cols = tuple(self.total_cols)
        phospho_cols = tuple(self.phospho_cols)
        corrected_cols = tuple(self.corrected_cols)

        validate_core_column_alignment(
            total_cols,
            phospho_cols,
            corrected_cols,
            context="Dataset schema",
        )

        object.__setattr__(self, "total_cols", total_cols)
        object.__setattr__(self, "phospho_cols", phospho_cols)
        object.__setattr__(self, "corrected_cols", corrected_cols)

    @property
    def group_to_corrected_col(self) -> Mapping[str, str]:
        return {
            f"group{index}": corrected_col
            for index, corrected_col in enumerate(self.corrected_cols, start=1)
        }

    @classmethod
    def from_groups(
        cls,
        *,
        total_cols: Sequence[str] | None = None,
        phospho_cols: Sequence[str] | None = None,
        corrected_cols: Sequence[str] | None = None,
    ) -> DatasetSchema:
        return cls(
            total_cols=tuple(total_cols or DEFAULT_TOTAL_COLS),
            phospho_cols=tuple(phospho_cols or DEFAULT_PHOSPHO_COLS),
            corrected_cols=tuple(corrected_cols or DEFAULT_CORRECTED_COLS),
        )
