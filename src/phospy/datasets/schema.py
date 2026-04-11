from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from ..errors import InputCompatibilityError
from ..internal.constants import (
    DEFAULT_CORRECTED_COLS,
    DEFAULT_PHOSPHO_COLS,
    DEFAULT_TOTAL_COLS,
)
from ..validation.compatibility import validate_core_column_alignment


@dataclass(frozen=True, slots=True)
class DatasetSchema:
    """Immutable aligned sample/value column groups used by dataset processing.

    DatasetSchema intentionally models only the aligned numeric/sample columns.
    Core structural identifier columns such as gene_names, gene_p_site, and
    localization_prob are fixed package constants rather than user-configurable
    schema fields.
    """

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
    def comparison_groups(self) -> tuple[str, ...]:
        return self.total_cols

    @property
    def group_to_corrected_col(self) -> Mapping[str, str]:
        return dict(zip(self.comparison_groups, self.corrected_cols, strict=True))

    def validate_comparisons(
        self,
        comparisons: Sequence[tuple[str, str]] | None,
        *,
        context: str = "Dataset schema",
    ) -> tuple[tuple[str, str], ...] | None:
        if comparisons is None:
            return None

        resolved = tuple(comparisons)
        valid_groups = frozenset(self.comparison_groups)
        seen: set[tuple[str, str]] = set()
        for left_group, right_group in resolved:
            if left_group not in valid_groups:
                msg = f"{context} contains Unknown comparison group: {left_group}"
                raise InputCompatibilityError(msg)
            if right_group not in valid_groups:
                msg = f"{context} contains Unknown comparison group: {right_group}"
                raise InputCompatibilityError(msg)
            pair = (left_group, right_group)
            if pair in seen:
                msg = (
                    f"{context} contains Duplicate comparison pair: "
                    f"{left_group!r}, {right_group!r}"
                )
                raise InputCompatibilityError(msg)
            seen.add(pair)
        return resolved

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
