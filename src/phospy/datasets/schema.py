from __future__ import annotations

from collections.abc import Mapping, Sequence

from ..internal.constants import (
    DEFAULT_CORRECTED_COLS,
    DEFAULT_PHOSPHO_COLS,
    DEFAULT_TOTAL_COLS,
)
from ..validation.compatibility import validate_core_column_alignment
from ..validation.domain.comparisons import validate_comparison_specs


class DatasetSchema:
    """Immutable aligned sample/value column groups used by dataset processing.

    DatasetSchema intentionally models only the aligned numeric/sample columns.
    Core structural identifier columns such as gene_names, gene_p_site, and
    localization_prob are fixed package constants rather than user-configurable
    schema fields.
    """

    __slots__ = ("_total_cols", "_phospho_cols", "_corrected_cols", "_frozen")

    def __init__(
        self,
        total_cols: Sequence[str] = DEFAULT_TOTAL_COLS,
        phospho_cols: Sequence[str] = DEFAULT_PHOSPHO_COLS,
        corrected_cols: Sequence[str] = DEFAULT_CORRECTED_COLS,
    ) -> None:
        normalized_total_cols = tuple(total_cols)
        normalized_phospho_cols = tuple(phospho_cols)
        normalized_corrected_cols = tuple(corrected_cols)

        validate_core_column_alignment(
            normalized_total_cols,
            normalized_phospho_cols,
            normalized_corrected_cols,
            context="Dataset schema",
        )

        self._frozen = False
        self._total_cols = normalized_total_cols
        self._phospho_cols = normalized_phospho_cols
        self._corrected_cols = normalized_corrected_cols
        self._frozen = True

    def __setattr__(self, name: str, value: object) -> None:
        if getattr(self, "_frozen", False):
            msg = "DatasetSchema is immutable"
            raise AttributeError(msg)
        super().__setattr__(name, value)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, DatasetSchema):
            return NotImplemented
        return (
            self.total_cols,
            self.phospho_cols,
            self.corrected_cols,
        ) == (
            other.total_cols,
            other.phospho_cols,
            other.corrected_cols,
        )

    def __hash__(self) -> int:
        return hash((self.total_cols, self.phospho_cols, self.corrected_cols))

    def __repr__(self) -> str:
        return (
            "DatasetSchema("
            f"total_cols={self.total_cols!r}, "
            f"phospho_cols={self.phospho_cols!r}, "
            f"corrected_cols={self.corrected_cols!r})"
        )

    @property
    def total_cols(self) -> tuple[str, ...]:
        return self._total_cols

    @property
    def phospho_cols(self) -> tuple[str, ...]:
        return self._phospho_cols

    @property
    def corrected_cols(self) -> tuple[str, ...]:
        return self._corrected_cols

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
        return validate_comparison_specs(
            comparison_groups=self.comparison_groups,
            comparisons=comparisons,
            context=context,
        )

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
