from __future__ import annotations

from typing import TypeAlias

ComparisonSpec: TypeAlias = tuple[str, str]

DEFAULT_TOTAL_COLS: tuple[str, ...] = tuple(f"group{i}" for i in range(1, 7))
DEFAULT_PHOSPHO_COLS: tuple[str, ...] = tuple(f"p_group{i}" for i in range(1, 7))
DEFAULT_CORRECTED_COLS: tuple[str, ...] = tuple(
    f"phospho_corrected_{i}" for i in range(1, 7)
)
DEFAULT_TOTAL_SENTINEL: float = 10.0
DEFAULT_PHOSPHO_SENTINEL: float = 12.0
