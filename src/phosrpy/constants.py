from __future__ import annotations

from typing import TypeAlias

ComparisonSpec: TypeAlias = tuple[str, str]

DEFAULT_TOTAL_COLS = [f"group{i}" for i in range(1, 7)]
DEFAULT_PHOSPHO_COLS = [f"p_group{i}" for i in range(1, 7)]
DEFAULT_CORRECTED_COLS = [f"phospho_corrected_{i}" for i in range(1, 7)]
