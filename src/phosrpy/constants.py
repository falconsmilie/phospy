from __future__ import annotations

from typing import TypeAlias

ComparisonPair: TypeAlias = tuple[str, str]
LegacyComparisonSpec: TypeAlias = tuple[str, str, str, str]
ComparisonSpec: TypeAlias = ComparisonPair | LegacyComparisonSpec

DEFAULT_TOTAL_COLS = [f"group{i}" for i in range(1, 7)]
DEFAULT_PHOSPHO_COLS = [f"p_group{i}" for i in range(1, 7)]
DEFAULT_CORRECTED_COLS = [f"phospho_corrected_{i}" for i in range(1, 7)]

# Core defaults stay structural rather than embedding any study-specific labels.
DEFAULT_COMPARISON_PAIRS: list[ComparisonPair] = [
    ("group1", "group4"),
    ("group2", "group5"),
    ("group3", "group6"),
    ("group1", "group2"),
    ("group1", "group3"),
    ("group2", "group3"),
    ("group4", "group5"),
    ("group4", "group6"),
    ("group5", "group6"),
]

# Backward-compatible alias for older internal code paths or callers.
DEFAULT_COMPARISONS = DEFAULT_COMPARISON_PAIRS
