from __future__ import annotations

from typing import TypeAlias

ComparisonSpec: TypeAlias = tuple[str, str]

DEFAULT_TOTAL_COLS = [f"group{i}" for i in range(1, 7)]
DEFAULT_PHOSPHO_COLS = [f"p_group{i}" for i in range(1, 7)]
DEFAULT_CORRECTED_COLS = [f"phospho_corrected_{i}" for i in range(1, 7)]

# Core defaults stay structural rather than embedding any study-specific labels.
DEFAULT_COMPARISONS: list[ComparisonSpec] = [
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

