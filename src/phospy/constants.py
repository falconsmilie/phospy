from __future__ import annotations

from typing import Final, TypeAlias

ComparisonSpec: TypeAlias = tuple[str, str]

# Fixed structural dataset columns. These are intentionally canonical internal
# constants rather than user-configurable schema fields. DatasetSchema only
# models aligned sample/value column groups.
TOTAL_GENE_COLUMN: Final[str] = "genes"
PHOSPHO_UID_COLUMN: Final[str] = "uid"
PHOSPHO_GENE_COLUMN: Final[str] = "gene_names"
GENE_P_SITE_COLUMN: Final[str] = "gene_p_site"
LOCALIZATION_PROB_COLUMN: Final[str] = "localization_prob"
CENTRALIZED_SEQUENCE_COLUMN: Final[str] = "centralized_sequence"
SITE_MATRIX_GENE_COLUMN: Final[str] = "gene"
SITE_MATRIX_P_SITE_COLUMN: Final[str] = "p_site"
SITE_MATRIX_ID_COLUMN: Final[str] = "site_id"

PHOSPHO_REQUIRED_METADATA_COLUMNS: Final[tuple[str, ...]] = (
    PHOSPHO_UID_COLUMN,
    PHOSPHO_GENE_COLUMN,
    GENE_P_SITE_COLUMN,
    LOCALIZATION_PROB_COLUMN,
    CENTRALIZED_SEQUENCE_COLUMN,
)

DEFAULT_TOTAL_COLS: tuple[str, ...] = tuple(f"group{i}" for i in range(1, 7))
DEFAULT_PHOSPHO_COLS: tuple[str, ...] = tuple(f"p_group{i}" for i in range(1, 7))
DEFAULT_CORRECTED_COLS: tuple[str, ...] = tuple(
    f"phospho_corrected_{i}" for i in range(1, 7)
)
DEFAULT_TOTAL_SENTINEL: float = 10.0
DEFAULT_PHOSPHO_SENTINEL: float = 12.0
