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

# Persisted core preprocessing output basenames.
CORE_TOTAL_UNIQUE_BASENAME: Final[str] = "df_total_unique"
CORE_TOTAL_FILTERED_BASENAME: Final[str] = "df_total_filtered"
CORE_PHOSPHO_FILTERED_BASENAME: Final[str] = "df_phospho_filtered"
CORE_PHOSPHO_CORRECTED_BASENAME: Final[str] = "df_phospho_corrected"
CORE_PHOSR_INPUT_BASENAME: Final[str] = "phosr_input"
CORE_SITE_MATRIX_BASENAME: Final[str] = "mat_phospho_corrected"
CORE_SITE_SEQUENCES_BASENAME: Final[str] = "site_sequences"
CORE_OUTPUT_ARTIFACT_BASENAMES: Final[tuple[str, ...]] = (
    CORE_TOTAL_UNIQUE_BASENAME,
    CORE_TOTAL_FILTERED_BASENAME,
    CORE_PHOSPHO_FILTERED_BASENAME,
    CORE_PHOSPHO_CORRECTED_BASENAME,
    CORE_PHOSR_INPUT_BASENAME,
    CORE_SITE_MATRIX_BASENAME,
    CORE_SITE_SEQUENCES_BASENAME,
)

# Persisted kinase activity output basenames and filenames.
KINASE_ACTIVITY_MATRIX_BASENAME: Final[str] = "kinase_activity_matrix"
KSEA_SCORES_BASENAME: Final[str] = "ksea_scores"
KSEA_COUNTS_BASENAME: Final[str] = "ksea_counts"
KINASE_TARGET_COUNTS_BASENAME: Final[str] = "kinase_target_counts"
KINASE_TARGET_TABLE_BASENAME: Final[str] = "kinase_target_table"
KINASE_OUTPUT_ARTIFACT_BASENAMES: Final[tuple[str, ...]] = (
    KINASE_ACTIVITY_MATRIX_BASENAME,
    KSEA_SCORES_BASENAME,
    KSEA_COUNTS_BASENAME,
    KINASE_TARGET_COUNTS_BASENAME,
    KINASE_TARGET_TABLE_BASENAME,
)
KINASE_ACTIVITY_MATRIX_FILENAME: Final[str] = f"{KINASE_ACTIVITY_MATRIX_BASENAME}.csv"
KSEA_SCORES_FILENAME: Final[str] = f"{KSEA_SCORES_BASENAME}.csv"
KSEA_COUNTS_FILENAME: Final[str] = f"{KSEA_COUNTS_BASENAME}.csv"
KINASE_TARGET_COUNTS_FILENAME: Final[str] = f"{KINASE_TARGET_COUNTS_BASENAME}.csv"
KINASE_TARGET_TABLE_FILENAME: Final[str] = f"{KINASE_TARGET_TABLE_BASENAME}.csv"
KINASE_OUTPUT_FILENAMES: Final[tuple[str, ...]] = (
    KINASE_ACTIVITY_MATRIX_FILENAME,
    KSEA_SCORES_FILENAME,
    KSEA_COUNTS_FILENAME,
    KINASE_TARGET_COUNTS_FILENAME,
    KINASE_TARGET_TABLE_FILENAME,
)

RUN_MANIFEST_FILENAME: Final[str] = "run_manifest.json"

DEFAULT_TOTAL_COLS: tuple[str, ...] = tuple(f"group{i}" for i in range(1, 7))
DEFAULT_PHOSPHO_COLS: tuple[str, ...] = tuple(f"p_group{i}" for i in range(1, 7))
DEFAULT_CORRECTED_COLS: tuple[str, ...] = tuple(
    f"phospho_corrected_{i}" for i in range(1, 7)
)
DEFAULT_TOTAL_SENTINEL: float = 10.0
DEFAULT_PHOSPHO_SENTINEL: float = 12.0


# Bundled kinase reference provider constants.
BUNDLED_REFERENCE_PROVIDER_NAME: Final[str] = "BundledReferenceProvider"
BUNDLED_REFERENCE_SOURCE: Final[str] = "phospy-bundled"
BUNDLED_REFERENCE_VERSION: Final[str] = "2026.04"
BUNDLED_REFERENCE_AUTO: Final[str] = "auto"
BUNDLED_REFERENCE_SPECIES_ALIASES: Final[dict[str, str]] = {
    "rat": "rat",
    "rattus_norvegicus": "rat",
    "rattus norvegicus": "rat",
}
BUNDLED_REFERENCE_DEFAULTS: Final[dict[str, str]] = {
    "rat": "l6_native",
}
BUNDLED_REFERENCE_ALIASES: Final[dict[str, dict[str, str]]] = {
    "rat": {
        "auto": "l6_native",
        "l6": "l6_native",
        "l6_native": "l6_native",
    },
}
