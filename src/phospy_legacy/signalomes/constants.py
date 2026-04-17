from __future__ import annotations

from typing import Final

MODULE_ID_COLUMN: Final[str] = "module_id"
PROTEIN_ID_COLUMN: Final[str] = "protein_id"
SITE_ID_COLUMN: Final[str] = "site_id"
KINASE_COLUMN: Final[str] = "kinase"
SHARE_PERCENT_COLUMN: Final[str] = "share_percent"

TOP_KINASE_CANDIDATES_COLUMN: Final[str] = "top_kinase_candidates"
TOP_KINASE_WEIGHTS_COLUMN: Final[str] = "top_kinase_weights"
TOP_KINASE_TIE_COUNT_COLUMN: Final[str] = "top_kinase_tie_count"
TOP_KINASE_IS_AMBIGUOUS_COLUMN: Final[str] = "top_kinase_is_ambiguous"
TOP_SCORE_COLUMN: Final[str] = "top_score"

DEGREE_COLUMN: Final[str] = "degree"
N_SUBSTRATES_COLUMN: Final[str] = "n_substrates"
MODULE_COUNT_COLUMN: Final[str] = "module_count"
TOTAL_SHARE_PERCENT_COLUMN: Final[str] = "total_share_percent"
IS_KINASE_OF_INTEREST_COLUMN: Final[str] = "is_kinase_of_interest"

SOURCE_KINASE_COLUMN: Final[str] = "source_kinase"
TARGET_KINASE_COLUMN: Final[str] = "target_kinase"
CORRELATION_COLUMN: Final[str] = "correlation"
WEIGHT_COLUMN: Final[str] = "weight"

__all__ = [
    "CORRELATION_COLUMN",
    "DEGREE_COLUMN",
    "IS_KINASE_OF_INTEREST_COLUMN",
    "KINASE_COLUMN",
    "MODULE_COUNT_COLUMN",
    "MODULE_ID_COLUMN",
    "N_SUBSTRATES_COLUMN",
    "PROTEIN_ID_COLUMN",
    "SHARE_PERCENT_COLUMN",
    "SITE_ID_COLUMN",
    "SOURCE_KINASE_COLUMN",
    "TARGET_KINASE_COLUMN",
    "TOP_KINASE_CANDIDATES_COLUMN",
    "TOP_KINASE_IS_AMBIGUOUS_COLUMN",
    "TOP_KINASE_TIE_COUNT_COLUMN",
    "TOP_KINASE_WEIGHTS_COLUMN",
    "TOP_SCORE_COLUMN",
    "TOTAL_SHARE_PERCENT_COLUMN",
    "WEIGHT_COLUMN",
]
