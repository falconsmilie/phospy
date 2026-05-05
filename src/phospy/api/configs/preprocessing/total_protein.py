"""Total-protein-correction preprocessing policy configuration."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import pandas as pd

from phospy.errors.input import PhosPyInputError

DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE = "none"
DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL = "subtract_log_total"
DatasetTotalProteinCorrectionPolicy = Literal[
    "none",
    "subtract_log_total",
]
DATASET_TOTAL_PROTEIN_CORRECTION_POLICIES = frozenset(
    {
        DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
        DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL,
    }
)
DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_DIRECT = "direct"
DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_MAPPING_TABLE = "mapping_table"
DatasetTotalProteinCorrectionIdentityMode = Literal["direct", "mapping_table"]
DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODES = frozenset(
    {
        DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_DIRECT,
        DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_MAPPING_TABLE,
    }
)
DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICY_STRICT = "strict"
DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICY_GENE_SYMBOL_NORMALISED = (
    "gene_symbol_normalised"
)
DatasetTotalProteinCorrectionIdentityMatchingPolicy = Literal[
    "strict",
    "gene_symbol_normalised",
]
DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICIES = frozenset(
    {
        DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICY_STRICT,
        DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICY_GENE_SYMBOL_NORMALISED,
    }
)
DATASET_TOTAL_PROTEIN_CORRECTION_DUPLICATE_POLICY_ERROR = "error"
DatasetTotalProteinCorrectionDuplicatePolicy = Literal["error"]
DATASET_TOTAL_PROTEIN_CORRECTION_DUPLICATE_POLICIES = frozenset(
    {DATASET_TOTAL_PROTEIN_CORRECTION_DUPLICATE_POLICY_ERROR}
)
DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ERROR = "error"
DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ALLOW_UNCORRECTED = (
    "allow_uncorrected"
)
DatasetTotalProteinCorrectionUnmatchedPolicy = Literal["error", "allow_uncorrected"]
DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICIES = frozenset(
    {
        DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ERROR,
        DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ALLOW_UNCORRECTED,
    }
)


@dataclass(frozen=True, slots=True)
class DatasetTotalProteinCorrectionIdentityConfig:
    """Identity-mapping policy for total/protein correction matching.

    Supported modes:

    - `"direct"`: map `site_metadata[phosphosite_key]` directly to a key in the
      total-protein table (currently resolved from `total.index`).
    - `"mapping_table"`: map phosphosite keys to total-protein keys through an
      explicit two-column mapping table.

    Supported matching policies:

    - `"strict"`: compare identity keys exactly (after trimming surrounding
      whitespace only).
    - `"gene_symbol_normalised"`: compare gene-symbol identity keys after
      uppercasing, which is biologically lossy and must be explicitly chosen.
    """

    mode: DatasetTotalProteinCorrectionIdentityMode = (
        DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_DIRECT
    )
    phosphosite_key: str = "gene_symbol"
    total_protein_key: str = "__index__"
    mapping_table: pd.DataFrame | None = None
    mapping_phosphosite_key: str | None = None
    mapping_total_protein_key: str | None = None
    matching_policy: DatasetTotalProteinCorrectionIdentityMatchingPolicy = (
        DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICY_STRICT
    )
    duplicate_policy: DatasetTotalProteinCorrectionDuplicatePolicy = (
        DATASET_TOTAL_PROTEIN_CORRECTION_DUPLICATE_POLICY_ERROR
    )
    unmatched_policy: DatasetTotalProteinCorrectionUnmatchedPolicy = (
        DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ERROR
    )

    def __post_init__(self) -> None:
        mode = self.mode
        if mode not in DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODES:
            supported = ", ".join(
                sorted(DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODES)
            )
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                f"identity.mode must be one of: {supported}"
            )

        phosphosite_key = self.phosphosite_key
        if not isinstance(phosphosite_key, str) or not phosphosite_key.strip():
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                "identity.phosphosite_key must be a non-empty string"
            )
        total_protein_key = self.total_protein_key
        if not isinstance(total_protein_key, str) or not total_protein_key.strip():
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                "identity.total_protein_key must be a non-empty string"
            )
        matching_policy = self.matching_policy
        if (
            matching_policy
            not in DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICIES
        ):
            supported = ", ".join(
                sorted(DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICIES)
            )
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                f"identity.matching_policy must be one of: {supported}"
            )

        duplicate_policy = self.duplicate_policy
        if duplicate_policy not in DATASET_TOTAL_PROTEIN_CORRECTION_DUPLICATE_POLICIES:
            supported = ", ".join(
                sorted(DATASET_TOTAL_PROTEIN_CORRECTION_DUPLICATE_POLICIES)
            )
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                f"identity.duplicate_policy must be one of: {supported}"
            )
        unmatched_policy = self.unmatched_policy
        if unmatched_policy not in DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICIES:
            supported = ", ".join(
                sorted(DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICIES)
            )
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                f"identity.unmatched_policy must be one of: {supported}"
            )

        mapping_table = self.mapping_table
        mapping_phosphosite_key = self.mapping_phosphosite_key
        mapping_total_protein_key = self.mapping_total_protein_key
        uses_gene_symbol_keys = "gene_symbol" in {
            str(phosphosite_key).strip().lower(),
            str(total_protein_key).strip().lower(),
            (
                ""
                if mapping_phosphosite_key is None
                else str(mapping_phosphosite_key).strip().lower()
            ),
            (
                ""
                if mapping_total_protein_key is None
                else str(mapping_total_protein_key).strip().lower()
            ),
        }
        if (
            matching_policy
            == DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICY_GENE_SYMBOL_NORMALISED
            and not uses_gene_symbol_keys
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                "identity.matching_policy='gene_symbol_normalised' requires at "
                "least one gene_symbol identity key "
                "(phosphosite_key/total_protein_key/mapping keys)"
            )

        if mode == DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_DIRECT:
            if mapping_table is not None:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.total_protein_correction."
                    "identity.mapping_table must be None when identity.mode='direct'"
                )
            if mapping_phosphosite_key is not None:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.total_protein_correction."
                    "identity.mapping_phosphosite_key must be None when "
                    "identity.mode='direct'"
                )
            if mapping_total_protein_key is not None:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.total_protein_correction."
                    "identity.mapping_total_protein_key must be None when "
                    "identity.mode='direct'"
                )
            return

        if mode != DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_MAPPING_TABLE:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                "identity contains an unsupported mode"
            )
        if mapping_table is None:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                "identity.mapping_table is required when identity.mode='mapping_table'"
            )
        if not isinstance(mapping_table, pd.DataFrame):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                "identity.mapping_table must be a pandas DataFrame"
            )
        if (
            not isinstance(mapping_phosphosite_key, str)
            or not mapping_phosphosite_key.strip()
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                "identity.mapping_phosphosite_key must be a non-empty string when "
                "identity.mode='mapping_table'"
            )
        if (
            not isinstance(mapping_total_protein_key, str)
            or not mapping_total_protein_key.strip()
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                "identity.mapping_total_protein_key must be a non-empty string when "
                "identity.mode='mapping_table'"
            )


@dataclass(frozen=True, slots=True)
class DatasetTotalProteinCorrectionConfig:
    """Public total/protein correction policy options for dataset building.

    - `"none"`: do not apply total/protein correction.
    - `"subtract_log_total"`: subtract matched log-scale total-protein abundance
      from log-scale phosphosite abundance in the builder preprocessing lane:
      `log2_phospho - log2_total`.
    """

    policy: DatasetTotalProteinCorrectionPolicy = (
        DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE
    )
    identity: DatasetTotalProteinCorrectionIdentityConfig = field(
        default_factory=DatasetTotalProteinCorrectionIdentityConfig
    )

    def __post_init__(self) -> None:
        policy = self.policy
        if policy not in DATASET_TOTAL_PROTEIN_CORRECTION_POLICIES:
            supported = ", ".join(sorted(DATASET_TOTAL_PROTEIN_CORRECTION_POLICIES))
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                f"policy must be one of: {supported}"
            )
        if not isinstance(self.identity, DatasetTotalProteinCorrectionIdentityConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                "identity must be a DatasetTotalProteinCorrectionIdentityConfig"
            )


__all__ = [
    "DATASET_TOTAL_PROTEIN_CORRECTION_DUPLICATE_POLICIES",
    "DATASET_TOTAL_PROTEIN_CORRECTION_DUPLICATE_POLICY_ERROR",
    "DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICIES",
    "DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICY_GENE_SYMBOL_NORMALISED",
    "DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICY_STRICT",
    "DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_DIRECT",
    "DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_MAPPING_TABLE",
    "DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODES",
    "DATASET_TOTAL_PROTEIN_CORRECTION_POLICIES",
    "DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE",
    "DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL",
    "DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICIES",
    "DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ALLOW_UNCORRECTED",
    "DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ERROR",
    "DatasetTotalProteinCorrectionConfig",
    "DatasetTotalProteinCorrectionDuplicatePolicy",
    "DatasetTotalProteinCorrectionIdentityConfig",
    "DatasetTotalProteinCorrectionIdentityMatchingPolicy",
    "DatasetTotalProteinCorrectionIdentityMode",
    "DatasetTotalProteinCorrectionPolicy",
    "DatasetTotalProteinCorrectionUnmatchedPolicy",
]
