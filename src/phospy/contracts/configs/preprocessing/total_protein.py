"""Total-protein-correction preprocessing policy configuration."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import pandas as pd

from phospy.validation.configs.preprocessing import (
    validate_protein_aware_preparation_config,
    validate_total_protein_correction_config,
    validate_total_protein_correction_identity_config,
)

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
DATASET_PROTEIN_AWARE_PREPARATION_POLICY_DISABLED = "disabled"
DATASET_PROTEIN_AWARE_PREPARATION_POLICY_PREPARE_MODEL_INPUTS = "prepare_model_inputs"
DatasetProteinAwarePreparationPolicy = Literal[
    "disabled",
    "prepare_model_inputs",
]
DATASET_PROTEIN_AWARE_PREPARATION_POLICIES = frozenset(
    {
        DATASET_PROTEIN_AWARE_PREPARATION_POLICY_DISABLED,
        DATASET_PROTEIN_AWARE_PREPARATION_POLICY_PREPARE_MODEL_INPUTS,
    }
)
DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS = (
    "require_unambiguous"
)
DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_ALLOW_MISSING_WITH_REPORT = (
    "allow_missing_with_report"
)
DatasetProteinAwarePreparationMappingPolicy = Literal[
    "require_unambiguous",
    "allow_missing_with_report",
]
DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICIES = frozenset(
    {
        DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS,
        DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_ALLOW_MISSING_WITH_REPORT,
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
        validate_total_protein_correction_identity_config(
            mode=self.mode,
            phosphosite_key=self.phosphosite_key,
            total_protein_key=self.total_protein_key,
            mapping_table=self.mapping_table,
            mapping_phosphosite_key=self.mapping_phosphosite_key,
            mapping_total_protein_key=self.mapping_total_protein_key,
            matching_policy=self.matching_policy,
            duplicate_policy=self.duplicate_policy,
            unmatched_policy=self.unmatched_policy,
            supported_modes=DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODES,
            supported_matching_policies=(
                DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICIES
            ),
            supported_duplicate_policies=(
                DATASET_TOTAL_PROTEIN_CORRECTION_DUPLICATE_POLICIES
            ),
            supported_unmatched_policies=(
                DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICIES
            ),
            mode_direct=DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_DIRECT,
            mode_mapping_table=(
                DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_MAPPING_TABLE
            ),
            matching_policy_gene_symbol_normalised=(
                DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MATCHING_POLICY_GENE_SYMBOL_NORMALISED
            ),
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
        validate_total_protein_correction_config(
            policy=self.policy,
            identity=self.identity,
            supported_policies=DATASET_TOTAL_PROTEIN_CORRECTION_POLICIES,
            identity_type=DatasetTotalProteinCorrectionIdentityConfig,
        )


@dataclass(frozen=True, slots=True)
class DatasetProteinAwarePreparationConfig:
    """Config-only intent for future protein-aware model-input preparation.

    This contract is separate from `DatasetTotalProteinCorrectionConfig`: it
    does not subtract total protein, align matrices, decide site eligibility, or
    run full joint PTM/protein modelling. `policy="prepare_model_inputs"` only
    declares that a future preparation lane should build aligned
    phosphosite/protein inputs.
    """

    policy: DatasetProteinAwarePreparationPolicy = (
        DATASET_PROTEIN_AWARE_PREPARATION_POLICY_DISABLED
    )
    protein_mapping_policy: DatasetProteinAwarePreparationMappingPolicy = (
        DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS
    )

    def __post_init__(self) -> None:
        validate_protein_aware_preparation_config(
            policy=self.policy,
            protein_mapping_policy=self.protein_mapping_policy,
            supported_policies=DATASET_PROTEIN_AWARE_PREPARATION_POLICIES,
            supported_mapping_policies=(
                DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICIES
            ),
        )


__all__ = [
    "DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICIES",
    "DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_ALLOW_MISSING_WITH_REPORT",
    "DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS",
    "DATASET_PROTEIN_AWARE_PREPARATION_POLICIES",
    "DATASET_PROTEIN_AWARE_PREPARATION_POLICY_DISABLED",
    "DATASET_PROTEIN_AWARE_PREPARATION_POLICY_PREPARE_MODEL_INPUTS",
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
    "DatasetProteinAwarePreparationConfig",
    "DatasetProteinAwarePreparationMappingPolicy",
    "DatasetProteinAwarePreparationPolicy",
    "DatasetTotalProteinCorrectionConfig",
    "DatasetTotalProteinCorrectionDuplicatePolicy",
    "DatasetTotalProteinCorrectionIdentityConfig",
    "DatasetTotalProteinCorrectionIdentityMatchingPolicy",
    "DatasetTotalProteinCorrectionIdentityMode",
    "DatasetTotalProteinCorrectionPolicy",
    "DatasetTotalProteinCorrectionUnmatchedPolicy",
]
