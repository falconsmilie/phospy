"""Site-sequence and readiness preprocessing policy configuration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from phospy.errors.input import PhosPyInputError
from phospy.validation.configs.preprocessing import (
    validate_site_sequence_resolution_config,
)

DATASET_SITE_SEQUENCE_RESOLUTION_MODE_VALIDATE_EXISTING_AND_FILL_MISSING = (
    "validate_existing_and_fill_missing"
)
DATASET_SITE_SEQUENCE_RESOLUTION_MODE_FILL_MISSING_ONLY = "fill_missing_only"
DATASET_SITE_SEQUENCE_RESOLUTION_MODE_VALIDATE_EXISTING_ONLY = "validate_existing_only"
DATASET_SITE_SEQUENCE_RESOLUTION_MODE_REPLACE_EXISTING = "replace_existing"
DatasetSiteSequenceResolutionMode = Literal[
    "validate_existing_and_fill_missing",
    "fill_missing_only",
    "validate_existing_only",
    "replace_existing",
]
DATASET_SITE_SEQUENCE_RESOLUTION_MODES = frozenset(
    {
        DATASET_SITE_SEQUENCE_RESOLUTION_MODE_VALIDATE_EXISTING_AND_FILL_MISSING,
        DATASET_SITE_SEQUENCE_RESOLUTION_MODE_FILL_MISSING_ONLY,
        DATASET_SITE_SEQUENCE_RESOLUTION_MODE_VALIDATE_EXISTING_ONLY,
        DATASET_SITE_SEQUENCE_RESOLUTION_MODE_REPLACE_EXISTING,
    }
)
DATASET_SITE_SEQUENCE_CONFLICT_POLICY_ERROR = "error"
DATASET_SITE_SEQUENCE_CONFLICT_POLICY_PRESERVE_EXISTING = "preserve_existing"
DATASET_SITE_SEQUENCE_CONFLICT_POLICY_REPLACE_EXISTING = "replace_existing"
DatasetSiteSequenceConflictPolicy = Literal[
    "error",
    "preserve_existing",
    "replace_existing",
]
DATASET_SITE_SEQUENCE_CONFLICT_POLICIES = frozenset(
    {
        DATASET_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
        DATASET_SITE_SEQUENCE_CONFLICT_POLICY_PRESERVE_EXISTING,
        DATASET_SITE_SEQUENCE_CONFLICT_POLICY_REPLACE_EXISTING,
    }
)


@dataclass(frozen=True, slots=True)
class DatasetRuvReadinessConfig:
    """RUV-readiness reporting configuration for future correction support."""

    enabled: bool = False
    control_feature_column: str = "is_control_feature"
    replicate_group_column: str = "replicate_group"
    batch_column: str | None = "batch"

    def __post_init__(self) -> None:
        enabled = self.enabled
        if not isinstance(enabled, bool):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.ruv_readiness.enabled "
                "must be a bool"
            )
        control_feature_column = self.control_feature_column
        if (
            not isinstance(control_feature_column, str)
            or not control_feature_column.strip()
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.ruv_readiness."
                "control_feature_column must be a non-empty string"
            )
        replicate_group_column = self.replicate_group_column
        if (
            not isinstance(replicate_group_column, str)
            or not replicate_group_column.strip()
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.ruv_readiness."
                "replicate_group_column must be a non-empty string"
            )
        batch_column = self.batch_column
        if batch_column is None:
            return
        if not isinstance(batch_column, str) or not batch_column.strip():
            raise PhosPyInputError(
                "dataset build request preprocessing_config.ruv_readiness."
                "batch_column must be None or a non-empty string"
            )


@dataclass(frozen=True, slots=True)
class DatasetSiteSequenceResolutionConfig:
    """Optional local-FASTA-backed site-sequence resolution policy."""

    fasta_path: str | None = None
    mode: DatasetSiteSequenceResolutionMode = (
        DATASET_SITE_SEQUENCE_RESOLUTION_MODE_VALIDATE_EXISTING_AND_FILL_MISSING
    )
    conflict_policy: DatasetSiteSequenceConflictPolicy | None = None
    flank_size: int = 7
    accession_column: str = "protein_accession"
    site_column: str = "site"

    def __post_init__(self) -> None:
        validate_site_sequence_resolution_config(
            mode=self.mode,
            conflict_policy=self.conflict_policy,
            flank_size=self.flank_size,
            accession_column=self.accession_column,
            site_column=self.site_column,
            fasta_path=self.fasta_path,
            supported_modes=DATASET_SITE_SEQUENCE_RESOLUTION_MODES,
            supported_conflict_policies=DATASET_SITE_SEQUENCE_CONFLICT_POLICIES,
        )


__all__ = [
    "DATASET_SITE_SEQUENCE_CONFLICT_POLICIES",
    "DATASET_SITE_SEQUENCE_CONFLICT_POLICY_ERROR",
    "DATASET_SITE_SEQUENCE_CONFLICT_POLICY_PRESERVE_EXISTING",
    "DATASET_SITE_SEQUENCE_CONFLICT_POLICY_REPLACE_EXISTING",
    "DATASET_SITE_SEQUENCE_RESOLUTION_MODES",
    "DATASET_SITE_SEQUENCE_RESOLUTION_MODE_FILL_MISSING_ONLY",
    "DATASET_SITE_SEQUENCE_RESOLUTION_MODE_REPLACE_EXISTING",
    "DATASET_SITE_SEQUENCE_RESOLUTION_MODE_VALIDATE_EXISTING_AND_FILL_MISSING",
    "DATASET_SITE_SEQUENCE_RESOLUTION_MODE_VALIDATE_EXISTING_ONLY",
    "DatasetRuvReadinessConfig",
    "DatasetSiteSequenceConflictPolicy",
    "DatasetSiteSequenceResolutionConfig",
    "DatasetSiteSequenceResolutionMode",
]
