"""Total-protein correction identity-policy resolution."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.provenance.hashing import hash_table_tolerance
from phospy.science.configs.preprocessing import (
    DATASET_TOTAL_PROTEIN_CORRECTION_DUPLICATE_POLICIES,
    DATASET_TOTAL_PROTEIN_CORRECTION_DUPLICATE_POLICY_ERROR,
    DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_DIRECT,
    DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_MAPPING_TABLE,
    DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODES,
    DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICIES,
    DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ERROR,
    DatasetTotalProteinCorrectionDuplicatePolicy,
    DatasetTotalProteinCorrectionIdentityConfig,
    DatasetTotalProteinCorrectionIdentityMode,
    DatasetTotalProteinCorrectionUnmatchedPolicy,
)
from phospy.science.datasets.preprocessing.policy_models import (
    TotalProteinCorrectionIdentityMatchingPolicy,
)


@dataclass(frozen=True, slots=True)
class TotalProteinCorrectionIdentityPolicy:
    """Resolved identity policy consumed by total/protein correction stages."""

    mode: DatasetTotalProteinCorrectionIdentityMode
    phosphosite_key: str
    total_protein_key: str
    matching_policy: TotalProteinCorrectionIdentityMatchingPolicy
    duplicate_policy: DatasetTotalProteinCorrectionDuplicatePolicy
    unmatched_policy: DatasetTotalProteinCorrectionUnmatchedPolicy
    mapping_table: tuple[tuple[str, str], ...] | None = None
    mapping_phosphosite_key: str | None = None
    mapping_total_protein_key: str | None = None
    mapping_table_fingerprint: str | None = None

    def __post_init__(self) -> None:
        mode = str(self.mode).strip()
        if mode not in DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODES:
            raise PhosPyInputError(
                "dataset preprocessing plan total_protein_correction "
                "identity.mode (internal model) contains an unsupported mode"
            )
        object.__setattr__(self, "mode", mode)
        object.__setattr__(
            self,
            "phosphosite_key",
            _require_non_empty_string(
                self.phosphosite_key,
                field_name=(
                    "dataset preprocessing plan total_protein_correction "
                    "identity.phosphosite_key (internal model)"
                ),
            ),
        )
        object.__setattr__(
            self,
            "total_protein_key",
            _require_non_empty_string(
                self.total_protein_key,
                field_name=(
                    "dataset preprocessing plan total_protein_correction "
                    "identity.total_protein_key (internal model)"
                ),
            ),
        )
        object.__setattr__(
            self,
            "matching_policy",
            TotalProteinCorrectionIdentityMatchingPolicy.parse(
                self.matching_policy,
                field_name=(
                    "dataset preprocessing plan total_protein_correction "
                    "identity.matching_policy (internal model)"
                ),
            ),
        )
        duplicate_policy = str(self.duplicate_policy).strip()
        if duplicate_policy not in DATASET_TOTAL_PROTEIN_CORRECTION_DUPLICATE_POLICIES:
            raise PhosPyInputError(
                "dataset preprocessing plan total_protein_correction "
                "identity.duplicate_policy (internal model) contains an "
                "unsupported policy"
            )
        object.__setattr__(self, "duplicate_policy", duplicate_policy)
        unmatched_policy = str(self.unmatched_policy).strip()
        if unmatched_policy not in DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICIES:
            raise PhosPyInputError(
                "dataset preprocessing plan total_protein_correction "
                "identity.unmatched_policy (internal model) contains an "
                "unsupported policy"
            )
        object.__setattr__(self, "unmatched_policy", unmatched_policy)
        if mode == DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_DIRECT:
            _reject_mapping_fields_for_direct_mode(self)
            return
        _validate_mapping_table_identity(self)


class TotalProteinCorrectionIdentityResolver:
    """Resolve public total-protein identity config into an internal policy."""

    def run(
        self,
        config: DatasetTotalProteinCorrectionIdentityConfig,
    ) -> TotalProteinCorrectionIdentityPolicy:
        if not isinstance(config, DatasetTotalProteinCorrectionIdentityConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                "identity must be a DatasetTotalProteinCorrectionIdentityConfig"
            )
        if config.mode == DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_DIRECT:
            return self._resolve_direct(config)
        if config.mode == DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_MAPPING_TABLE:
            return self._resolve_mapping_table(config)
        raise PhosPyInputError(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity contains an unsupported mode"
        )

    def _resolve_direct(
        self,
        config: DatasetTotalProteinCorrectionIdentityConfig,
    ) -> TotalProteinCorrectionIdentityPolicy:
        return TotalProteinCorrectionIdentityPolicy(
            mode=config.mode,
            phosphosite_key=str(config.phosphosite_key).strip(),
            total_protein_key=str(config.total_protein_key).strip(),
            matching_policy=TotalProteinCorrectionIdentityMatchingPolicy.parse(
                config.matching_policy,
                field_name=(
                    "preprocessing_config.total_protein_correction.identity."
                    "matching_policy"
                ),
            ),
            duplicate_policy=config.duplicate_policy,
            unmatched_policy=config.unmatched_policy,
            mapping_table=None,
            mapping_phosphosite_key=None,
            mapping_total_protein_key=None,
            mapping_table_fingerprint=None,
        )

    def _resolve_mapping_table(
        self,
        config: DatasetTotalProteinCorrectionIdentityConfig,
    ) -> TotalProteinCorrectionIdentityPolicy:
        mapping_table = config.mapping_table
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
        mapping_phosphosite_key = str(config.mapping_phosphosite_key).strip()
        mapping_total_protein_key = str(config.mapping_total_protein_key).strip()
        _validate_mapping_key_exists(
            mapping_table=mapping_table,
            key=mapping_phosphosite_key,
        )
        _validate_mapping_key_exists(
            mapping_table=mapping_table,
            key=mapping_total_protein_key,
        )
        normalized_table = _normalize_mapping_table(
            mapping_table=mapping_table,
            mapping_phosphosite_key=mapping_phosphosite_key,
            mapping_total_protein_key=mapping_total_protein_key,
        )
        return TotalProteinCorrectionIdentityPolicy(
            mode=config.mode,
            phosphosite_key=str(config.phosphosite_key).strip(),
            total_protein_key=str(config.total_protein_key).strip(),
            matching_policy=TotalProteinCorrectionIdentityMatchingPolicy.parse(
                config.matching_policy,
                field_name=(
                    "preprocessing_config.total_protein_correction.identity."
                    "matching_policy"
                ),
            ),
            duplicate_policy=config.duplicate_policy,
            unmatched_policy=config.unmatched_policy,
            mapping_table=_mapping_rows(normalized_table),
            mapping_phosphosite_key=mapping_phosphosite_key,
            mapping_total_protein_key=mapping_total_protein_key,
            mapping_table_fingerprint=_fingerprint_mapping_table(normalized_table),
        )


def _require_non_empty_string(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return value.strip()


def _reject_mapping_fields_for_direct_mode(
    policy: TotalProteinCorrectionIdentityPolicy,
) -> None:
    if policy.mapping_table is not None:
        raise PhosPyInputError(
            "dataset preprocessing plan total_protein_correction identity.mapping_table "
            "(internal model) must be None when identity.mode='direct'"
        )
    if policy.mapping_phosphosite_key is not None:
        raise PhosPyInputError(
            "dataset preprocessing plan total_protein_correction "
            "identity.mapping_phosphosite_key (internal model) must be None "
            "when identity.mode='direct'"
        )
    if policy.mapping_total_protein_key is not None:
        raise PhosPyInputError(
            "dataset preprocessing plan total_protein_correction "
            "identity.mapping_total_protein_key (internal model) must be None "
            "when identity.mode='direct'"
        )
    if policy.mapping_table_fingerprint is not None:
        raise PhosPyInputError(
            "dataset preprocessing plan total_protein_correction "
            "identity.mapping_table_fingerprint (internal model) must be None "
            "when identity.mode='direct'"
        )


def _validate_mapping_table_identity(
    policy: TotalProteinCorrectionIdentityPolicy,
) -> None:
    if policy.mapping_table is None:
        raise PhosPyInputError(
            "dataset preprocessing plan total_protein_correction identity.mapping_table "
            "(internal model) is required when identity.mode='mapping_table'"
        )
    if not isinstance(policy.mapping_table, tuple):
        raise PhosPyInputError(
            "dataset preprocessing plan total_protein_correction identity.mapping_table "
            "(internal model) must be an immutable tuple"
        )
    for row in policy.mapping_table:
        if (
            not isinstance(row, tuple)
            or len(row) != 2
            or not all(isinstance(value, str) for value in row)
        ):
            raise PhosPyInputError(
                "dataset preprocessing plan total_protein_correction "
                "identity.mapping_table (internal model) must contain "
                "(phosphosite_id, total_protein_id) string pairs"
            )
    mapping_phosphosite_key = _require_non_empty_string(
        policy.mapping_phosphosite_key,
        field_name=(
            "dataset preprocessing plan total_protein_correction "
            "identity.mapping_phosphosite_key (internal model)"
        ),
    )
    mapping_total_protein_key = _require_non_empty_string(
        policy.mapping_total_protein_key,
        field_name=(
            "dataset preprocessing plan total_protein_correction "
            "identity.mapping_total_protein_key (internal model)"
        ),
    )
    object.__setattr__(
        policy,
        "mapping_phosphosite_key",
        mapping_phosphosite_key,
    )
    object.__setattr__(
        policy,
        "mapping_total_protein_key",
        mapping_total_protein_key,
    )
    _require_non_empty_string(
        policy.mapping_table_fingerprint,
        field_name=(
            "dataset preprocessing plan total_protein_correction "
            "identity.mapping_table_fingerprint (internal model)"
        ),
    )


def _validate_mapping_key_exists(
    *,
    mapping_table: pd.DataFrame,
    key: str,
) -> None:
    if key in mapping_table.columns:
        return
    raise PhosPyInputError(
        "dataset build request preprocessing_config.total_protein_correction."
        "identity.mapping_table is missing column "
        f"{key!r}"
    )


def _normalize_mapping_table(
    *,
    mapping_table: pd.DataFrame,
    mapping_phosphosite_key: str,
    mapping_total_protein_key: str,
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "phosphosite_id": mapping_table.loc[:, mapping_phosphosite_key]
            .astype("string")
            .str.strip(),
            "total_protein_id": mapping_table.loc[:, mapping_total_protein_key]
            .astype("string")
            .str.strip(),
        }
    )


def _is_missing_mapping_value(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _mapping_rows(
    normalized_table: pd.DataFrame,
) -> tuple[tuple[str, str], ...]:
    return tuple(
        (
            ""
            if _is_missing_mapping_value(record.get("phosphosite_id"))
            else str(record.get("phosphosite_id")),
            ""
            if _is_missing_mapping_value(record.get("total_protein_id"))
            else str(record.get("total_protein_id")),
        )
        for record in normalized_table.to_dict(orient="records")
    )


def _fingerprint_mapping_table(normalized_table: pd.DataFrame) -> str:
    fingerprint_table = (
        normalized_table.fillna("<MISSING>")
        .sort_values(by=["phosphosite_id", "total_protein_id"], kind="mergesort")
        .reset_index(drop=True)
    )
    return hash_table_tolerance(
        fingerprint_table,
        name="total_protein_correction.identity.mapping_table",
    )


DEFAULT_TOTAL_PROTEIN_CORRECTION_IDENTITY_POLICY = TotalProteinCorrectionIdentityPolicy(
    mode=DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_DIRECT,
    phosphosite_key="gene_symbol",
    total_protein_key="__index__",
    matching_policy=TotalProteinCorrectionIdentityMatchingPolicy.STRICT,
    duplicate_policy=DATASET_TOTAL_PROTEIN_CORRECTION_DUPLICATE_POLICY_ERROR,
    unmatched_policy=DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ERROR,
    mapping_table=None,
    mapping_phosphosite_key=None,
    mapping_total_protein_key=None,
    mapping_table_fingerprint=None,
)


__all__ = [
    "DEFAULT_TOTAL_PROTEIN_CORRECTION_IDENTITY_POLICY",
    "TotalProteinCorrectionIdentityPolicy",
    "TotalProteinCorrectionIdentityResolver",
]
