"""Models for peptide-to-site differential aggregation."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

import pandas as pd

from phospy._frame_ownership import export_dataframe, own_dataframe
from phospy.errors.input import PhosPyInputError
from phospy.scientific_policies import ScientificPolicyRecord

PEPTIDE_TO_SITE_STRATEGY_COMPAT_BEST_P_VALUE = "compat_best_p_value"
PEPTIDE_TO_SITE_STRATEGY_INVERSE_VARIANCE_WEIGHTED = "inverse_variance_weighted"
PEPTIDE_TO_SITE_STRATEGY_FIXED_EFFECT_META = "fixed_effect_meta"
PEPTIDE_TO_SITE_STRATEGY_STOUFFER_Z = "stouffer_z"
PEPTIDE_TO_SITE_STRATEGY_RANDOM_EFFECT_META = "random_effect_meta"
SUPPORTED_PEPTIDE_TO_SITE_STRATEGIES: tuple[str, ...] = (
    PEPTIDE_TO_SITE_STRATEGY_COMPAT_BEST_P_VALUE,
    PEPTIDE_TO_SITE_STRATEGY_INVERSE_VARIANCE_WEIGHTED,
    PEPTIDE_TO_SITE_STRATEGY_FIXED_EFFECT_META,
    PEPTIDE_TO_SITE_STRATEGY_STOUFFER_Z,
    PEPTIDE_TO_SITE_STRATEGY_RANDOM_EFFECT_META,
)
STOUFFER_WEIGHTING_EQUAL = "equal"
STOUFFER_WEIGHTING_INVERSE_VARIANCE = "inverse_variance"
SUPPORTED_STOUFFER_WEIGHTING: tuple[str, ...] = (
    STOUFFER_WEIGHTING_EQUAL,
    STOUFFER_WEIGHTING_INVERSE_VARIANCE,
)
MISSING_VARIANCE_POLICY_DROP = "drop"
SUPPORTED_MISSING_VARIANCE_POLICIES: tuple[str, ...] = (MISSING_VARIANCE_POLICY_DROP,)


@dataclass(frozen=True, slots=True)
class PeptideToSiteAggregationConfig:
    """Configuration for peptide-to-site differential aggregation."""

    strategy: str = PEPTIDE_TO_SITE_STRATEGY_FIXED_EFFECT_META
    min_peptides_per_site: int = 1
    missing_variance_policy: str = MISSING_VARIANCE_POLICY_DROP
    stouffer_weighting: str = STOUFFER_WEIGHTING_EQUAL
    random_effect_tau2_floor: float = 0.0

    def __post_init__(self) -> None:
        if self.strategy not in SUPPORTED_PEPTIDE_TO_SITE_STRATEGIES:
            supported = ", ".join(
                repr(value) for value in SUPPORTED_PEPTIDE_TO_SITE_STRATEGIES
            )
            raise PhosPyInputError(
                f"peptide_to_site_aggregation.strategy must be one of: {supported}"
            )
        if (
            not isinstance(self.min_peptides_per_site, int)
            or self.min_peptides_per_site < 1
        ):
            raise PhosPyInputError(
                "peptide_to_site_aggregation.min_peptides_per_site must be an int >= 1"
            )
        if self.missing_variance_policy not in SUPPORTED_MISSING_VARIANCE_POLICIES:
            supported = ", ".join(
                repr(value) for value in SUPPORTED_MISSING_VARIANCE_POLICIES
            )
            raise PhosPyInputError(
                "peptide_to_site_aggregation.missing_variance_policy must be one of: "
                f"{supported}"
            )
        if self.stouffer_weighting not in SUPPORTED_STOUFFER_WEIGHTING:
            supported = ", ".join(repr(value) for value in SUPPORTED_STOUFFER_WEIGHTING)
            raise PhosPyInputError(
                f"peptide_to_site_aggregation.stouffer_weighting must be one of: {supported}"
            )
        if (
            isinstance(self.random_effect_tau2_floor, bool)
            or not isinstance(self.random_effect_tau2_floor, int | float)
            or float(self.random_effect_tau2_floor) < 0.0
        ):
            raise PhosPyInputError(
                "peptide_to_site_aggregation.random_effect_tau2_floor must be a "
                "numeric value >= 0.0"
            )
        object.__setattr__(
            self, "random_effect_tau2_floor", float(self.random_effect_tau2_floor)
        )


@dataclass(frozen=True, slots=True, init=False)
class PeptideToSiteAggregationResult:
    """Site-level differential result from peptide-level aggregation."""

    contrast_name: str
    table: pd.DataFrame
    warnings: tuple[str, ...]
    provenance: Mapping[str, object]
    scientific_policies: tuple[ScientificPolicyRecord, ...]

    def __init__(
        self,
        *,
        contrast_name: str,
        table: pd.DataFrame,
        warnings: tuple[str, ...] = (),
        provenance: Mapping[str, object] | None = None,
        scientific_policies: tuple[ScientificPolicyRecord, ...] = (),
        _assume_owned: bool = False,
    ) -> None:
        if not isinstance(contrast_name, str) or not contrast_name.strip():
            raise PhosPyInputError(
                "peptide_to_site_aggregation_result.contrast_name must be a non-empty string"
            )
        table = own_dataframe(
            table,
            field_name="peptide_to_site_aggregation_result.table",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        required_columns = (
            "logFC",
            "uncertainty_statistic",
            "P.Value",
            "adj.P.Val",
            "n_peptide_observations",
            "n_peptides_used",
        )
        missing = [column for column in required_columns if column not in table.columns]
        if missing:
            joined = ", ".join(missing)
            raise PhosPyInputError(
                "peptide_to_site_aggregation_result.table is missing required columns: "
                f"{joined}"
            )
        warnings_tuple = tuple(str(value) for value in warnings)
        for policy in scientific_policies:
            if not isinstance(policy, ScientificPolicyRecord):
                raise PhosPyInputError(
                    "peptide_to_site_aggregation_result.scientific_policies must "
                    "contain ScientificPolicyRecord values"
                )
        object.__setattr__(self, "contrast_name", contrast_name.strip())
        object.__setattr__(self, "table", table)
        object.__setattr__(self, "warnings", warnings_tuple)
        object.__setattr__(
            self,
            "provenance",
            MappingProxyType(
                {}
                if provenance is None
                else {str(key): value for key, value in provenance.items()}
            ),
        )
        object.__setattr__(self, "scientific_policies", tuple(scientific_policies))

    def to_dataframe(self) -> pd.DataFrame:
        return export_dataframe(self.table)
