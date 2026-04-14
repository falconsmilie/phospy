from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd

from .internal.types import (
    KINASE_PROFILE_MISSING_VALUE_STRATEGY_PROPAGATE_ANY_MISSING,
    KinaseProfileMissingValueStrategy,
)
from .validation.values.enums import validate_missing_value_strategy
from .validation.values.numeric import validate_positive_int

__all__ = [
    "DEFAULT_KINASE_PROFILE_POLICY",
    "KinaseProfilePolicy",
    "KinaseProfileResult",
    "build_kinase_substrate_profiles",
]


@dataclass(frozen=True, slots=True)
class KinaseProfilePolicy:
    """Explicit policy for aggregating multiple quantified kinase substrates.

    ``missing_value_strategy`` controls how column-wise aggregation handles
    missing phosphosite values when more than one quantified substrate is
    available for a kinase:

    - ``"propagate_any_missing"`` keeps the current strict behaviour and returns
      ``NaN`` for a sample whenever any contributing substrate is missing
    - ``"median_skipna"`` ignores missing substrate values when computing the
      column-wise median
    """

    missing_value_strategy: KinaseProfileMissingValueStrategy = (
        KINASE_PROFILE_MISSING_VALUE_STRATEGY_PROPAGATE_ANY_MISSING
    )

    def __post_init__(self) -> None:
        validate_missing_value_strategy(self.missing_value_strategy)

    @classmethod
    def from_value(cls, value: object) -> KinaseProfilePolicy:
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if isinstance(value, Mapping):
            return cls(**dict(value))
        msg = "profile_policy must be a KinaseProfilePolicy or mapping"
        raise TypeError(msg)


DEFAULT_KINASE_PROFILE_POLICY = KinaseProfilePolicy()


@dataclass(slots=True)
class KinaseProfileResult:
    """Kinase substrate profiles built from one phosphosite matrix."""

    profile_matrix: pd.DataFrame
    substrate_counts: pd.Series
    quantified_substrates: dict[str, list[str]]


def build_kinase_substrate_profiles(
    substrate_map: Mapping[str, Sequence[str]],
    phospho_matrix: pd.DataFrame,
    min_substrates: int = 1,
    *,
    policy: KinaseProfilePolicy | None = None,
    missing_value_strategy: KinaseProfileMissingValueStrategy | None = None,
) -> KinaseProfileResult:
    """Build kinase substrate profiles from quantified phosphosite values.

    This mirrors the core behaviour of PhosR's ``kinaseSubstrateProfile()``:
    for each kinase, intersect its known substrates with the phosphosite matrix,
    use the single quantified row directly when exactly one site is available,
    and otherwise summarise quantified substrates column-wise with the median.

    The aggregation policy is explicit so callers can review and override the
    missing-value behaviour rather than relying on a silent hardcoded default.
    """

    validate_positive_int(min_substrates, name="min_substrates")
    resolved_policy = _resolve_profile_policy(
        policy=policy,
        missing_value_strategy=missing_value_strategy,
    )

    observed_sites = set(phospho_matrix.index)
    numeric_matrix = phospho_matrix.astype(float)

    profile_rows: dict[str, pd.Series] = {}
    quantified_substrates: dict[str, list[str]] = {}
    substrate_counts: dict[str, int] = {}

    for kinase, substrates in substrate_map.items():
        substrate_sequence = list(substrates)
        quantified_sites = _quantified_sites(
            substrates=substrate_sequence,
            observed_sites=observed_sites,
        )
        substrate_counts[kinase] = len(quantified_sites)
        if len(quantified_sites) < min_substrates:
            continue

        profile_rows[kinase] = _aggregate_quantified_sites(
            numeric_matrix.loc[quantified_sites, :],
            policy=resolved_policy,
        )
        quantified_substrates[kinase] = quantified_sites

    if profile_rows:
        profile_matrix = pd.DataFrame.from_dict(profile_rows, orient="index")
        profile_matrix = profile_matrix.loc[:, phospho_matrix.columns.copy()]
    else:
        profile_matrix = pd.DataFrame(
            columns=phospho_matrix.columns.copy(), dtype=float
        )

    profile_matrix.index.name = "kinase"

    count_series = pd.Series(substrate_counts, dtype="int64", name="NumSub")
    count_series.index.name = "kinase"

    return KinaseProfileResult(
        profile_matrix=profile_matrix,
        substrate_counts=count_series,
        quantified_substrates=quantified_substrates,
    )


def _resolve_profile_policy(
    *,
    policy: KinaseProfilePolicy | None,
    missing_value_strategy: KinaseProfileMissingValueStrategy | None,
) -> KinaseProfilePolicy:
    if policy is not None and missing_value_strategy is not None:
        msg = (
            "build_kinase_substrate_profiles() accepts either policy or "
            "missing_value_strategy, not both"
        )
        raise ValueError(msg)
    if policy is not None:
        return KinaseProfilePolicy.from_value(policy)
    if missing_value_strategy is not None:
        return KinaseProfilePolicy(missing_value_strategy=missing_value_strategy)
    return DEFAULT_KINASE_PROFILE_POLICY


def _aggregate_quantified_sites(
    quantified_matrix: pd.DataFrame,
    *,
    policy: KinaseProfilePolicy,
) -> pd.Series:
    if quantified_matrix.shape[0] == 1:
        return quantified_matrix.iloc[0].astype(float)

    if (
        policy.missing_value_strategy
        == KINASE_PROFILE_MISSING_VALUE_STRATEGY_PROPAGATE_ANY_MISSING
    ):
        return quantified_matrix.median(axis=0, skipna=False).astype(float)

    return quantified_matrix.median(axis=0, skipna=True).astype(float)


def _quantified_sites(
    substrates: Sequence[str],
    observed_sites: set[str],
) -> list[str]:
    quantified: list[str] = []
    seen: set[str] = set()

    for substrate in substrates:
        if substrate in observed_sites and substrate not in seen:
            quantified.append(substrate)
            seen.add(substrate)

    return quantified
